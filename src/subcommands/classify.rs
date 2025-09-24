use anyhow::Result;
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;

use crate::io::{bam, fasta, paf, runfiles};
use crate::model::{ClassifyParams, Weights};
use crate::scoring;
use crate::util::{logging, mapping};
use crate::{model, summary};

#[derive(Args, Debug)]
pub struct CmdClassify {
    #[arg(long)]
    pub mito: PathBuf,
    #[arg(long)]
    pub nuclear: PathBuf,
    #[arg(
        long,
        value_delimiter = ',',
        num_args = 1..,           // ← at least one value required
        help = "One or more reads files (FASTQ/FASTA; .gz ok). Repeat or comma-separate."
    )]
    pub reads: Vec<PathBuf>,
    #[arg(long, value_parser=["hifi","ont"])]
    pub platform: String,
    #[arg(long)]
    pub out: PathBuf,

    #[arg(long, help = "Path to minimap2 (else PATH)")]
    pub minimap2: Option<PathBuf>,
    #[arg(long, help = "Path to samtools (else PATH)")]
    pub samtools: Option<PathBuf>,
    #[arg(long, help = "Threads (default: logical CPUs, capped at 16)")]
    pub threads: Option<usize>,
    #[arg(long, help = "Keep tmp/ outputs so they can be reused")]
    pub keep_tmp: bool,
}

impl CmdClassify {
    pub fn run(self) -> Result<()> {
        // 0) Preflight
        fs::create_dir_all(&self.out)?;
        logging::init_logging(&self.out)?;
        log::info!("onsm classify started");

        fasta::validate_fasta(&self.mito)?;
        fasta::validate_fasta(&self.nuclear)?;

        for r in &self.reads {
            runfiles::ensure_exists(r)?;
        }

        // Resolve binaries once
        let (mm2_bin, sam_bin) =
            mapping::resolve_bins(self.minimap2.as_deref(), self.samtools.as_deref())?;
        log::info!("Using minimap2 at {}", mm2_bin.display());
        log::info!("Using samtools at {}", sam_bin.display());
        if let Ok(v) = mapping::get_version(&mm2_bin) {
            log::info!("minimap2: {v}");
        }
        if let Ok(v) = mapping::get_version(&sam_bin) {
            log::info!("samtools: {v}");
        }

        let threads = self.threads.unwrap_or_else(|| {
            let n = std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(4);
            n.min(16)
        });
        log::info!("Threads: {threads}");

        let manifest = model::RunManifest::new(
            &self.mito,
            &self.nuclear,
            &self.reads,
            &self.platform,
            threads,
            model::MIN_ID,
            model::MIN_LEN,
            model::MERGE_GAP,
            model::FLANK_BP,
            model::WIN_BP,
        );
        model::RunManifest::save_to(&self.out, &manifest)?;

        // 1) Asm↔Asm → PAF
        let tmp = self.out.join("tmp");
        fs::create_dir_all(&tmp)?;
        let paf_m2n = tmp.join("mito_to_nuc.paf");
        let paf_n2m = tmp.join("nuc_to_mito.paf");
        mapping::map_asm_to_asm(&mm2_bin, &self.mito, &self.nuclear, &paf_m2n, threads)?;
        mapping::map_asm_to_asm(&mm2_bin, &self.nuclear, &self.mito, &paf_n2m, threads)?;

        // 2) reads→ref → BAM
        let bam_r2n = tmp.join("reads_to_nuc.bam");
        let bam_r2m = tmp.join("reads_to_mito.bam");
        mapping::map_reads_to_ref(
            &mm2_bin,
            &sam_bin,
            &self.platform,
            &self.reads,
            &self.nuclear,
            &bam_r2n,
            threads,
        )?;
        mapping::map_reads_to_ref(
            &mm2_bin,
            &sam_bin,
            &self.platform,
            &self.reads,
            &self.mito,
            &bam_r2m,
            threads,
        )?;

        // 3) Parse PAF + pair
        let m2n = paf::read_paf(&paf_m2n, model::MIN_ID, model::MIN_LEN)?;
        let n2m = paf::read_paf(&paf_n2m, model::MIN_ID, model::MIN_LEN)?;
        let pairs = paf::pair_and_merge(&m2n, n2m, model::MERGE_GAP)?;
        log::info!("paired {} candidate loci", pairs.len());

        // 4) Coverage & spans (samtools)
        let (coverage, spans) = bam::compute_coverage_and_spans_with_tools(
            &bam_r2n,
            &bam_r2m,
            &pairs,
            model::FLANK_BP,
            model::WIN_BP,
            &sam_bin,
        )?;

        // 5) Score & classify (fixed params)
        let weights = Weights::default();
        let params = ClassifyParams::default();
        let (pairs_tsv, classes_tsv) =
            scoring::classify_pairs(&pairs, &coverage, &spans, weights, params)?;

        // 6) Write outputs
        fs::write(self.out.join("pairs.tsv"), pairs_tsv)?;
        fs::write(self.out.join("classification.tsv"), classes_tsv.clone())?;
        serde_json::to_writer_pretty(fs::File::create(self.out.join("coverage.json"))?, &coverage)?;

        let calls = summary::parse_calls_tsv_str(&classes_tsv);

        let summary_tbl = summary::compute_percentages(&self.mito, &self.nuclear, &pairs, &calls)?;
        summary::write_summary_tsv(&self.out.join("summary.tsv"), &summary_tbl)?;

        // 7) Cleanup
        if !self.keep_tmp {
            let _ = fs::remove_dir_all(&tmp);
        } else {
            log::info!("Keeping tmp/ for reuse: {}", tmp.display());
        }

        log::info!("done.");
        Ok(())
    }
}
