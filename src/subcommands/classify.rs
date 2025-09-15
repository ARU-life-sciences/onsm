//! `onsm classify` — main pipeline (HiFi/ONT long reads only).

use anyhow::{anyhow, Result};
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;

use crate::io::{bam, fasta, paf, runfiles};
use crate::model::{ClassifyParams, RunManifest, Weights};
use crate::scoring;
use crate::subcommands::reuse::parse_calls_tsv;
use crate::util::{logging, mapping};

#[derive(Args, Debug)]
pub struct CmdClassify {
    #[arg(long)]
    pub mito: PathBuf,
    #[arg(long)]
    pub nuclear: PathBuf,
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,
    #[arg(long, value_parser=["hifi","ont"])]
    pub platform: String,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(
        long,
        value_name = "PATH",
        help = "Path to minimap2 binary (else use PATH)"
    )]
    pub minimap2: Option<PathBuf>,
    #[arg(
        long,
        value_name = "PATH",
        help = "Path to samtools binary (else use PATH)"
    )]
    pub samtools: Option<PathBuf>,
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    #[arg(long, default_value_t = 0.90)]
    pub min_id: f32,
    #[arg(long, default_value_t = 100)]
    pub min_len: u32,
    #[arg(long, default_value_t = 50)]
    pub merge_gap: u32,
    #[arg(long, default_value_t = 500)]
    pub flank: u32,
    #[arg(long, default_value_t = 250)]
    pub win: u32,
    #[arg(long, default_value_t = 1)]
    pub seed: u64,
    #[arg(long)]
    pub keep_tmp: bool,

    // scoring weights
    #[arg(long, default_value_t = 0.25)]
    pub w_a: f32,
    #[arg(long, default_value_t = 0.15)]
    pub w_l: f32,
    #[arg(long, default_value_t = 0.25)]
    pub w_d: f32,
    #[arg(long, default_value_t = 0.25)]
    pub w_s: f32,
    #[arg(long, default_value_t = 0.10)]
    pub w_f: f32,
    #[arg(long, default_value_t = 0.15)]
    pub call_threshold: f32,
    #[arg(long, default_value_t = 0.30)]
    pub highconf_threshold: f32,
}

impl CmdClassify {
    pub fn run(self) -> Result<()> {
        // 0) Preflight
        fs::create_dir_all(&self.out)?;
        logging::init_logging(&self.out)?;
        log::info!("onsm classify started");

        fasta::validate_fasta(&self.mito)?;
        fasta::validate_fasta(&self.nuclear)?;
        if self.reads.is_empty() {
            return Err(anyhow!("--reads: at least one FASTQ(.gz) required"));
        }
        for r in &self.reads {
            runfiles::ensure_exists(r)?;
        }

        let manifest = RunManifest::from_inputs(
            &self.mito,
            &self.nuclear,
            &self.reads,
            &self.platform,
            self.threads,
            self.min_id,
            self.min_len,
            self.merge_gap,
            self.flank,
            self.win,
            self.seed,
        )?;

        serde_json::to_writer_pretty(
            fs::File::create(self.out.join("run_manifest.json"))?,
            &manifest,
        )?;

        // Resolve binaries once from flags/env/PATH
        let (mm2_bin, sam_bin) =
            mapping::resolve_bins(self.minimap2.as_deref(), self.samtools.as_deref())?;
        log::info!("Using minimap2 at {}", mm2_bin.display());
        log::info!("Using samtools at {}", sam_bin.display());
        // Optionally:
        if let Ok(v) = mapping::get_version(&mm2_bin) {
            log::info!("minimap2: {v}");
        }
        if let Ok(v) = mapping::get_version(&sam_bin) {
            log::info!("samtools: {v}");
        }

        // 1) Assembly-to-assembly mappings → PAF (mito→nuclear, nuclear→mito)
        let tmp_dir = self.out.join("tmp");
        fs::create_dir_all(&tmp_dir)?;
        let paf_m2n = tmp_dir.join("mito_to_nuc.paf");
        let paf_n2m = tmp_dir.join("nuc_to_mito.paf");

        mapping::map_asm_to_asm(&mm2_bin, &self.mito, &self.nuclear, &paf_m2n, self.threads)?;
        mapping::map_asm_to_asm(&mm2_bin, &self.nuclear, &self.mito, &paf_n2m, self.threads)?;

        // 2) Reads-to-reference mappings → BAM (reads→nuclear, reads→mito)
        let bam_r2n = tmp_dir.join("reads_to_nuc.bam");
        let bam_r2m = tmp_dir.join("reads_to_mito.bam");
        mapping::map_reads_to_ref(
            &mm2_bin,
            &sam_bin,
            &self.platform,
            &self.reads,
            &self.nuclear,
            &bam_r2n,
            self.threads,
        )?;
        mapping::map_reads_to_ref(
            &mm2_bin,
            &sam_bin,
            &self.platform,
            &self.reads,
            &self.mito,
            &bam_r2m,
            self.threads,
        )?;

        // 3) Parse PAF, merge HSPs into loci
        let m2n_records = paf::read_paf(&paf_m2n, self.min_id, self.min_len)?;
        let n2m_records = paf::read_paf(&paf_n2m, self.min_id, self.min_len)?;
        let pairs = paf::pair_and_merge(m2n_records, n2m_records, self.merge_gap)?;

        // 4) Coverage & spans
        let (coverage, spans) = bam::compute_coverage_and_spans_with_tools(
            &bam_r2n, &bam_r2m, &pairs, self.flank, self.win, &sam_bin,
        )?;

        // 5) Flank embedding: identity via short asm-to-asm alignments (stubbed; fill later)
        let embed = paf::compute_embedding_identities(
            &pairs,
            &self.mito,
            &self.nuclear,
            self.flank,
            self.threads,
        )?;

        // 6) Score & classify
        let weights = Weights {
            w_a: self.w_a,
            w_l: self.w_l,
            w_d: self.w_d,
            w_s: self.w_s,
            w_f: self.w_f,
        };
        let params = ClassifyParams {
            call_threshold: self.call_threshold,
            highconf_threshold: self.highconf_threshold,
        };
        let (pairs_tsv, classes_tsv) =
            scoring::classify_pairs(&pairs, &coverage, &spans, &embed, weights, params)?;

        // 7) Write outputs
        fs::write(self.out.join("pairs.tsv"), pairs_tsv)?;
        fs::write(self.out.join("classification.tsv"), classes_tsv.clone())?;
        serde_json::to_writer_pretty(fs::File::create(self.out.join("coverage.json"))?, &coverage)?;

        // 8) Cleanup
        if !self.keep_tmp {
            let _ = fs::remove_dir_all(tmp_dir);
        }

        let calls = parse_calls_tsv(&classes_tsv);
        let summary =
            crate::summary::compute_percentages(&self.mito, &self.nuclear, &pairs, &calls)?;
        crate::summary::write_summary_tsv(&self.out.join("summary.tsv"), &summary)?;
        log::info!(
            "SUMMARY: NUMT in nuclear = {:.4}%, NIMT in mito = {:.4}% ; mito covered by NUMT homologs = {:.4}%",
            summary.nuc_pct_numt, summary.mito_pct_nimt, summary.mito_pct_covered_by_numt_homologs
        );

        log::info!("done.");
        Ok(())
    }
}
