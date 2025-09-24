use anyhow::Result;
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;

use crate::io::{bam, paf};
use crate::model::{self, ClassifyParams, Weights};
use crate::scoring;
use crate::summary;
use crate::util::{logging, mapping};

#[derive(Args, Debug)]
pub struct CmdReuse {
    /// Output directory from a previous `onsm classify`
    #[arg(long, value_name = "DIR")]
    pub from: PathBuf,

    /// Where to write re-scored outputs (default: reuse_<timestamp> inside --from)
    #[arg(long, num_args = 1)]
    pub out_dir: PathBuf,

    /// Override samtools (else PATH)
    #[arg(long)]
    pub samtools: Option<PathBuf>,

    /// Optional: override minimap2 for any future embedding features
    #[arg(long)]
    pub minimap2: Option<PathBuf>,
}

impl CmdReuse {
    pub fn run(self) -> Result<()> {
        logging::init_logging(&self.out_dir)?;

        // 1) Load manifest
        let m = model::RunManifest::load_from(&self.from)?;
        let tmp = self.from.join("tmp");

        // 2) Resolve tools (samtools used for coverage)
        let (_mm2_bin, sam_bin) =
            mapping::resolve_bins(self.minimap2.as_deref(), self.samtools.as_deref())?;
        log::info!("REUSE: using samtools at {}", sam_bin.display());

        // 3) Derive artifact paths from the previous run
        let paf_m2n = tmp.join("mito_to_nuc.paf");
        let paf_n2m = tmp.join("nuc_to_mito.paf");
        let bam_r2n = tmp.join("reads_to_nuc.bam");
        let bam_r2m = tmp.join("reads_to_mito.bam");

        for p in [&paf_m2n, &paf_n2m, &bam_r2n, &bam_r2m] {
            if !p.exists() {
                anyhow::bail!("Required artifact missing: {}", p.display());
            }
        }

        // 4) Prepare new out dir
        fs::create_dir_all(&self.out_dir)?;

        // 5) Parse & pair
        let m2n = paf::read_paf(&paf_m2n, m.min_id, m.min_len)?;
        let n2m = paf::read_paf(&paf_n2m, m.min_id, m.min_len)?;
        let pairs = paf::pair_and_merge(&m2n, n2m, m.merge_gap)?;
        log::info!("REUSE: paired {} candidate loci", pairs.len());

        // 6) Coverage & spans
        let (coverage, spans) = bam::compute_coverage_and_spans_with_tools(
            &bam_r2n, &bam_r2m, &pairs, m.flank_bp, m.win_bp, &sam_bin,
        )?;

        // 7) Score & classify (same defaults)
        let weights = Weights::default();
        let params = ClassifyParams::default();
        let (pairs_tsv, classes_tsv) =
            scoring::classify_pairs(&pairs, &coverage, &spans, weights, params)?;

        // 8) Write outputs
        fs::write(self.out_dir.join("pairs.tsv"), &pairs_tsv)?;
        fs::write(self.out_dir.join("classification.tsv"), &classes_tsv)?;
        serde_json::to_writer_pretty(
            fs::File::create(self.out_dir.join("coverage.json"))?,
            &coverage,
        )?;

        // 9) Summary (recomputed on the new outputs)
        let calls = summary::parse_calls_tsv_str(&classes_tsv);
        let summary_tbl = summary::compute_percentages(&m.mito, &m.nuclear, &pairs, &calls)?;
        summary::write_summary_tsv(&self.out_dir.join("summary.tsv"), &summary_tbl)?;

        log::info!("REUSE: done → {}", self.out_dir.display());
        Ok(())
    }
}
