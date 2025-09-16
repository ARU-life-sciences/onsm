use anyhow::Result;
use clap::Args;
use fs_err as fs;
use log::info;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::time::Instant;

use crate::io::{bam, paf};
use crate::model::{ClassifyParams, Weights};
use crate::scoring;
use crate::summary;

pub fn parse_calls_tsv(s: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for line in s.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(pid), Some(call)) = (it.next(), it.next()) {
            m.insert(pid.to_string(), call.to_string());
        }
    }
    m
}

#[derive(Args, Debug)]
pub struct CmdReuse {
    #[arg(long)]
    pub mito: PathBuf,
    #[arg(long)]
    pub nuclear: PathBuf,
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,
    #[arg(long, value_parser=["hifi","ont"])]
    pub platform: String,
    #[arg(long)]
    pub samtools_bin: Option<PathBuf>,
    #[arg(long)]
    pub paf_mito_to_nuc: PathBuf,
    #[arg(long)]
    pub paf_nuc_to_mito: PathBuf,
    #[arg(long)]
    pub bam_reads_to_nuc: PathBuf,
    #[arg(long)]
    pub bam_reads_to_mito: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, default_value_t = 0.15)]
    pub call_threshold: f32,
    #[arg(long)]
    pub skip_hash_check: bool,

    // PAF filters
    #[arg(long, default_value_t = 0.90)]
    pub min_id: f32,
    #[arg(long, default_value_t = 100)]
    pub min_len: u32,

    #[arg(
        long,
        default_value_t = false,
        help = "Keep nuclear contigs that look like full mito assemblies"
    )]
    pub allow_organelle_in_nuclear: bool,
}

impl CmdReuse {
    pub fn run(self) -> Result<()> {
        let t0 = Instant::now();
        info!("REUSE: start → out={}", self.out.display());
        info!(
            "REUSE: inputs: mito={}, nuclear={}, paf(m2n)={}, paf(n2m)={}, bam(r2n)={}, bam(r2m)={}",
            self.mito.display(),
            self.nuclear.display(),
            self.paf_mito_to_nuc.display(),
            self.paf_nuc_to_mito.display(),
            self.bam_reads_to_nuc.display(),
            self.bam_reads_to_mito.display()
        );

        fs::create_dir_all(&self.out)?;

        // resolve samtools path
        let samtools_path: PathBuf = if let Some(p) = &self.samtools_bin {
            p.clone()
        } else {
            which::which("samtools").map_err(|_| {
                anyhow::anyhow!("samtools not found: pass --samtools-bin or add to PATH")
            })?
        };
        log::info!("REUSE: using samtools at {}", samtools_path.display());

        // 1) Load PAFs
        let t_paf = Instant::now();
        // TODO: min_id and min_len should be cli params
        let m2n = paf::read_paf(&self.paf_mito_to_nuc, self.min_id, self.min_len)?;
        let n2m = paf::read_paf(&self.paf_nuc_to_mito, self.min_id, self.min_len)?;
        info!(
            "REUSE: PAFs loaded in {:.2}s (m2n={}, n2m={})",
            t_paf.elapsed().as_secs_f32(),
            m2n.len(),
            n2m.len()
        );

        // 2) Pair/merge loci
        let t_pair = Instant::now();
        let pairs0 = paf::pair_and_merge(&m2n, n2m, 50)?;
        info!(
            "REUSE: paired {} candidate loci in {:.2}s",
            pairs0.len(),
            t_pair.elapsed().as_secs_f32()
        );

        // 2.5) Exclude nuclear contigs that are essentially the mito assembly duplicated inside nuclear
        // (unless explicitly allowed via --allow-organelle-in-nuclear).
        let pairs = if self.allow_organelle_in_nuclear {
            log::warn!(
                "REUSE: keeping organelle-like nuclear contigs (--allow-organelle-in-nuclear)."
            );
            pairs0
        } else {
            let mito_lens = summary::contig_lengths(&self.mito)?;
            let nuc_lens = summary::contig_lengths(&self.nuclear)?;
            // Detect from mito→nuc hits (m2n) with strict criteria (e.g., ≥99.5% ident, ≥95% mito coverage, ~size parity)
            let flagged: HashSet<String> =
                paf::detect_organelle_in_nuclear(&m2n, &mito_lens, &nuc_lens);

            if !flagged.is_empty() {
                log::warn!(
            "REUSE: excluding {} nuclear contig(s) that look like full mito assemblies: {}{}",
            flagged.len(),
            flagged.iter().take(5).cloned().collect::<Vec<_>>().join(", "),
            if flagged.len() > 5 { ", …" } else { "" }
        );
                // Write for provenance
                fs::write(
                    self.out.join("excluded_organelle_like_nuclear_contigs.txt"),
                    flagged.iter().cloned().collect::<Vec<_>>().join("\n"),
                )?;
            }

            pairs0
                .into_iter()
                .filter(|p| !flagged.contains(&p.nuc_contig))
                .collect()
        };

        info!(
            "REUSE: {} candidate loci after filtering organelle-like nuclear contigs.",
            pairs.len()
        );

        // 3) Coverage + spans
        let t_cov = Instant::now();
        info!(
            "REUSE: coverage+spans for {} pairs (flank=500, win=250)…",
            pairs.len()
        );

        let (coverage, spans) = bam::compute_coverage_and_spans_with_tools(
            &self.bam_reads_to_nuc,
            &self.bam_reads_to_mito,
            &pairs,
            500,
            250,
            &samtools_path,
        )?;

        info!(
            "REUSE: coverage+spans done in {:.2}s (median depth: nuclear={:.3}, mito={:.3})",
            t_cov.elapsed().as_secs_f32(),
            coverage.nuclear_median,
            coverage.mito_median
        );

        // 4) Scoring
        let t_score = Instant::now();
        let embed = paf::default_embed_stub(&pairs); // placeholder
        let weights = Weights::default();
        let params = ClassifyParams {
            call_threshold: self.call_threshold,
            highconf_threshold: 0.30,
        };
        info!(
            "REUSE: scoring {} pairs (threshold={:.3})…",
            pairs.len(),
            params.call_threshold
        );
        let (pairs_tsv, classes_tsv) =
            scoring::classify_pairs(&pairs, &coverage, &spans, &embed, weights, params)?;
        info!(
            "REUSE: scoring finished in {:.2}s",
            t_score.elapsed().as_secs_f32()
        );

        // 5) Write outputs
        let t_io = Instant::now();
        fs::write(self.out.join("pairs.tsv"), &pairs_tsv)?;
        fs::write(self.out.join("classification.tsv"), &classes_tsv)?;
        info!(
            "REUSE: wrote outputs in {:.2}s",
            t_io.elapsed().as_secs_f32()
        );

        // Optional tiny summary
        let counts = count_calls_tsv(&classes_tsv);
        info!(
            "REUSE: calls — NUMT: {}, NIMT: {}, Ambiguous: {}",
            counts.numt, counts.nimt, counts.ambig
        );

        info!("REUSE: done in {:.2}s total", t0.elapsed().as_secs_f32());

        let calls = parse_calls_tsv(&classes_tsv);
        let summary = summary::compute_percentages(&self.mito, &self.nuclear, &pairs, &calls)?;
        summary::write_summary_tsv(&self.out.join("summary.tsv"), &summary)?;
        log::info!(
            "SUMMARY: NUMT in nuclear = {:.4}%, NIMT in mito = {:.4}% ; mito covered by NUMT homologs = {:.4}%",
            summary.nuc_pct_numt, summary.mito_pct_nimt, summary.mito_pct_covered_by_numt_homologs
        );

        Ok(())
    }
}

/// Minimal helper so we can log a call summary without parsing files again.
struct CallCounts {
    numt: usize,
    nimt: usize,
    ambig: usize,
}
fn count_calls_tsv(s: &str) -> CallCounts {
    let mut c = CallCounts {
        numt: 0,
        nimt: 0,
        ambig: 0,
    };
    for line in s.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        // columns: pair_id, call, confidence, reason
        let mut it = line.split('\t');
        let _ = it.next(); // pair_id
        if let Some(call) = it.next() {
            match call {
                "Likely_NUMT" => c.numt += 1,
                "Likely_NIMT" => c.nimt += 1,
                _ => c.ambig += 1,
            }
        }
    }
    c
}
