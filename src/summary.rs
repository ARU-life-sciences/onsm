//! Summaries: genome-size-normalized coverage of called NUMTs/NIMTs.

use anyhow::{anyhow, Result};
use fs_err as fs;
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::path::Path;

use crate::io::paf::PairedLocus;

/// Load total sequence length for each contig from a FASTA (supports .gz).
pub fn contig_lengths(fasta: &Path) -> Result<HashMap<String, u64>> {
    if !fasta.exists() {
        return Err(anyhow!("FASTA not found: {}", fasta.display()));
    }
    let mut rdr =
        parse_fastx_file(fasta).map_err(|e| anyhow!("open FASTA {}: {}", fasta.display(), e))?;
    let mut lens = HashMap::new();
    while let Some(rec) = rdr.next() {
        let rec = rec.map_err(|e| anyhow!("read FASTA {}: {}", fasta.display(), e))?;
        let name = String::from_utf8(rec.id().to_vec()).unwrap_or_else(|_| "contig".to_string());
        lens.entry(name)
            .and_modify(|v| *v += rec.seq().len() as u64)
            .or_insert(rec.seq().len() as u64);
    }
    Ok(lens)
}

/// Merge intervals on a contig and return merged bp.
fn merged_bp(intervals: &[(u32, u32)]) -> u64 {
    if intervals.is_empty() {
        return 0;
    }
    // sort by start
    let mut v = intervals.to_vec();
    v.sort_unstable_by_key(|x| x.0);
    let mut merged: u64 = 0;
    let mut cur = v[0];
    for &(s, e) in v.iter().skip(1) {
        if s <= cur.1 {
            // overlap/touch
            if e > cur.1 {
                cur.1 = e;
            }
        } else {
            merged += (cur.1 - cur.0) as u64;
            cur = (s, e);
        }
    }
    merged += (cur.1 - cur.0) as u64;
    merged
}

/// Genome-level summary of called NUMTs/NIMTs and homolog coverage.
pub struct PercentSummary {
    pub n_pairs: usize,
    pub n_numt: usize,
    pub n_nimt: usize,

    pub nuclear_total_bp: u64,
    pub mito_total_bp: u64,

    // Occupancy (how much of each genome is *occupied* by called insertions)
    pub nuc_bp_numt: u64,
    pub nuc_pct_numt: f64,
    pub mito_bp_nimt: u64,
    pub mito_pct_nimt: f64,

    // Homolog coverage (NUCmer-style: how much of one genome is represented in the other)
    pub mito_bp_covered_by_numt_homologs: u64,
    pub mito_pct_covered_by_numt_homologs: f64,
    pub nuc_bp_covered_by_nimt_homologs: u64,
    pub nuc_pct_covered_by_nimt_homologs: f64,
}

pub fn compute_percentages(
    mito_fa: &Path,
    nuc_fa: &Path,
    pairs: &[PairedLocus],
    calls: &HashMap<String, String>, // pair_id -> "Likely_NUMT"/"Likely_NIMT"/...
) -> Result<PercentSummary> {
    let mito_lens = contig_lengths(mito_fa)?;
    let nuc_lens = contig_lengths(nuc_fa)?;
    let mito_tot = mito_lens.values().copied().sum::<u64>().max(1);
    let nuc_tot = nuc_lens.values().copied().sum::<u64>().max(1);

    // Collect intervals
    let mut nuc_numt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut mito_nimt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut mito_from_numt_side: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut nuc_from_nimt_side: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    let mut n_numt = 0usize;
    let mut n_nimt = 0usize;

    for p in pairs {
        match calls.get(&p.pair_id).map(|s| s.as_str()) {
            Some("Likely_NUMT") => {
                n_numt += 1;
                nuc_numt
                    .entry(p.nuc_contig.clone())
                    .or_default()
                    .push((p.nuc_start, p.nuc_end));
                mito_from_numt_side
                    .entry(p.mito_contig.clone())
                    .or_default()
                    .push((p.mito_start, p.mito_end));
            }
            Some("Likely_NIMT") => {
                n_nimt += 1;
                mito_nimt
                    .entry(p.mito_contig.clone())
                    .or_default()
                    .push((p.mito_start, p.mito_end));
                nuc_from_nimt_side
                    .entry(p.nuc_contig.clone())
                    .or_default()
                    .push((p.nuc_start, p.nuc_end));
            }
            _ => {}
        }
    }

    let sum_merged = |m: &HashMap<String, Vec<(u32, u32)>>| -> u64 {
        m.iter().map(|(_ctg, ivs)| merged_bp(ivs)).sum()
    };

    let nuc_bp_numt = sum_merged(&nuc_numt);
    let mito_bp_nimt = sum_merged(&mito_nimt);
    let mito_bp_cov_by_numt = sum_merged(&mito_from_numt_side);
    let nuc_bp_cov_by_nimt = sum_merged(&nuc_from_nimt_side);

    Ok(PercentSummary {
        n_pairs: pairs.len(),
        n_numt,
        n_nimt,

        nuclear_total_bp: nuc_tot,
        mito_total_bp: mito_tot,

        nuc_bp_numt,
        nuc_pct_numt: (nuc_bp_numt as f64) * 100.0 / (nuc_tot as f64),
        mito_bp_nimt,
        mito_pct_nimt: (mito_bp_nimt as f64) * 100.0 / (mito_tot as f64),

        mito_bp_covered_by_numt_homologs: mito_bp_cov_by_numt,
        mito_pct_covered_by_numt_homologs: (mito_bp_cov_by_numt as f64) * 100.0 / (mito_tot as f64),
        nuc_bp_covered_by_nimt_homologs: nuc_bp_cov_by_nimt,
        nuc_pct_covered_by_nimt_homologs: (nuc_bp_cov_by_nimt as f64) * 100.0 / (nuc_tot as f64),
    })
}

pub fn write_summary_tsv(path: &Path, s: &PercentSummary) -> Result<()> {
    let mut out = String::new();
    out.push_str("metric\tvalue\n");
    out.push_str(&format!("n_pairs\t{}\n", s.n_pairs));
    out.push_str(&format!("n_numt\t{}\n", s.n_numt));
    out.push_str(&format!("n_nimt\t{}\n", s.n_nimt));

    out.push_str(&format!("nuclear_bp_total\t{}\n", s.nuclear_total_bp));
    out.push_str(&format!("nuclear_bp_numt\t{}\n", s.nuc_bp_numt));
    out.push_str(&format!("nuclear_pct_numt\t{:.6}\n", s.nuc_pct_numt));

    out.push_str(&format!("mito_bp_total\t{}\n", s.mito_total_bp));
    out.push_str(&format!("mito_bp_nimt\t{}\n", s.mito_bp_nimt));
    out.push_str(&format!("mito_pct_nimt\t{:.6}\n", s.mito_pct_nimt));

    out.push_str(&format!(
        "mito_bp_covered_by_numt_homologs\t{}\n",
        s.mito_bp_covered_by_numt_homologs
    ));
    out.push_str(&format!(
        "mito_pct_covered_by_numt_homologs\t{:.6}\n",
        s.mito_pct_covered_by_numt_homologs
    ));

    out.push_str(&format!(
        "nuc_bp_covered_by_nimt_homologs\t{}\n",
        s.nuc_bp_covered_by_nimt_homologs
    ));
    out.push_str(&format!(
        "nuc_pct_covered_by_nimt_homologs\t{:.6}\n",
        s.nuc_pct_covered_by_nimt_homologs
    ));
    fs::write(path, out)?;
    Ok(())
}
