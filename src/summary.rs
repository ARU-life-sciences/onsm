//! Summaries: genome-level NUMT/NIMT percentages and “homologous bp” coverage.
//!
//! Given:
//!   - the mito and nuclear FASTA paths,
//!   - the paired loci (mito<->nuclear intervals),
//!   - and a call map {pair_id -> "Likely_NUMT" | "Likely_NIMT" | "Ambiguous"},
//!
//! we compute:
//!   * total assembly lengths (from FASTA),
//!   * union-lengths of loci on each side stratified by call type,
//!   * percentages (as PERCENT values; e.g., 0.0207 means 0.0207%).
//!
//! We treat PairedLocus coordinates as 0-based half-open [start, end).

use anyhow::Result;
use std::collections::HashMap;
use std::path::Path;

use crate::io::fasta;
use crate::model::PairedLocus;

/// Output struct that directly matches the `summary.tsv` rows you showed.
#[derive(Debug, Clone)]
pub struct Summary {
    pub n_pairs: usize,
    pub n_numt: usize,
    pub n_nimt: usize,

    pub nuclear_bp_total: u64,
    pub nuclear_bp_numt: u64,
    pub nuclear_pct_numt: f64,

    pub mito_bp_total: u64,
    pub mito_bp_nimt: u64,
    pub mito_pct_nimt: f64,

    // “Homologous coverage” on the opposite genome, stratified by call
    pub mito_bp_covered_by_numt_homologs: u64,
    pub mito_pct_covered_by_numt_homologs: f64,

    pub nuc_bp_covered_by_nimt_homologs: u64,
    pub nuc_pct_covered_by_nimt_homologs: f64,
}

/// Compute the summary for a run.
///
/// - `pairs`: candidate loci (reciprocal mapping + merging).
/// - `calls`: map of pair_id → call string ("Likely_NUMT", "Likely_NIMT", or other).
pub fn compute_percentages(
    mito_fa: &Path,
    nuc_fa: &Path,
    pairs: &[PairedLocus],
    calls: &HashMap<String, String>,
) -> Result<Summary> {
    // Assembly lengths
    let mito_lens = fasta::contig_lengths(mito_fa)?;
    let nuc_lens = fasta::contig_lengths(nuc_fa)?;
    let mito_bp_total: u64 = mito_lens.values().sum();
    let nuclear_bp_total: u64 = nuc_lens.values().sum();

    // Counters & per-contig interval buckets
    let mut n_numt = 0usize;
    let mut n_nimt = 0usize;

    // Intervals to union later, keyed by contig
    let mut nuc_intervals_numt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut mito_intervals_nimt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    // “Homologous coverage on the opposite genome”:
    //   NUMT calls contribute their *mito* intervals (coverage of mito by NUMT homologs)
    //   NIMT calls contribute their *nuclear* intervals (coverage of nuclear by NIMT homologs)
    let mut mito_intervals_from_numt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut nuc_intervals_from_nimt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

    for p in pairs {
        let call = calls
            .get(&p.pair_id)
            .map(String::as_str)
            .unwrap_or("Ambiguous");
        match call {
            "Likely_NUMT" => {
                n_numt += 1;

                // nuclear bp that are NUMT (union across nuclear side of these loci)
                add_interval(
                    &mut nuc_intervals_numt,
                    &p.nuc_contig,
                    p.nuc_start,
                    p.nuc_end,
                );

                // homologous coverage on mito (the counterpart region)
                add_interval(
                    &mut mito_intervals_from_numt,
                    &p.mito_contig,
                    p.mito_start,
                    p.mito_end,
                );
            }
            "Likely_NIMT" => {
                n_nimt += 1;

                // mito bp that are NIMT (union across mito side of these loci)
                add_interval(
                    &mut mito_intervals_nimt,
                    &p.mito_contig,
                    p.mito_start,
                    p.mito_end,
                );

                // homologous coverage on nuclear side
                add_interval(
                    &mut nuc_intervals_from_nimt,
                    &p.nuc_contig,
                    p.nuc_start,
                    p.nuc_end,
                );
            }
            _ => { /* Ambiguous – ignored for summary */ }
        }
    }

    // Union-lengths
    let nuclear_bp_numt = union_len_all(&nuc_intervals_numt);
    let mito_bp_nimt = union_len_all(&mito_intervals_nimt);

    let mito_bp_covered_by_numt_homologs = union_len_all(&mito_intervals_from_numt);
    let nuc_bp_covered_by_nimt_homologs = union_len_all(&nuc_intervals_from_nimt);

    // Percentages (as *percent* values, e.g., 0.0207 means 0.0207%)
    let nuclear_pct_numt = pct(nuclear_bp_numt, nuclear_bp_total);
    let mito_pct_nimt = pct(mito_bp_nimt, mito_bp_total);

    let mito_pct_covered_by_numt_homologs = pct(mito_bp_covered_by_numt_homologs, mito_bp_total);
    let nuc_pct_covered_by_nimt_homologs = pct(nuc_bp_covered_by_nimt_homologs, nuclear_bp_total);

    Ok(Summary {
        n_pairs: pairs.len(),
        n_numt,
        n_nimt,

        nuclear_bp_total,
        nuclear_bp_numt,
        nuclear_pct_numt,

        mito_bp_total,
        mito_bp_nimt,
        mito_pct_nimt,

        mito_bp_covered_by_numt_homologs,
        mito_pct_covered_by_numt_homologs,

        nuc_bp_covered_by_nimt_homologs,
        nuc_pct_covered_by_nimt_homologs,
    })
}

/// Write the summary as a 2-column TSV (metric\tvalue), mirroring your examples.
pub fn write_summary_tsv(out_path: &Path, s: &Summary) -> Result<()> {
    use std::fmt::Write;
    let mut t = String::new();

    writeln!(&mut t, "metric\tvalue")?;
    writeln!(&mut t, "n_pairs\t{}", s.n_pairs)?;
    writeln!(&mut t, "n_numt\t{}", s.n_numt)?;
    writeln!(&mut t, "n_nimt\t{}", s.n_nimt)?;
    writeln!(&mut t, "nuclear_bp_total\t{}", s.nuclear_bp_total)?;
    writeln!(&mut t, "nuclear_bp_numt\t{}", s.nuclear_bp_numt)?;
    writeln!(&mut t, "nuclear_pct_numt\t{:.6}", s.nuclear_pct_numt)?;
    writeln!(&mut t, "mito_bp_total\t{}", s.mito_bp_total)?;
    writeln!(&mut t, "mito_bp_nimt\t{}", s.mito_bp_nimt)?;
    writeln!(&mut t, "mito_pct_nimt\t{:.6}", s.mito_pct_nimt)?;
    writeln!(
        &mut t,
        "mito_bp_covered_by_numt_homologs\t{}",
        s.mito_bp_covered_by_numt_homologs
    )?;
    writeln!(
        &mut t,
        "mito_pct_covered_by_numt_homologs\t{:.6}",
        s.mito_pct_covered_by_numt_homologs
    )?;
    writeln!(
        &mut t,
        "nuc_bp_covered_by_nimt_homologs\t{}",
        s.nuc_bp_covered_by_nimt_homologs
    )?;
    writeln!(
        &mut t,
        "nuc_pct_covered_by_nimt_homologs\t{:.6}",
        s.nuc_pct_covered_by_nimt_homologs
    )?;

    fs_err::write(out_path, t)?;
    Ok(())
}

/// Parse the contents of classification.tsv (string) into a call map:
/// pair_id -> "Likely_NUMT" | "Likely_NIMT" | "Ambiguous" (or whatever is present).
pub fn parse_calls_tsv_str(s: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for line in s.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        if let (Some(pid), Some(call)) = (it.next(), it.next()) {
            m.insert(pid.to_string(), call.to_string());
        }
    }
    m
}

/// Convenience: parse classification.tsv from a file path.
pub fn parse_calls_tsv_file(path: &Path) -> Result<HashMap<String, String>> {
    let txt = fs_err::read_to_string(path)?;
    Ok(parse_calls_tsv_str(&txt))
}

/* ------------------------- internal helpers ------------------------- */

fn pct(numer: u64, denom: u64) -> f64 {
    if denom == 0 {
        0.0
    } else {
        100.0 * (numer as f64) / (denom as f64)
    }
}

/// Push an interval (0-based half-open) into a per-contig map.
fn add_interval(map: &mut HashMap<String, Vec<(u32, u32)>>, contig: &str, start: u32, end: u32) {
    if start >= end {
        return;
    }
    map.entry(contig.to_string())
        .or_default()
        .push((start, end));
}

/// Merge overlapping intervals for one contig and return total length.
fn union_len(mut v: Vec<(u32, u32)>) -> u64 {
    if v.is_empty() {
        return 0;
    }
    v.sort_by_key(|x| (x.0, x.1));
    let mut total: u64 = 0;
    let mut cur = v[0];

    for (s, e) in v.into_iter().skip(1) {
        if s <= cur.1 {
            // overlap or touch
            cur.1 = cur.1.max(e);
        } else {
            total += (cur.1 - cur.0) as u64;
            cur = (s, e);
        }
    }
    total += (cur.1 - cur.0) as u64;
    total
}

/// Union length across all contigs.
fn union_len_all(m: &HashMap<String, Vec<(u32, u32)>>) -> u64 {
    m.values().map(|v| union_len(v.clone())).sum()
}

/* ------------------------------ tests ------------------------------ */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_merges_overlaps() {
        // intervals: [0,10), [5,20), [30,40) => total = 20 + 10 = 30
        let v = vec![(0, 10), (5, 20), (30, 40)];
        assert_eq!(union_len(v), 30);
    }

    #[test]
    fn percentages_are_percent_values() {
        assert_eq!(pct(0, 1000), 0.0);
        assert!((pct(1, 1000) - 0.1).abs() < 1e-9); // 0.1%
        assert!((pct(500, 1000) - 50.0).abs() < 1e-9);
    }

    #[test]
    fn add_and_union_per_contig() {
        let mut m: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        add_interval(&mut m, "chr1", 0, 10);
        add_interval(&mut m, "chr1", 5, 20);
        add_interval(&mut m, "chr2", 0, 5);
        assert_eq!(union_len_all(&m), 25); // chr1=20, chr2=5
    }

    #[test]
    fn small_summary_counts() {
        // Three loci: one NUMT, one NIMT, one ambiguous
        let pairs = vec![
            PairedLocus {
                pair_id: "P1".into(),
                nuc_contig: "chr1".into(),
                nuc_start: 100,
                nuc_end: 200,
                mito_contig: "m1".into(),
                mito_start: 100,
                mito_end: 300,
                aln_len: 200,
                aln_ident: 0.99,
            },
            PairedLocus {
                pair_id: "P2".into(),
                nuc_contig: "chr1".into(),
                nuc_start: 500,
                nuc_end: 600,
                mito_contig: "m1".into(),
                mito_start: 400,
                mito_end: 450,
                aln_len: 100,
                aln_ident: 0.95,
            },
            PairedLocus {
                pair_id: "P3".into(),
                nuc_contig: "chr2".into(),
                nuc_start: 0,
                nuc_end: 50,
                mito_contig: "m1".into(),
                mito_start: 800,
                mito_end: 900,
                aln_len: 100,
                aln_ident: 0.90,
            },
        ];
        let calls: HashMap<_, _> = [
            ("P1".to_string(), "Likely_NUMT".to_string()),
            ("P2".to_string(), "Likely_NIMT".to_string()),
            ("P3".to_string(), "Ambiguous".to_string()),
        ]
        .into_iter()
        .collect();

        // Fake assembly sizes: mito=1000, chr1=1000, chr2=1000
        // We can’t read real FASTA here; instead exercise the interval logic.
        // So we only test the union logic directly:
        let mut mito_from_numt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        let mut nuc_from_nimt: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

        // P1 NUMT → contributes mito[100,300)
        add_interval(&mut mito_from_numt, "m1", 100, 300);
        // P2 NIMT → contributes nuc[chr1:500,600)
        add_interval(&mut nuc_from_nimt, "chr1", 500, 600);

        assert_eq!(union_len_all(&mito_from_numt), 200);
        assert_eq!(union_len_all(&nuc_from_nimt), 100);
    }
}
