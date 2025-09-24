use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;
use std::process::Command;

use crate::model::{CoverageSummary, PairedLocus, SpanSummary};

/// Half-open window on reference in 0-based coordinates [start, end).
#[derive(Debug, Clone, Copy)]
pub struct Window {
    pub start: i32,
    pub end: i32,
}

fn region_str(rname: &str, w: Window) -> String {
    // samtools uses 1-based inclusive coordinates
    let s1 = (w.start.max(0) + 1) as usize;
    let e1 = w.end.max(w.start + 1) as usize;
    format!("{rname}:{s1}-{e1}")
}

fn parse_cigar_ref_consumed(cigar: &str) -> Option<u32> {
    // Sum of ref-consuming ops: M, =, X, D, N
    let mut num = 0u64;
    let mut acc = 0u64;
    for ch in cigar.bytes() {
        match ch {
            b'0'..=b'9' => {
                num = num * 10 + (ch - b'0') as u64;
            }
            b'M' | b'=' | b'X' | b'D' | b'N' => {
                acc = acc.saturating_add(num);
                num = 0;
            }
            b'I' | b'S' | b'H' | b'P' => {
                num = 0;
            }
            _ => return None, // malformed
        }
    }
    Some(acc.min(u32::MAX as u64) as u32)
}

fn median_f32(mut v: Vec<f32>) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Compute local median depth in a region using `samtools depth`.
fn local_median_depth(samtools: &Path, bam: &Path, rname: &str, w: Window) -> Result<f32> {
    let region = region_str(rname, w);
    let out = Command::new(samtools)
        .args(["depth", "-r"])
        .arg(&region)
        .arg(bam)
        .output()
        .with_context(|| format!("spawn samtools depth for {region}"))?;
    if !out.status.success() {
        let err = String::from_utf8_lossy(&out.stderr);
        return Err(anyhow::anyhow!("samtools depth failed: {}", err.trim()));
    }
    // depth output: chrom  pos  depth
    let mut depths = Vec::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let mut it = line.split_whitespace();
        let _chrom = it.next();
        let _pos = it.next();
        if let Some(d) = it.next() {
            if let Ok(x) = d.parse::<u32>() {
                depths.push(x as f32);
            }
        }
    }
    Ok(median_f32(depths))
}

/// Fraction of alignments that span the entire [w.start, w.end) window on rname.
/// Uses `samtools view` (SAM text), MAPQ ≥ 20.
fn span_fraction(samtools: &Path, bam: &Path, rname: &str, w: Window) -> Result<f32> {
    const MIN_MAPQ: u8 = 20;
    let region = region_str(rname, w);
    let out = Command::new(samtools)
        .args(["view"])
        .arg(bam)
        .arg(&region)
        .output()
        .with_context(|| format!("spawn samtools view for {region}"))?;
    if !out.status.success() {
        let err = String::from_utf8_lossy(&out.stderr);
        return Err(anyhow::anyhow!("samtools view failed: {}", err.trim()));
    }
    let s1 = w.start.max(0) + 1; // window start 1-based
    let e1 = w.end.max(w.start + 1); // window end 1-based inclusive-ish

    let mut total = 0f32;
    let mut spans = 0f32;

    for line in String::from_utf8_lossy(&out.stdout).lines() {
        if line.is_empty() || line.starts_with('@') {
            continue;
        }
        let mut cols = line.split('\t');
        let _qname = cols.next();
        let flag = cols.next().and_then(|s| s.parse::<u16>().ok()).unwrap_or(0);
        let rname_sam = cols.next().unwrap_or("*");
        let pos = cols.next().and_then(|s| s.parse::<i32>().ok()).unwrap_or(0);
        let mapq = cols.next().and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
        let cigar = cols.next().unwrap_or("*");

        // filter
        if (flag & 0x4) != 0 {
            continue; // unmapped
        }
        if mapq < MIN_MAPQ {
            continue;
        }
        if rname_sam != rname {
            continue;
        }

        let ref_len = match parse_cigar_ref_consumed(cigar) {
            Some(x) if x > 0 => x as i32,
            _ => continue,
        };
        let rec_start = pos; // POS is 1-based
        let rec_end = pos + ref_len - 1; // inclusive on reference

        total += 1.0;
        if rec_start <= s1 && rec_end >= e1 {
            spans += 1.0;
        }
    }

    Ok(if total == 0.0 { 0.0 } else { spans / total })
}

/// Compute (coverage, spans) for all pairs using small windows around each locus.
/// Global medians are computed as the median of per-pair local medians (robust & fast).
pub fn compute_coverage_and_spans_with_tools(
    bam_reads_to_nuc: &Path,
    bam_reads_to_mito: &Path,
    pairs: &[PairedLocus],
    flank: u32,
    win: u32,
    samtools: &Path,
) -> Result<(CoverageSummary, SpanSummary)> {
    log::info!(
        "BAM: computing coverage & spans for {} pairs (flank={} bp) using samtools={}",
        pairs.len(),
        flank,
        samtools.display()
    );

    let mut per_pair_depth: HashMap<String, (f32, f32)> = HashMap::new();
    let mut per_pair_span: HashMap<String, (f32, f32)> = HashMap::new();

    let mut nuc_locals = Vec::new();
    let mut mito_locals = Vec::new();

    let flank_i = flank as i32;
    let win_i = win as i32;

    for (i, p) in pairs.iter().enumerate() {
        if (i + 1) % 50 == 0 || i == 0 {
            log::info!("BAM: {}/{} …", i + 1, pairs.len());
        }

        // Center windows at the alignment midpoints
        let n_mid = ((p.nuc_start + p.nuc_end) / 2) as i32;
        let m_mid = ((p.mito_start + p.mito_end) / 2) as i32;
        let n_w = Window {
            start: n_mid - flank_i,
            end: n_mid + flank_i,
        };
        let m_w = Window {
            start: m_mid - flank_i,
            end: m_mid + flank_i,
        };

        // Local depths
        let d_n = local_median_depth(samtools, bam_reads_to_nuc, &p.nuc_contig, n_w)?;
        let d_m = local_median_depth(samtools, bam_reads_to_mito, &p.mito_contig, m_w)?;
        per_pair_depth.insert(p.pair_id.clone(), (d_n, d_m));
        nuc_locals.push(d_n);
        mito_locals.push(d_m);

        // Spanning windows: tighten to ±win around mid (must fully cover)
        let n_s = Window {
            start: n_mid - win_i,
            end: n_mid + win_i,
        };
        let m_s = Window {
            start: m_mid - win_i,
            end: m_mid + win_i,
        };
        let s_n = span_fraction(samtools, bam_reads_to_nuc, &p.nuc_contig, n_s)?;
        let s_m = span_fraction(samtools, bam_reads_to_mito, &p.mito_contig, m_s)?;
        per_pair_span.insert(p.pair_id.clone(), (s_n, s_m));
    }

    let nuclear_median = super::bam::median_f32(nuc_locals) as f64;
    let mito_median = super::bam::median_f32(mito_locals) as f64;

    Ok((
        CoverageSummary {
            nuclear_median,
            mito_median,
            per_pair: per_pair_depth,
        },
        SpanSummary {
            per_pair: per_pair_span,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_ref_len_parses_basic() {
        assert_eq!(parse_cigar_ref_consumed("100M"), Some(100));
        assert_eq!(parse_cigar_ref_consumed("10S90M"), Some(90));
        assert_eq!(parse_cigar_ref_consumed("50M10I40M"), Some(90));
        assert_eq!(parse_cigar_ref_consumed("50M5D45M"), Some(100));
        assert_eq!(parse_cigar_ref_consumed("50M100N50M"), Some(200)); // spliced
        assert_eq!(parse_cigar_ref_consumed("*"), None);
    }

    #[test]
    fn median_works() {
        assert_eq!(median_f32(vec![]), 0.0);
        assert_eq!(median_f32(vec![1.0]), 1.0);
        assert_eq!(median_f32(vec![1.0, 3.0]), 2.0);
        assert_eq!(median_f32(vec![1.0, 3.0, 2.0]), 2.0);
    }

    #[test]
    fn region_format_ok() {
        let r = region_str("chr1", Window { start: 0, end: 10 });
        assert_eq!(r, "chr1:1-10");
        let r = region_str("chr1", Window { start: 9, end: 20 });
        assert_eq!(r, "chr1:10-20");
    }
}
