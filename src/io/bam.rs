//! BAM coverage & spanning-read metrics using noodles (bam 0.83, sam 0.79).
//! Requires BAM + BAI next to the BAM.

use anyhow::Result;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::MappingQuality;
use std::collections::HashMap;
use std::path::Path;

use crate::io::paf::PairedLocus;

use std::process::Command;
use tempfile::tempdir;

const SPAN_CAP_BP: i32 = 20_000;

/// Overall coverage medians + per-pair local depths.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CoverageSummary {
    pub nuclear_median: f64,
    pub mito_median: f64,
    /// pair_id -> (D_nuc_loc, D_mito_loc)
    pub per_pair: HashMap<String, (f32, f32)>,
}

/// Per-pair spanning support fractions.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, Default)]
pub struct SpanSummary {
    /// pair_id -> (S_nuc_span, S_mito_span)
    pub per_pair: HashMap<String, (f32, f32)>,
}

#[derive(Clone, Copy)]
struct Window {
    start: i32, // 0-based inclusive
    end: i32,   // 0-based exclusive
}

fn avg_depth_for_region_path_with_samtools(
    bam_path: &Path,
    rname: &str,
    window: Window,
    samtools_bin: &Path,
) -> Result<f32> {
    let s1 = (window.start.max(0) + 1) as usize;
    let e1 = window.end.max(1) as usize;
    let region = format!("{rname}:{s1}-{e1}");

    let dir = tempdir()?;
    let out_path = dir.path().join("slice.bam");

    let output = Command::new(samtools_bin)
        .args(["view", "-b", "-o"])
        .arg(&out_path)
        .arg(bam_path)
        .arg(&region)
        .output()
        .map_err(|e| anyhow::anyhow!("failed to spawn {} view: {}", samtools_bin.display(), e))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow::anyhow!(
            "samtools view failed for region {} on {}: {}",
            region,
            bam_path.display(),
            stderr.trim()
        ));
    }

    let meta = fs_err::metadata(&out_path)
        .map_err(|e| anyhow::anyhow!("slice not created for {}: {}", region, e))?;
    if meta.len() == 0 {
        return Ok(0.0);
    }

    let mut reader = noodles_bam::io::Reader::new(fs_err::File::open(&out_path)?);
    let _header = reader.read_header()?;

    let win_len = (window.end - window.start).max(1) as i64;
    let mut covered: i64 = 0;

    for result in reader.records() {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }

        let rec_start0 = match record.alignment_start() {
            Some(Ok(p)) => (usize::from(p) as i32 - 1).max(0),
            _ => continue,
        };

        let mut ref_len = 0i32;
        for op in record.cigar().iter() {
            let op = op?;
            match op.kind() {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => ref_len += i32::try_from(op.len()).unwrap_or(0),
                _ => {}
            }
        }
        if ref_len <= 0 {
            continue;
        }
        let rec_end0 = rec_start0 + ref_len;

        let ov_start = rec_start0.max(window.start);
        let ov_end = rec_end0.min(window.end);
        if ov_end <= ov_start {
            continue;
        }

        let mut ref_pos = rec_start0;
        for op in record.cigar().iter() {
            let op = op?;
            let oplen = i32::try_from(op.len()).unwrap_or(0);
            let consumes_ref = matches!(
                op.kind(),
                Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion
                    | Kind::Skip
            );
            let contributes = matches!(
                op.kind(),
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
            );
            if consumes_ref {
                let seg_start = ref_pos;
                let seg_end = ref_pos + oplen;
                if contributes {
                    let s = seg_start.max(ov_start);
                    let e = seg_end.min(ov_end);
                    if e > s {
                        covered += (e - s) as i64;
                    }
                }
                ref_pos += oplen;
            }
        }
    }

    Ok((covered as f32) / (win_len as f32))
}

fn span_fraction_for_region_path_with_samtools(
    bam_path: &Path,
    rname: &str,
    window: Window, // full candidate window (can be huge)
    samtools_bin: &Path,
    span_cap_bp: i32, // NEW: cap for span check, e.g. 20_000
) -> Result<f32> {
    const MIN_FRAC: f32 = 0.80;
    const MIN_MAPQ: u8 = 20;

    // Build a capped window around the midpoint for the span calculation
    let mid = (window.start + window.end) / 2;
    let half = (span_cap_bp / 2).max(1);
    let span_w = Window {
        start: (mid - half).max(0),
        end: mid + half,
    };

    // Slice only the capped region (faster than slicing the full mega-window)
    let s1 = (span_w.start.max(0) + 1) as usize; // 1-based inclusive
    let e1 = span_w.end.max(1) as usize;
    let region = format!("{rname}:{s1}-{e1}");

    let dir = tempdir()?;
    let out_path = dir.path().join("slice.bam");

    let output = Command::new(samtools_bin)
        .args(["view", "-b", "-o"])
        .arg(&out_path)
        .arg(bam_path)
        .arg(&region)
        .output()
        .map_err(|e| anyhow::anyhow!("failed to spawn {} view: {}", samtools_bin.display(), e))?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow::anyhow!(
            "samtools view failed for region {} on {}: {}",
            region,
            bam_path.display(),
            stderr.trim()
        ));
    }

    let meta = fs_err::metadata(&out_path)
        .map_err(|e| anyhow::anyhow!("slice not created for {}: {}", region, e))?;
    if meta.len() == 0 {
        return Ok(0.0);
    }

    let mut reader = noodles_bam::io::Reader::new(fs_err::File::open(&out_path)?);
    let _header = reader.read_header()?;

    let win_len = (span_w.end - span_w.start).max(1) as f32; // use capped window length
    let mut total = 0f32;
    let mut spans = 0f32;

    for result in reader.records() {
        let record = result?;
        if record.flags().is_unmapped() {
            continue;
        }
        if let Some(q) = record.mapping_quality() {
            if q < MappingQuality::new(MIN_MAPQ).unwrap() {
                continue;
            }
        }

        // 1-based -> 0-based
        let rec_start0 = match record.alignment_start() {
            Some(Ok(p)) => (usize::from(p) as i32 - 1).max(0),
            _ => continue,
        };

        // Reference span from ref-consuming CIGAR ops
        let mut ref_len = 0i32;
        for op in record.cigar().iter() {
            let op = op?;
            match op.kind() {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => ref_len += i32::try_from(op.len()).unwrap_or(0),
                _ => {}
            }
        }
        if ref_len <= 0 {
            continue;
        }
        let rec_end0 = rec_start0 + ref_len;

        // Overlap against the **capped** span window
        let ov_start = rec_start0.max(span_w.start);
        let ov_end = rec_end0.min(span_w.end);
        let overlap = (ov_end - ov_start).max(0) as f32;

        total += 1.0;
        if overlap / win_len >= MIN_FRAC {
            spans += 1.0;
        }
    }

    Ok(if total == 0.0 { 0.0 } else { spans / total })
}

pub fn compute_coverage_and_spans_with_tools(
    bam_r2n: &Path,
    bam_r2m: &Path,
    pairs: &[PairedLocus],
    flank: u32,
    _win_bp: u32,
    samtools_bin: &Path,
) -> Result<(CoverageSummary, SpanSummary)> {
    use std::time::Instant;
    let t0 = Instant::now();
    log::info!(
        "BAM: computing coverage & spans for {} pairs (flank={} bp) using samtools={}",
        pairs.len(),
        flank,
        samtools_bin.display()
    );

    let mut per_pair = std::collections::HashMap::new();
    let mut per_span = std::collections::HashMap::new();
    let mut nuc_locals = Vec::with_capacity(pairs.len());
    let mut mito_locals = Vec::with_capacity(pairs.len());

    let n = pairs.len().max(1);
    let mut last_tick = Instant::now();

    for (i, p) in pairs.iter().enumerate() {
        let nuc_w = Window {
            start: p.nuc_start as i32 - flank as i32,
            end: p.nuc_end as i32 + flank as i32,
        };
        let mito_w = Window {
            start: p.mito_start as i32 - flank as i32,
            end: p.mito_end as i32 + flank as i32,
        };

        let d_nuc =
            avg_depth_for_region_path_with_samtools(bam_r2n, &p.nuc_contig, nuc_w, samtools_bin)?;
        let d_mito =
            avg_depth_for_region_path_with_samtools(bam_r2m, &p.mito_contig, mito_w, samtools_bin)?;
        let s_nuc = span_fraction_for_region_path_with_samtools(
            bam_r2n,
            &p.nuc_contig,
            nuc_w,
            samtools_bin,
            SPAN_CAP_BP,
        )?;
        let s_mito = span_fraction_for_region_path_with_samtools(
            bam_r2m,
            &p.mito_contig,
            mito_w,
            samtools_bin,
            SPAN_CAP_BP,
        )?;

        per_pair.insert(p.pair_id.clone(), (d_nuc, d_mito));
        per_span.insert(p.pair_id.clone(), (s_nuc, s_mito));
        nuc_locals.push(d_nuc);
        mito_locals.push(d_mito);

        if last_tick.elapsed().as_secs_f32() > 2.0 {
            let done = i + 1;
            let frac = (done as f32) / (n as f32);
            let elapsed = t0.elapsed().as_secs_f32();
            let eta = if frac > 0.0 {
                elapsed * (1.0 - frac) / frac
            } else {
                f32::INFINITY
            };
            log::info!(
                "BAM: {}/{} pairs ({:.0}%) processed, elapsed {:.1}s, ETA ~{:.1}s",
                done,
                n,
                100.0 * frac,
                elapsed,
                eta
            );
            last_tick = Instant::now();
        }
    }

    fn median(v: &mut [f32]) -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let m = v.len();
        if m % 2 == 1 {
            v[m / 2] as f64
        } else {
            (v[m / 2 - 1] as f64 + v[m / 2] as f64) / 2.0
        }
    }
    let mut nt = nuc_locals.clone();
    let mut mt = mito_locals.clone();
    let cov = CoverageSummary {
        nuclear_median: median(&mut nt),
        mito_median: median(&mut mt),
        per_pair,
    };
    let spans = SpanSummary { per_pair: per_span };

    log::info!(
        "BAM: finished {} pairs in {:.1}s (medians: nuclear={:.3}, mito={:.3})",
        n,
        t0.elapsed().as_secs_f32(),
        cov.nuclear_median,
        cov.mito_median
    );
    Ok((cov, spans))
}
