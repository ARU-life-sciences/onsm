use anyhow::Result;
use std::collections::HashMap;

use crate::model::{ClassifyParams, PairedLocus, Weights};
use crate::model::{CoverageSummary, SpanSummary};
use std::fmt::Write as _;

#[derive(Debug, Clone)]
#[allow(clippy::upper_case_acronyms)]
enum Call {
    NUMT,
    NIMT,
    Ambiguous,
}

impl Call {
    fn as_str_and_reason(&self) -> (&'static str, &'static str) {
        match self {
            Call::NUMT => ("Likely_NUMT", "score_difference"),
            Call::NIMT => ("Likely_NIMT", "score_difference"),
            Call::Ambiguous => ("Ambiguous", "delta_below_threshold"),
        }
    }
}

fn clamp01(x: f32) -> f32 {
    x.clamp(0.0, 1.0)
}

// Soft saturating map of length (bp) to [0..1].
fn scale_len(len_bp: u32) -> f32 {
    let l = len_bp as f32;
    let l50 = 5_000.0f32;
    clamp01(l / (l + l50))
}

pub fn classify_pairs(
    pairs: &[PairedLocus],
    coverage: &CoverageSummary,
    spans: &SpanSummary,
    w: Weights,
    params: ClassifyParams,
) -> Result<(String, String)> {
    // lookups
    let depth_map: HashMap<&str, (f32, f32)> = coverage
        .per_pair
        .iter()
        .map(|(k, v)| (k.as_str(), *v))
        .collect();
    let span_map: HashMap<&str, (f32, f32)> = spans
        .per_pair
        .iter()
        .map(|(k, v)| (k.as_str(), *v))
        .collect();

    let dn_med = coverage.nuclear_median as f32;
    let dm_med = coverage.mito_median as f32;

    let mut pairs_tsv = String::from(
        "pair_id\tnuc_contig\tnuc_start\tnuc_end\tmito_contig\tmito_start\tmito_end\taln_len\taln_ident\trnuc\trmito\ts_nuc\ts_mito\tscore_numt\tscore_nimt\n"
    );
    let mut class_tsv = String::from("pair_id\tcall\tconfidence\treason_codes\n");

    for p in pairs {
        let a = clamp01(p.aln_ident);
        let l = scale_len(p.aln_len);
        let base = w.w_a * a + w.w_l * l;

        let (d_n_loc, d_m_loc) = depth_map
            .get(p.pair_id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));
        // normalized local medians
        let rnuc = if dn_med > 0.0 { d_n_loc / dn_med } else { 0.0 };
        let rmito = if dm_med > 0.0 { d_m_loc / dm_med } else { 0.0 };

        // Depth consistency terms (favor ~1.0)
        let d_numt = clamp01(1.0 - (rnuc - 1.0).abs());
        let d_nimt = clamp01(1.0 - (rmito - 1.0).abs());

        // Spanning
        let (s_nuc, s_mito) = span_map
            .get(p.pair_id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));

        // Contrast boosters (signed): + favors NUMT, − favors NIMT
        let eps = 1e-3_f32;
        let log2_ratio = ((rnuc + eps) / (rmito + eps)).ln() / std::f32::consts::LN_2;
        let depth_contrast = (1.25 * log2_ratio).tanh(); // (-1..1)
        let span_contrast = s_nuc - s_mito; // (-1..1)

        // Build scores
        let pro_numt = w.w_d * d_numt + w.w_s * s_nuc;
        let pro_nimt = w.w_d * d_nimt + w.w_s * s_mito;
        let pen_numt = w.w_d * d_nimt + w.w_s * s_mito;
        let pen_nimt = w.w_d * d_numt + w.w_s * s_nuc;
        let boost_numt = w.w_d * depth_contrast + w.w_s * span_contrast;
        let boost_nimt = -w.w_d * depth_contrast - w.w_s * span_contrast;

        let score_numt = base + pro_numt - pen_numt + boost_numt;
        let score_nimt = base + pro_nimt - pen_nimt + boost_nimt;

        let diff = score_numt - score_nimt;
        let delta = diff.abs();
        let call = if diff >= params.call_threshold {
            Call::NUMT
        } else if -diff >= params.call_threshold {
            Call::NIMT
        } else {
            Call::Ambiguous
        };

        let (call_str, reason) = call.as_str_and_reason();

        let _ = writeln!(
            &mut pairs_tsv,
            "{pid}\t{nc}\t{ns}\t{ne}\t{mc}\t{ms}\t{me}\t{al}\t{ai:.4}\t{rn:.3}\t{rm:.3}\t{sn:.3}\t{sm:.3}\t{snmt:.4}\t{simt:.4}",
            pid = p.pair_id,
            nc = p.nuc_contig, ns = p.nuc_start, ne = p.nuc_end,
            mc = p.mito_contig, ms = p.mito_start, me = p.mito_end,
            al = p.aln_len, ai = a,
            rn = rnuc, rm = rmito,
            sn = s_nuc, sm = s_mito,
            snmt = score_numt, simt = score_nimt
        );
        let _ = writeln!(
            &mut class_tsv,
            "{pid}\t{call}\t{conf:.4}\t{reason}",
            pid = p.pair_id,
            call = call_str,
            conf = delta,
            reason = reason
        );
    }

    Ok((pairs_tsv, class_tsv))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{CoverageSummary, SpanSummary};

    #[test]
    fn favors_numt_when_nuclear_support_strong() {
        let pairs = vec![PairedLocus {
            pair_id: "P1".into(),
            nuc_contig: "chr1".into(),
            nuc_start: 100,
            nuc_end: 200,
            mito_contig: "m1".into(),
            mito_start: 50,
            mito_end: 150,
            aln_len: 5000,
            aln_ident: 0.98,
        }];
        let cov = CoverageSummary {
            nuclear_median: 30.0,
            mito_median: 30.0,
            per_pair: [("P1".into(), (30.0, 10.0))].into_iter().collect(), // rnuc=1.0, rmito=0.33
        };
        let spans = SpanSummary {
            per_pair: [("P1".into(), (0.8, 0.1))].into_iter().collect(),
        };
        let (pairs_tsv, class_tsv) = classify_pairs(
            &pairs,
            &cov,
            &spans,
            Weights::default(),
            ClassifyParams::default(),
        )
        .unwrap();
        assert!(pairs_tsv.contains("score_numt"));
        assert!(class_tsv.contains("Likely_NUMT"));
    }
}
