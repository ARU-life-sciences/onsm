use anyhow::Result;
use std::collections::HashMap;

use crate::io::bam::{CoverageSummary, SpanSummary};
use crate::io::paf::EmbedFeatures; // if EmbedFeatures lives in paf.rs
use crate::io::paf::PairedLocus;
use crate::model::{ClassifyParams, Weights};

#[derive(Debug, Clone)]
#[allow(clippy::upper_case_acronyms)]
enum Call {
    NUMT,
    NIMT,
    Ambiguous,
}

fn clamp01(x: f32) -> f32 {
    x.clamp(0.0, 1.0)
}

// Map length (bp) to ~[0..1] with a soft saturation.
// Tunable: L50 = 5kb gives 0.5 contribution at 5kb, 0.8 at ~12kb, ~1 above ~25kb.
fn scale_len(len_bp: u32) -> f32 {
    let l = len_bp as f32;
    let l50 = 5_000.0f32;
    clamp01(l / (l + l50))
}

/// Classify each PairedLocus using alignment, coverage ratios, and spanning-read support.
/// Returns (pairs_tsv, classification_tsv).
///
/// Transparency columns added to pairs.tsv: rnuc, rmito, s_nuc, s_mito.
/// Drop them if you prefer the original minimal schema.
pub fn classify_pairs(
    pairs: &[PairedLocus],
    coverage: &CoverageSummary,
    spans: &SpanSummary,
    _embed: &[(String, EmbedFeatures)],
    w: Weights,
    params: ClassifyParams,
) -> Result<(String, String)> {
    // fast lookups
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
    // let embed_map: HashMap<&str, EmbedFeatures> =
    //     embed.iter().map(|(k, v)| (k.as_str(), v.clone())).collect();

    let dn_med = coverage.nuclear_median as f32;
    let dm_med = coverage.mito_median as f32;

    // Build TSV headers (with transparency columns)
    let mut pairs_tsv = String::from(
        "pair_id\tnuc_contig\tnuc_start\tnuc_end\tmito_contig\tmito_start\tmito_end\taln_len\taln_ident\trnuc\trmito\ts_nuc\ts_mito\tscore_numt\tscore_nimt\n"
    );
    let mut class_tsv = String::from("pair_id\tcall\tconfidence\treason_codes\n");

    for p in pairs {
        // Alignment terms
        let _a = clamp01(p.aln_ident); // already [0..1]
        let _l = scale_len(p.aln_len); // soft-saturated length

        // Depth ratios (local/median), guarded for zero medians
        let (d_nuc_loc, d_mito_loc) = depth_map
            .get(p.pair_id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));
        let rnuc = if dn_med > 0.0 {
            d_nuc_loc / dn_med
        } else {
            0.0
        };
        let rmito = if dm_med > 0.0 {
            d_mito_loc / dm_med
        } else {
            0.0
        };

        // Direction-specific depth terms
        // NUMT: locus on *nuclear* looks nuclear-like depth (rnuc ~ 1)
        // NIMT: locus on *mito*    looks mito-like depth    (rmito ~ 1)
        // Direction-specific depth terms around 1.0
        let d_numt = clamp01(1.0f32 - (rnuc - 1.0f32).abs());
        let d_nimt = clamp01(1.0f32 - (rmito - 1.0f32).abs());

        // Spanning fractions (already in [0,1])
        let (s_nuc, s_mito) = span_map
            .get(p.pair_id.as_str())
            .copied()
            .unwrap_or((0.0, 0.0));

        // --- NEW: depth contrast boost (log2 ratio + tanh gain) ---
        let eps = 1e-3_f32;
        let log2_ratio = ((rnuc + eps) / (rmito + eps)).ln() / std::f32::consts::LN_2; // log2
        let k_depth = 1.25f32; // gain: 0.8–2.0 is reasonable
        let depth_contrast = (k_depth * log2_ratio).tanh(); // in (-1, 1)
                                                            // positive => favors NUMT, negative => favors NIMT

        // --- NEW: explicit contrast on spans too ---
        let span_contrast = s_nuc - s_mito; // in [-1, 1]

        // Base alignment terms shared by both directions
        let a = clamp01(p.aln_ident);
        let l = scale_len(p.aln_len);

        // Common alignment contribution
        let base = w.w_a * a + w.w_l * l;

        // Old “pro” terms for each side
        let pro_numt = w.w_d * d_numt + w.w_s * s_nuc;
        let pro_nimt = w.w_d * d_nimt + w.w_s * s_mito;

        // --- NEW: add anti-evidence penalties & contrast boosts ---
        // Penalties: subtract the other side’s depth and spans.
        // Boosts: add depth_contrast (positive -> NUMT, negative -> NIMT)
        //         and span_contrast likewise.
        let pen_numt = w.w_d * d_nimt + w.w_s * s_mito;
        let pen_nimt = w.w_d * d_numt + w.w_s * s_nuc;

        // You can reuse w.w_d and w.w_s as the weights for the contrast boosts,
        // or introduce separate knobs. Here we reuse them (keeps CLI the same).
        let boost_numt = w.w_d * depth_contrast + w.w_s * span_contrast;
        let boost_nimt = -w.w_d * depth_contrast - w.w_s * span_contrast;

        // Final scores
        let score_numt = base + pro_numt - pen_numt + boost_numt; // + w.w_f * e.f_nuc_embed;
        let score_nimt = base + pro_nimt - pen_nimt + boost_nimt; // + w.w_f * e.f_mito_embed;

        // Decision
        let delta = (score_numt - score_nimt).abs();
        let call = if (score_numt - score_nimt) >= params.call_threshold {
            Call::NUMT
        } else if (score_nimt - score_numt) >= params.call_threshold {
            Call::NIMT
        } else {
            Call::Ambiguous
        };

        // Emit rows
        // pairs.tsv (add rnuc/rmito/s_* for transparency)
        use std::fmt::Write as _;
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

        // classification.tsv
        let (call_str, reason) = match call {
            Call::NUMT => ("Likely_NUMT", "score_difference"),
            Call::NIMT => ("Likely_NIMT", "score_difference"),
            Call::Ambiguous => ("Ambiguous", "delta_below_threshold"),
        };
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
