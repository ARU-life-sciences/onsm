use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Default algorithm constants (few knobs, sensible defaults)
pub const MIN_ID: f32 = 0.90;
pub const MIN_LEN: u32 = 100;
pub const MERGE_GAP: u32 = 50;
pub const FLANK_BP: u32 = 500; // window half-width
pub const WIN_BP: u32 = 250; // “spanning” sub-window half-width
pub const CALL_THRESHOLD: f32 = 0.15;
pub const HIGHCONF_THRESHOLD: f32 = 0.30;

// Scoring weights
pub const W_A: f32 = 0.25; // alignment identity
pub const W_L: f32 = 0.15; // alignment length (soft-saturated)
pub const W_D: f32 = 0.25; // depth consistency
pub const W_S: f32 = 0.25; // spanning support

/// A paired locus after reciprocal mapping/merging.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedLocus {
    pub pair_id: String,
    pub nuc_contig: String,
    pub nuc_start: u32,
    pub nuc_end: u32,
    pub mito_contig: String,
    pub mito_start: u32,
    pub mito_end: u32,
    pub aln_len: u32,
    pub aln_ident: f32, // [0,1]
}

/// Depth/coverage summary.
/// `per_pair[pid] = (nuc_local_median_depth, mito_local_median_depth)`
/// Medians are *absolute* here; scoring will normalize by the genome-wide medians below.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageSummary {
    pub nuclear_median: f64,
    pub mito_median: f64,
    pub per_pair: HashMap<String, (f32, f32)>,
}

/// Spanning-read support summary.
/// `per_pair[pid] = (frac_spanning_nuc_window, frac_spanning_mito_window)` in [0,1].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpanSummary {
    pub per_pair: HashMap<String, (f32, f32)>,
}

/// Immutable scoring params (constants exposed here).
#[derive(Debug, Clone, Copy)]
pub struct ClassifyParams {
    pub call_threshold: f32,
    pub highconf_threshold: f32,
}

/// Weights (pulled from constants)
#[derive(Debug, Clone, Copy)]
pub struct Weights {
    pub w_a: f32,
    pub w_l: f32,
    pub w_d: f32,
    pub w_s: f32,
}

impl Default for Weights {
    fn default() -> Self {
        Self {
            w_a: crate::model::W_A,
            w_l: crate::model::W_L,
            w_d: crate::model::W_D,
            w_s: crate::model::W_S,
        }
    }
}

impl Default for ClassifyParams {
    fn default() -> Self {
        Self {
            call_threshold: CALL_THRESHOLD,
            highconf_threshold: HIGHCONF_THRESHOLD,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub mito: PathBuf,
    pub nuclear: PathBuf,
    pub reads: Vec<PathBuf>,
    pub platform: String,
    pub threads: usize,

    // fixed thresholds used by your simplified pipeline
    pub min_id: f32,
    pub min_len: u32,
    pub merge_gap: u32,
    pub flank_bp: u32,
    pub win_bp: u32,
}

impl RunManifest {
    pub fn new(
        mito: &Path,
        nuclear: &Path,
        reads: &[PathBuf],
        platform: &str,
        threads: usize,
        min_id: f32,
        min_len: u32,
        merge_gap: u32,
        flank_bp: u32,
        win_bp: u32,
    ) -> Self {
        Self {
            mito: mito.to_path_buf(),
            nuclear: nuclear.to_path_buf(),
            reads: reads.to_vec(),
            platform: platform.to_string(),
            threads,
            min_id,
            min_len,
            merge_gap,
            flank_bp,
            win_bp,
        }
    }

    pub fn save_to(out_dir: &Path, m: &Self) -> anyhow::Result<()> {
        fs_err::create_dir_all(out_dir)?;
        let f = fs_err::File::create(out_dir.join("run_manifest.json"))?;
        serde_json::to_writer_pretty(f, m)?;
        Ok(())
    }

    pub fn load_from(out_dir: &Path) -> anyhow::Result<Self> {
        let f = fs_err::File::open(out_dir.join("run_manifest.json"))?;
        let m: Self = serde_json::from_reader(f)?;
        Ok(m)
    }
}
