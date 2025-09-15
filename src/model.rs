//! Shared data types used across subcommands.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Weights {
    pub w_a: f32, // identity
    pub w_l: f32, // length
    pub w_d: f32, // depth conformity
    pub w_s: f32, // spanning reads
    pub w_f: f32, // flank embedding
}

impl Default for Weights {
    fn default() -> Self {
        Self {
            w_a: 0.25,
            w_l: 0.15,
            w_d: 0.40,
            w_s: 0.20,
            w_f: 0.00,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ClassifyParams {
    pub call_threshold: f32,
    pub highconf_threshold: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub mito_md5: String,
    pub nuclear_md5: String,
    pub reads_md5: Vec<String>,
    pub platform: String,
    pub threads: usize,
    pub min_id: f32,
    pub min_len: u32,
    pub merge_gap: u32,
    pub flank: u32,
    pub win: u32,
    pub seed: u64,
}

impl RunManifest {
    pub fn from_inputs(
        mito: &Path,
        nuclear: &Path,
        reads: &[std::path::PathBuf],
        platform: &str,
        threads: usize,
        min_id: f32,
        min_len: u32,
        merge_gap: u32,
        flank: u32,
        win: u32,
        seed: u64,
    ) -> Result<Self> {
        use crate::io::runfiles::md5_of_file;
        Ok(Self {
            mito_md5: md5_of_file(mito)?,
            nuclear_md5: md5_of_file(nuclear)?,
            reads_md5: reads
                .iter()
                .map(|r| md5_of_file(r))
                .collect::<Result<_>>()?,
            platform: platform.to_string(),
            threads,
            min_id,
            min_len,
            merge_gap,
            flank,
            win,
            seed,
        })
    }
}
