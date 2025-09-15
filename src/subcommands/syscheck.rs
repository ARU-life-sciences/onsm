//! `onsm syscheck` — environment & versions.

use anyhow::Result;
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;
use sysinfo::System;

use crate::util::mapping;

#[derive(Args, Debug)]
pub struct CmdSyscheck {
    #[arg(long)]
    pub out: Option<PathBuf>,
    /// Optional explicit minimap2 binary
    #[arg(long, value_name = "PATH")]
    pub minimap2: Option<PathBuf>,
    /// Optional explicit samtools binary
    #[arg(long, value_name = "PATH")]
    pub samtools: Option<PathBuf>,
}

impl CmdSyscheck {
    pub fn run(self) -> Result<()> {
        let mut s = System::new_all();
        s.refresh_all();

        // Resolve binaries (allow CLI flags / env / PATH)
        let (mm2_bin, sam_bin) =
            mapping::resolve_bins(self.minimap2.as_deref(), self.samtools.as_deref())?;

        // Try to grab versions
        let mm2_version = mapping::get_version(&mm2_bin).unwrap_or_else(|e| format!("error: {e}"));
        let sam_version = mapping::get_version(&sam_bin).unwrap_or_else(|e| format!("error: {e}"));

        let obj = serde_json::json!({
            "onsm_version": env!("CARGO_PKG_VERSION"),
            "rustc": rustc_version_runtime::version().to_string(),
            "cpus": s.cpus().len(),
            "total_memory_mb": s.total_memory() / 1024 / 1024,
            "executables": {
                "minimap2": {
                    "path": mm2_bin,
                    "version": mm2_version,
                },
                "samtools": {
                    "path": sam_bin,
                    "version": sam_version,
                }
            },
        });

        if let Some(path) = self.out {
            serde_json::to_writer_pretty(fs::File::create(path)?, &obj)?;
        } else {
            println!("{}", serde_json::to_string_pretty(&obj)?);
        }
        Ok(())
    }
}

// Tiny helper to read rustc version at runtime
mod rustc_version_runtime {
    pub fn version() -> String {
        option_env!("RUSTC_VERSION")
            .unwrap_or("unknown")
            .to_string()
    }
}
