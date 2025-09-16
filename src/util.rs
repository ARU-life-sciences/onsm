//! Logging, external process helpers (minimap2 wrappers), misc utilities.

use anyhow::{anyhow, Result};
use std::path::Path;
use std::process::Command;

pub mod logging {
    use anyhow::Result;
    use std::path::Path;

    /// Initialize env_logger once. Idempotent: subsequent calls are no-ops.
    ///
    /// Sets a default RUST_LOG if none provided by the user.
    pub fn init_logging(_outdir: &Path) -> Result<()> {
        // If the user didn't set RUST_LOG, use a sensible default.
        if std::env::var_os("RUST_LOG").is_none() {
            // Only onsm logs at info by default; users can override.
            std::env::set_var("RUST_LOG", "onsm=info,info");
        }
        // Try to init; if already initialized elsewhere, ignore the error.
        let _ = env_logger::Builder::from_default_env()
            .format_timestamp_millis()
            .try_init();
        Ok(())
    }
}

// Shell out to minimap2 in a controlled way
pub mod mapping {
    use super::*;

    pub fn get_version(bin: &Path) -> Result<String> {
        let output = Command::new(bin)
            .arg("--version")
            .output()
            .map_err(|e| anyhow!("failed to spawn {} --version: {}", bin.display(), e))?;

        if !output.status.success() {
            return Err(anyhow!(
                "{} --version exited with non-zero status: {}",
                bin.display(),
                output.status
            ));
        }

        let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
        if !stdout.is_empty() {
            return Ok(stdout);
        }

        let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
        if !stderr.is_empty() {
            return Ok(stderr);
        }

        Err(anyhow!(
            "{} --version produced no output on stdout/stderr",
            bin.display()
        ))
    }

    /// Resolve minimap2 and samtools executables.
    /// Priority: CLI override > environment variable > PATH search
    pub fn resolve_bins(
        mm2_opt: Option<&Path>,
        sam_opt: Option<&Path>,
    ) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
        // minimap2
        let mm2 = if let Some(p) = mm2_opt {
            if p.exists() {
                p.to_path_buf()
            } else {
                return Err(anyhow!("minimap2 not found at {:?}", p));
            }
        } else if let Ok(envp) = std::env::var("ONSM_MINIMAP2") {
            std::path::PathBuf::from(envp)
        } else {
            which::which("minimap2").map_err(|_| {
                anyhow!("minimap2 not found (set --minimap2, ONSM_MINIMAP2, or PATH)")
            })?
        };

        // samtools
        let sam = if let Some(p) = sam_opt {
            if p.exists() {
                p.to_path_buf()
            } else {
                return Err(anyhow!("samtools not found at {:?}", p));
            }
        } else if let Ok(envp) = std::env::var("ONSM_SAMTOOLS") {
            std::path::PathBuf::from(envp)
        } else {
            which::which("samtools").map_err(|_| {
                anyhow!("samtools not found (set --samtools, ONSM_SAMTOOLS, or PATH)")
            })?
        };

        Ok((mm2, sam))
    }

    pub fn map_asm_to_asm(
        mm2_bin: &Path,
        query_fa: &Path,
        target_fa: &Path,
        out_paf: &Path,
        threads: usize,
    ) -> Result<()> {
        let status = Command::new(mm2_bin)
            .args(["-x", "asm20", "-t", &threads.to_string()])
            .arg(target_fa)
            .arg(query_fa)
            .arg("-c") // output CIGAR for potential later use
            .arg("-o")
            .arg(out_paf)
            .status()?;
        if !status.success() {
            return Err(anyhow!("minimap2 asm→asm failed"));
        }
        Ok(())
    }

    pub fn map_reads_to_ref(
        mm2_bin: &Path,
        samtools_bin: &Path,
        platform: &str,
        reads: &[std::path::PathBuf],
        ref_fa: &Path,
        out_bam: &Path,
        threads: usize,
    ) -> Result<()> {
        // Ensure parent dir exists
        if let Some(parent) = out_bam.parent() {
            fs_err::create_dir_all(parent)?;
        }

        let preset = match platform {
            "hifi" => "map-hifi",
            "ont" => "map-ont",
            _ => "map-hifi",
        };

        let reads_str = reads
            .iter()
            .map(|p| format!("'{}'", p.display()))
            .collect::<Vec<_>>()
            .join(" ");

        // Temp prefix for sorting lives alongside the output BAM
        let tmp_prefix = out_bam.with_extension("tmp");

        // Emit SAM (-a), sort to BAM, then index the BAM.
        let pipe = format!(
            "set -euo pipefail; \
         '{}' -x {} -a -t {} '{}' {} \
         | '{}' sort -@ {} -O BAM -T '{}' -o '{}' ; \
         '{}' index -@ {} '{}'",
            mm2_bin.display(),
            preset,
            threads,
            ref_fa.display(),
            reads_str,
            samtools_bin.display(),
            threads,
            tmp_prefix.display(),
            out_bam.display(),
            samtools_bin.display(),
            threads,
            out_bam.display()
        );

        // Run
        let status = Command::new("bash").args(["-c", &pipe]).status()?;
        if !status.success() {
            return Err(anyhow!(
                "reads→ref mapping failed (minimap2/samtools pipeline non-zero exit)"
            ));
        }

        // Sanity: ensure BAM (and BAI) exist
        if !out_bam.exists() {
            return Err(anyhow!(
                "reads→ref mapping failed: BAM not created at {}",
                out_bam.display()
            ));
        }
        let bai = out_bam.with_extension(format!(
            "{}.bai",
            out_bam
                .extension()
                .and_then(|e| e.to_str())
                .unwrap_or("bam")
        ));
        // samtools writes either out.bam.bai or out.bai depending on version; check both
        let bai_alt = out_bam.with_extension("bai");
        if !bai.exists() && !bai_alt.exists() {
            return Err(anyhow!(
                "reads→ref mapping failed: BAM index not created for {}",
                out_bam.display()
            ));
        }

        // Clean sort temp if any remains (best-effort)
        let _ = fs_err::remove_file(tmp_prefix);

        Ok(())
    }
}

pub mod hashing {
    // reserved for future helpers; using runfiles::md5_of_file for now
}
