use anyhow::{anyhow, Context, Result};
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn resolve_bins(
    minimap2: Option<&Path>,
    samtools: Option<&Path>,
) -> Result<(PathBuf, PathBuf)> {
    let mm2 = match minimap2 {
        Some(p) => p.to_path_buf(),
        None => which::which("minimap2")
            .context("minimap2 not found in PATH. Install or pass --minimap2")?,
    };
    let sam = match samtools {
        Some(p) => p.to_path_buf(),
        None => which::which("samtools")
            .context("samtools not found in PATH. Install or pass --samtools")?,
    };
    Ok((mm2, sam))
}

pub fn get_version(bin: &Path) -> Result<String> {
    let out = Command::new(bin)
        .arg("--version")
        .output()
        .with_context(|| format!("spawn {} --version", bin.display()))?;
    let s = String::from_utf8_lossy(&out.stdout);
    Ok(s.lines().next().unwrap_or_default().to_string())
}

/// Run minimap2 assembly→assembly mapping with preset `-x asm10` to PAF.
pub fn map_asm_to_asm(
    mm2: &Path,
    query_fa: &Path,
    target_fa: &Path,
    out_paf: &Path,
    threads: usize,
) -> Result<()> {
    log::info!(
        "minimap2 asm-asm: {} → {} → {}",
        query_fa.display(),
        target_fa.display(),
        out_paf.display()
    );
    let status = Command::new(mm2)
        .args(["-x", "asm10", "-c", "-t"])
        .arg(threads.to_string())
        .arg(target_fa)
        .arg(query_fa)
        .arg("-o")
        .arg(out_paf)
        .status()
        .context("failed to spawn minimap2 for asm-asm")?;
    if !status.success() {
        return Err(anyhow!("minimap2 (asm-asm) failed with status {}", status));
    }
    Ok(())
}

/// Map reads→reference, convert to sorted BAM + index.
/// Presets: `map-hifi` or `map-ont`.
pub fn map_reads_to_ref(
    mm2: &Path,
    sam: &Path,
    platform: &str,
    reads: &[PathBuf],
    reference: &Path,
    out_bam: &Path,
    threads: usize,
) -> Result<()> {
    let preset = match platform {
        "hifi" => "map-hifi",
        "ont" => "map-ont",
        other => return Err(anyhow!("unknown --platform {other}; use hifi|ont")),
    };
    log::info!(
        "minimap2 reads→{} ({preset}): {} reads → {}",
        reference.display(),
        reads.len(),
        out_bam.display()
    );

    // minimap2 -x PRESET -a -t N ref.fa reads... | samtools sort -o out.bam
    let mut mm2_cmd = Command::new(mm2);
    mm2_cmd
        .args(["-x", preset, "-a", "-t"])
        .arg(threads.to_string());
    mm2_cmd.arg(reference);
    for r in reads {
        mm2_cmd.arg(r);
    }

    // Pipe to samtools sort
    let mut sort_cmd = Command::new(sam);
    sort_cmd.args(["sort", "-o"]).arg(out_bam);

    // Spawn with pipe
    let mut mm2_child = mm2_cmd
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("spawn minimap2 for reads→ref")?;
    let mm2_out = mm2_child.stdout.take().unwrap();
    let sort_status = sort_cmd
        .stdin(mm2_out)
        .status()
        .context("spawn samtools sort")?;
    let mm2_status = mm2_child.wait().context("wait minimap2")?;

    if !mm2_status.success() || !sort_status.success() {
        return Err(anyhow!(
            "reads→ref pipeline failed (minimap2={mm2_status}, sort={sort_status})"
        ));
    }

    // Index
    let status = Command::new(sam)
        .args(["index", out_bam.to_str().unwrap()])
        .status()
        .context("samtools index")?;
    if !status.success() {
        return Err(anyhow!("samtools index failed with {status}"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_bins_errors_when_missing() {
        // We can't guarantee PATH here; just ensure error messages are informative by
        // calling with obviously bad paths.
        let err = resolve_bins(Some(Path::new("/definitely/not/here")), None).unwrap_err();
        assert!(
            err.to_string().contains("minimap2"),
            "message mentions minimap2"
        );
    }

    #[test]
    fn preset_selection() {
        // no actual spawn, just exercise error branch
        let e = map_reads_to_ref(
            Path::new("minimap2"),
            Path::new("samtools"),
            "bad",
            &[],
            Path::new("ref.fa"),
            Path::new("out.bam"),
            1,
        )
        .unwrap_err();
        assert!(e.to_string().contains("unknown --platform"));
    }
}
