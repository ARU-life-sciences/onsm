use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::collections::HashMap;
use std::path::Path;

/// Checks file exists, is readable, and looks like FASTA by reading first record.
pub fn validate_fasta(p: &Path) -> Result<()> {
    if !p.exists() {
        return Err(anyhow::anyhow!("FASTA not found: {}", p.display()));
    }
    let mut rdr = parse_fastx_file(p).with_context(|| format!("open fasta {}", p.display()))?;
    let _ = rdr
        .next()
        .transpose()
        .with_context(|| format!("read first record in {}", p.display()))?;
    Ok(())
}

/// Returns map of contig name -> length (bp).
pub fn contig_lengths(p: &Path) -> Result<HashMap<String, u64>> {
    let mut m = HashMap::new();
    let mut rdr = parse_fastx_file(p).with_context(|| format!("open fasta {}", p.display()))?;
    while let Some(rec) = rdr
        .next()
        .transpose()
        .with_context(|| format!("read fasta {}", p.display()))?
    {
        let id = String::from_utf8_lossy(rec.id()).to_string();
        let len = rec.seq().len() as u64;
        m.insert(id, len);
    }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn lengths_ok() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, ">c1\nAAAAAA\n>c2\nACGTACGTAC\n>c3\nA\n").unwrap();
        let m = contig_lengths(f.path()).unwrap();
        assert_eq!(m.get("c1"), Some(&6));
        assert_eq!(m.get("c2"), Some(&10));
        assert_eq!(m.get("c3"), Some(&1));
    }
}
