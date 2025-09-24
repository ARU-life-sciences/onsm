use anyhow::Result;
use std::path::Path;

pub fn ensure_exists(p: &Path) -> Result<()> {
    if !p.exists() {
        return Err(anyhow::anyhow!("input not found: {}", p.display()));
    }
    if !p.is_file() {
        return Err(anyhow::anyhow!("not a file: {}", p.display()));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn exists_ok() {
        let f = NamedTempFile::new().unwrap();
        ensure_exists(f.path()).unwrap();
    }

    #[test]
    fn missing_err() {
        let e = ensure_exists(Path::new("/nope/nope/nope")).unwrap_err();
        assert!(e.to_string().contains("input not found"));
    }
}
