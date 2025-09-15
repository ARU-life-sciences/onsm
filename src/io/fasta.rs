//! FASTA validation with transparent .gz support via needletail.
//!
//! Accepts .fa/.fasta and .fa.gz/.fasta.gz. We just sanity-check that the file
//! exists and contains at least one sequence record parsable as FASTA-like
//! (needletail's parser handles gz internally).

use anyhow::{anyhow, Result};
use fs_err as fs;
use needletail::parse_fastx_file;
use std::path::Path;

pub fn validate_fasta(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(anyhow!("FASTA not found: {:?}", path));
    }

    // Quick guard against empty files (saves a confusing parser error).
    let meta = fs::metadata(path)?;
    if meta.len() == 0 {
        return Err(anyhow!("FASTA is empty: {:?}", path));
    }

    let mut reader =
        parse_fastx_file(path).map_err(|e| anyhow!("failed to open FASTA {:?}: {}", path, e))?;

    // Read the first record to verify it's a FASTA-like record (parser succeeds).
    let mut saw_record = false;
    while let Some(res) = reader.next() {
        let rec = res.map_err(|e| anyhow!("failed to read FASTA record in {:?}: {}", path, e))?;
        // A sequence of nonzero length is good enough for validation here.
        if !rec.seq().is_empty() {
            saw_record = true;
            break;
        }
    }

    if !saw_record {
        return Err(anyhow!("FASTA appears to contain no sequences: {:?}", path));
    }

    Ok(())
}
