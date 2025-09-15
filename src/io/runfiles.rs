//! Misc small helpers for files & hashing.

use anyhow::{anyhow, Result};
use fs_err as fs;
use std::io::Read;
use std::path::Path;

pub fn ensure_exists(path: &Path) -> Result<()> {
    if !path.exists() {
        Err(anyhow!("Path not found: {:?}", path))
    } else {
        Ok(())
    }
}

pub fn md5_of_file(path: &Path) -> Result<String> {
    let mut f = fs::File::open(path)?;
    let mut ctx = md5::Context::new();
    let mut buf = [0u8; 64 * 1024];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        ctx.consume(&buf[..n]);
    }
    Ok(format!("{:x}", ctx.compute()))
}
