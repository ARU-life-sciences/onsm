//! `onsm prep` — build & cache minimap2 indexes.

use crate::util::mapping;
use anyhow::Result;
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct CmdPrep {
    #[arg(long)]
    pub mito: PathBuf,
    #[arg(long)]
    pub nuclear: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long)]
    pub minimap2: Option<PathBuf>,
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    #[arg(long)]
    pub force: bool,
}

impl CmdPrep {
    pub fn run(self) -> Result<()> {
        fs::create_dir_all(&self.out)?;
        let mm2_bin = mapping::resolve_minimap2(self.minimap2.as_deref())?;
        log::info!("Using minimap2 at {mm2_bin:?}");
        if let Ok(v) = mapping::get_version(&mm2_bin) {
            log::info!("minimap2: {v}");
        }

        mapping::build_index(
            &mm2_bin,
            &self.mito,
            &self.out.join("mito.mmi"),
            self.force,
            self.threads,
        )?;
        mapping::build_index(
            &mm2_bin,
            &self.nuclear,
            &self.out.join("nuclear.mmi"),
            self.force,
            self.threads,
        )?;

        Ok(())
    }
}
