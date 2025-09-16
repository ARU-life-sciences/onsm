//! CLI definition and top-level dispatch.

use anyhow::Result;
use clap::{Parser, Subcommand};

use crate::subcommands::{classify::CmdClassify, reuse::CmdReuse, syscheck::CmdSyscheck};

#[derive(Parser, Debug)]
#[command(
    name = "onsm",
    version,
    about = "Organelle–Nuclear Similarity Mapper (long reads only)"
)]
pub struct Cli {
    #[command(subcommand)]
    cmd: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// One-shot pipeline: detect & score NUMTs vs NIMTs
    Classify(CmdClassify),

    /// Reuse existing PAF/BAM outputs to rescore without remapping
    Reuse(CmdReuse),

    /// Check environment, external tools, and features
    Syscheck(CmdSyscheck),
}

impl Cli {
    pub fn run(self) -> Result<()> {
        match self.cmd {
            Commands::Classify(cmd) => cmd.run(),
            Commands::Reuse(cmd) => cmd.run(),
            Commands::Syscheck(cmd) => cmd.run(),
        }
    }
}
