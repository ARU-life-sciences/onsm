//! `onsm dump` — extract a single pair's context.

use anyhow::{anyhow, Result};
use clap::Args;
use fs_err as fs;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct CmdDump {
    #[arg(long)]
    pub run: PathBuf,
    #[arg(long)]
    pub pair: String,
    #[arg(long)]
    pub fasta_out: Option<PathBuf>,
    #[arg(long)]
    pub json: Option<PathBuf>,
}

impl CmdDump {
    pub fn run(self) -> Result<()> {
        let pairs_path = self.run.join("pairs.tsv");
        if !pairs_path.exists() {
            return Err(anyhow!("pairs.tsv not found under {:?}", self.run));
        }
        let txt = fs::read_to_string(&pairs_path)?;
        let header = txt.lines().next().unwrap_or_default().to_string();
        let rec = txt
            .lines()
            .skip(1)
            .find(|l| l.starts_with(&self.pair))
            .ok_or_else(|| anyhow!("pair_id {} not found", self.pair))?
            .to_string();

        if let Some(js) = self.json {
            let obj = serde_json::json!({ "header": header, "record": rec });
            serde_json::to_writer_pretty(fs::File::create(js)?, &obj)?;
        }
        if let Some(fa) = self.fasta_out {
            // Placeholder: you’ll emit flanked sequences here in v1.1
            fs::write(
                fa,
                format!(">TODO_extract_sequences_for_{}\nNNNN\n", self.pair),
            )?;
        }
        Ok(())
    }
}
