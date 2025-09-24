use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "onsm", version, about = "Minimal NUMT/NIMT detector")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    Classify(onsm::subcommands::classify::CmdClassify),
    Reuse(onsm::subcommands::reuse::CmdReuse),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Classify(cmd) => cmd.run(),
        Cmd::Reuse(cmd) => cmd.run(),
    }
}
