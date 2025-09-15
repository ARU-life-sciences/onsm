//! onsm: Organelle–Nuclear Similarity Mapper
//! Entry point only; see `cli` and `subcommands/*`.

mod cli;
mod io;
mod model;
mod scoring;
mod subcommands;
mod summary;
mod util;

use anyhow::Result;
use cli::Cli;

fn main() -> Result<()> {
    // Initialize env_logger with a default filter of "onsm=info"
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("onsm=info"))
        .format_timestamp_millis()
        .init();
    let cli = <Cli as clap::Parser>::parse();
    cli.run()
}
