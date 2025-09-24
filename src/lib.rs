pub mod model;
pub mod scoring;
pub mod summary;

pub mod io {
    pub mod bam;
    pub mod fasta;
    pub mod paf;
    pub mod runfiles;
}

pub mod util {
    pub mod logging;
    pub mod mapping;
}

pub mod subcommands {
    pub mod classify;
    pub mod reuse;
}
