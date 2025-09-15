# onsm

The Organelle Nuclear Similarity Mapper.

The idea is to attempt to classify NUMTs vs NIMTs with some certainty. NUMTs are nuclear insertions of mitochondrial DNA, and NIMTs are the reverse, nuclear insertions into the mitochondrial genome.

Some features so far:

- Detects candidate homologous loci between nuclear and mitochondrial assemblies.
- Scores each locus for NUMT vs NIMT signal using long read support.
- Reports classification (Likely_NUMT, Likely_NIMT, Ambiguous) with confidence scores.
- Provides genome-level summary metrics (percent nuclear genome NUMT, percent mitochondrial genome NIMT, etc.).
- Works with HiFi or ONT long reads.

## Install

Not on crates.io and no releases yet, so compile:

```bash
git clone https://github.com/YOURNAME/onsm.git
cd onsm
cargo build --release
```

## Dependencies

Requires:

- minimap2 >2.24
- samtools >1.16

## Usage

One shot classification:

```bash
onsm classify \
  --mito data/Arabidopsis_thaliana.fasta \
  --nuclear data/GCA_933208065.1.fasta.gz \
  --reads data/Arabidopsis_thaliana.ccs.fastq.gz \
  --platform hifi \
  --out results_dir
  --keep-tmp
```

If you kept temp files (as above) you can use `reuse` (mainly for dev):

```bash
onsm reuse \
  --mito mito.fasta \
  --nuclear nuclear.fasta.gz \
  --platform hifi \
  --paf-mito-to-nuc mito_to_nuc.paf \
  --paf-nuc-to-mito nuc_to_mito.paf \
  --bam-reads-to-nuc reads_to_nuc.bam \
  --bam-reads-to-mito reads_to_mito.bam \
  --out reuse_results
```

