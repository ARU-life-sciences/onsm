# onsm

The Organelle Nuclear Similarity Mapper.

The idea is to attempt to classify NUMTs vs NIMTs with some certainty. NUMTs are nuclear insertions of mitochondrial DNA, and NIMTs are the reverse, nuclear insertions into the mitochondrial genome.

Some features so far:

- Detects candidate homologous loci between nuclear and mitochondrial assemblies.
- Scores each locus for NUMT vs NIMT signal using long read support.
- Reports classification (Likely_NUMT, Likely_NIMT, Ambiguous) with confidence scores.
- Provides genome-level summary metrics (percent nuclear genome NUMT, percent mitochondrial genome NIMT, etc.).
- Works with HiFi or ONT long reads (I've not tested with ONT).

## Install

Not on crates.io and no releases yet, so compile from this repo:

```bash
git clone https://github.com/ARU-life-sciences/onsm
cd onsm
cargo build --release
```

## Dependencies

Requires:

- minimap2 >2.24
- samtools >1.16

These are super easy to install via:
- https://github.com/lh3/minimap2
- https://www.htslib.org/download/

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

You may also want to run `onsm syscheck` which will dump some info about your system and make sure you're ready to run. You'll either need `minimap2` and `samtools` in PATH, or specify them in the `onsm syscheck`. Otherwise you'll get an error.

## Outputs

When you run onsm classify or onsm reuse, three main result files are created in the output directory:

### `classification.tsv`

This is the main call table. Each row corresponds to one candidate locus (a mito↔nuclear alignment pair).
Columns:

- pair_id – unique identifier for the candidate locus.
- call – classification of the locus:
  - Likely_NUMT: locus looks like a mitochondrial sequence inserted into the nuclear genome.
  - Likely_NIMT: locus looks like a nuclear sequence inserted into the mitochondrial genome.
  - Ambiguous: insufficient evidence to decide.
- confidence – a scaled score difference between the NUMT and NIMT models (higher = more confident).
- reason_codes – why a call was made (e.g. score_difference, delta_below_threshold).

### `pairs.tsv`

Detailed per-locus statistics from alignments and read support.
Columns:

- pair_id – matches classification.tsv.
- nuc_contig / nuc_start / nuc_end – coordinates of the nuclear locus.
- mito_contig / mito_start / mito_end – coordinates of the mitochondrial locus.
- aln_len / aln_ident – alignment length and identity (fraction).
- rnuc / rmito – normalized read depths (coverage ratios).
- s_nuc / s_mito – span fractions: proportion of reads spanning the locus in nuclear vs. mito references.
- score_numt / score_nimt – composite scores used by the classifier for NUMT vs NIMT hypotheses.

```
P000004   OZ173161.1  0 43942   u104  0 43942   43942   1.0000  0.768   0.703   0.001   0.001   0.4747   0.2947
```

Here, the alignment covers ~44 kb with ~77% identity; nuclear and mito coverages are close, but the scoring leaned toward NUMT (score_numt > score_nimt).

### `summary.tsv`

A high-level overview across all pairs:

- n_pairs – number of candidate loci found.
- n_numt / n_nimt – number of loci called as NUMT / NIMT.
- nuclear_bp_total – size of the nuclear assembly.
- nuclear_bp_numt – number of nuclear bases overlapping called NUMTs.
- nuclear_pct_numt – % of nuclear genome spanned by NUMTs.
- mito_bp_total – size of the mitochondrial assembly.
- mito_bp_nimt – number of mitochondrial bases overlapping called NIMTs.
- mito_pct_nimt – % of mitochondrial genome spanned by NIMTs.

## How are calls made?

The classifier compares evidence from both the **assembly alignments** (mito to nuclear) and the **read support** (long reads mapped to each assembly). For each candidate locus (row in `pairs.tsv`), the following features are considered:

1. **Alignment features**  
   - `aln_len` and `aln_ident` are combined into an **alignment score** (long, high-identity alignments are weighted higher).  
   - Alignments are paired in both directions (mito→nuc and nuc→mito) to define a candidate locus.

2. **Coverage ratios**  
   - `rnuc` = coverage in nuclear locus ÷ median nuclear coverage.  
   - `rmito` = coverage in mitochondrial locus ÷ median mitochondrial coverage.  
   - Intuition:  
     - If the nuclear copy has depth consistent with the nuclear genome (`rnuc ≈ 1`), and the mitochondrial locus is depleted (`rmito << 1`), this supports a **NUMT** (mito → nuc).  
     - If the mitochondrial copy has depth consistent with mitochondria, and the nuclear copy is depleted, this supports a **NIMT** (nuc → mito).

3. **Span fractions**  
   - `s_nuc` and `s_mito` = fraction of reads spanning across the candidate locus in each reference.  
   - Guards against partial alignments or collapsed repeats — true insertions are well-spanned in the “host” genome, but not in the “donor” genome.

4. **Scoring model**  
   Each locus is scored under two hypotheses:

   - **NUMT score**  
     ```
     score_numt = w_a * aln_score
                + w_l * length_score
                + w_d * depth_term(rnuc, rmito)
                + w_s * span_term(s_nuc, s_mito)
     ```

   - **NIMT score**  
     ```
     score_nimt = w_a * aln_score
                + w_l * length_score
                + w_d * depth_term(rmito, rnuc)
                + w_s * span_term(s_mito, s_nuc)
     ```

   where  
   - `w_a, w_l, w_d, w_s` are weights (default: 0.25, 0.15, 0.25, 0.25; configurable via CLI).  
   - `depth_term` and `span_term` penalize deviations from expected coverage/span in the host genome and low values in the donor genome.

5. **Decision rule**  
   - Compute the difference:  
     ```
     Δ = score_numt – score_nimt
     ```
   - If `Δ >= call_threshold` → **Likely_NUMT**  
   - If `Δ <= –call_threshold` → **Likely_NIMT**  
   - Else → **Ambiguous**  
   - Default `call_threshold = 0.15`.  
   - A stricter cutoff (`highconf_threshold = 0.30`) highlights particularly confident calls.

6. **Confidence value**  
   - Reported in `classification.tsv` as  
     ```
     confidence = |Δ|
     ```  
   - The `reason_codes` column records *why* the call was made:  
     - `score_difference`: one hypothesis clearly scored higher.  
     - `delta_below_threshold`: both scores too close → ambiguous.

