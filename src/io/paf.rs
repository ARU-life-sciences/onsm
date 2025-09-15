//! PAF I/O using the `paf` crate (v0.2.1) + pairing/merging helpers.
//
// Public API kept identical to the previous scaffold so other modules don’t change.

use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Thin, crate-internal alignment record used by downstream code.
/// Constructed from `paf::PafRecord`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PafRecord {
    pub qname: String,
    pub qstart: u32,
    pub qend: u32,
    pub tname: String,
    pub tstart: u32,
    pub tend: u32,
    pub matches: u32,
    pub alnlen: u32,
    pub mapq: u8,
    pub identity: f32, // computed as residue_matches / alignment_block_len
    pub strand: char,  // '+' or '-'
}

impl From<paf::PafRecord> for PafRecord {
    fn from(r: paf::PafRecord) -> Self {
        let matches = r.residue_matches();
        let alnlen = r.alignment_block_len();
        let ident = if alnlen > 0 {
            matches as f32 / alnlen as f32
        } else {
            0.0
        };
        Self {
            qname: r.query_name().to_string(),
            qstart: r.query_start(),
            qend: r.query_end(),
            tname: r.target_name().to_string(),
            tstart: r.target_start(),
            tend: r.target_end(),
            matches,
            alnlen,
            mapq: r.mapping_quality(),
            identity: ident,
            strand: r.strand(),
        }
    }
}

/// Read & filter a PAF file using `paf::Reader`.
pub fn read_paf(path: &Path, min_id: f32, min_len: u32) -> Result<Vec<PafRecord>> {
    if !path.exists() {
        return Err(anyhow!("PAF not found: {:?}", path));
    }
    let mut reader = paf::Reader::from_path(path)
        .map_err(|e| anyhow!("Failed to open PAF {:?}: {}", path, e))?;
    let mut out = Vec::new();
    for rec in reader.records() {
        let r = rec.map_err(|e| anyhow!("Error reading PAF record: {}", e))?;
        let ours: PafRecord = r.into();
        if ours.identity >= min_id && ours.alnlen >= min_len {
            out.push(ours);
        }
    }
    Ok(out)
}

/// Paired/merged locus (still minimal for v1).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedLocus {
    pub pair_id: String,
    pub nuc_contig: String,
    pub nuc_start: u32,
    pub nuc_end: u32,
    pub mito_contig: String,
    pub mito_start: u32,
    pub mito_end: u32,
    pub aln_len: u32,
    pub aln_ident: f32,
}

/// Very simple v1 pairing: use mito→nuclear hits as the driver; try to find a reciprocal
/// nuclear→mito hit with swapped contig names; pick the best by identity.
/// (You can enrich this later with chaining/merging and gap handling.)
pub fn pair_and_merge(
    m2n: Vec<PafRecord>,
    n2m: Vec<PafRecord>,
    _merge_gap: u32,
) -> Result<Vec<PairedLocus>> {
    let mut loci = Vec::new();
    for (i, rec) in m2n.into_iter().enumerate() {
        let best = n2m
            .iter()
            .filter(|r| r.qname == rec.tname && r.tname == rec.qname)
            .max_by(|a, b| a.identity.partial_cmp(&b.identity).unwrap());

        let ident = best
            .map(|b| b.identity.max(rec.identity))
            .unwrap_or(rec.identity);
        let (mito_s, mito_e) = (rec.qstart.min(rec.qend), rec.qstart.max(rec.qend));
        let (nuc_s, nuc_e) = (rec.tstart.min(rec.tend), rec.tstart.max(rec.tend));

        loci.push(PairedLocus {
            pair_id: format!("P{:06}", i + 1),
            nuc_contig: rec.tname.clone(),
            nuc_start: nuc_s,
            nuc_end: nuc_e,
            mito_contig: rec.qname.clone(),
            mito_start: mito_s,
            mito_end: mito_e,
            aln_len: rec.alnlen,
            aln_ident: ident,
        });
    }
    Ok(loci)
}

/// Embedding features placeholder (to be filled by short flank alignments later).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EmbedFeatures {
    pub f_nuc_embed: f32,
    pub f_mito_embed: f32,
}

pub fn compute_embedding_identities(
    pairs: &[PairedLocus],
    _mito: &std::path::Path,
    _nuclear: &std::path::Path,
    _flank: u32,
    _threads: usize,
) -> Result<Vec<(String, EmbedFeatures)>> {
    Ok(pairs
        .iter()
        .map(|p| (p.pair_id.clone(), EmbedFeatures::default()))
        .collect())
}

pub fn default_embed_stub(pairs: &[PairedLocus]) -> Vec<(String, EmbedFeatures)> {
    pairs
        .iter()
        .map(|p| (p.pair_id.clone(), EmbedFeatures::default()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use paf::{PafRecord as RPafRecord, Tag, Type, Writer};
    use std::collections::HashMap;
    use tempfile::NamedTempFile;

    #[test]
    fn roundtrip_read_with_paf_crate() {
        // Write a tiny PAF using the crate’s Writer, then read via our read_paf
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut w = Writer::from_path(tmp.path()).unwrap();
            let mut opt = HashMap::new();
            opt.insert("tp".to_string(), Tag::tp(Type::Char('P')));
            // 390/500 = 0.78 identity
            let rec = RPafRecord::new(
                "q1".to_string(),
                1000,
                10,
                510,
                '+',
                "t1".to_string(),
                2000,
                100,
                600,
                390,
                500,
                60,
                opt,
            );
            w.write_record(&rec).unwrap();
        }

        let v = read_paf(tmp.path(), 0.75, 100).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].qname, "q1");
        assert_eq!(v[0].tname, "t1");
        assert!((v[0].identity - 0.78).abs() < 1e-6);

        // Apply a stricter filter that drops it
        let v2 = read_paf(tmp.path(), 0.80, 100).unwrap();
        assert!(v2.is_empty());
    }
}
