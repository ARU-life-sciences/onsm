use anyhow::{anyhow, Context, Result};
use paf::Reader as PafReader;
use serde::{Deserialize, Serialize};
use std::path::Path;

use crate::model::PairedLocus;

/// Thin, crate-internal PAF record (we compute identity here).
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
    pub identity: f32,
    pub strand: char,
}

impl From<paf::PafRecord> for PafRecord {
    fn from(r: paf::PafRecord) -> Self {
        let matches = r.residue_matches();
        let alnlen = r.alignment_block_len();
        let identity = if alnlen > 0 {
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
            identity,
            strand: r.strand(),
        }
    }
}

/// Read & filter PAF: keep records with identity ≥ min_id and length ≥ min_len.
pub fn read_paf(path: &Path, min_id: f32, min_len: u32) -> Result<Vec<PafRecord>> {
    if !path.exists() {
        return Err(anyhow!("PAF not found: {}", path.display()));
    }
    let mut reader =
        PafReader::from_path(path).with_context(|| format!("open PAF {}", path.display()))?;
    let mut out = Vec::new();
    for rec in reader.records() {
        let r = rec.with_context(|| format!("read PAF record in {}", path.display()))?;
        let pr: PafRecord = r.into();
        if pr.identity >= min_id && pr.alnlen >= min_len {
            out.push(pr);
        }
    }
    Ok(out)
}

/// Very simple pairing:
/// drive by mito→nuclear records, look for best reciprocal nuclear→mito by swapped names.
pub fn pair_and_merge(
    m2n: &[PafRecord],
    n2m: Vec<PafRecord>,
    _merge_gap: u32,
) -> Result<Vec<PairedLocus>> {
    let mut loci = Vec::new();
    for (i, rec) in m2n.iter().enumerate() {
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

#[cfg(test)]
mod tests {
    use super::*;
    use paf::{PafRecord as R, Tag, Type, Writer};
    use std::collections::HashMap;
    use tempfile::NamedTempFile;

    #[test]
    fn read_filter_roundtrip() {
        let tmp = NamedTempFile::new().unwrap();
        {
            let mut w = Writer::from_path(tmp.path()).unwrap();
            let mut opt = HashMap::new();
            opt.insert("tp".to_string(), Tag::tp(Type::Char('P')));
            // matches=90, alnlen=100 → 0.9
            let rec = R::new(
                "mito1".to_string(),
                1000,
                0,
                100,
                '+',
                "chr1".to_string(),
                5000,
                1000,
                1100,
                90,
                100,
                60,
                opt,
            );
            w.write_record(&rec).unwrap();
        }
        let v = read_paf(tmp.path(), 0.90, 50).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].qname, "mito1");
        assert_eq!(v[0].tname, "chr1");
    }
}
