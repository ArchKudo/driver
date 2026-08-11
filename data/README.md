# Data Directory

## hg19.db
Single-file SQLite database containing the dNdScv reference data for human GRCh37/hg19 (182MB).

### Tables

| Table | Rows | Description |
|---|---|---|
| `refcds` | 20,091 | One row per gene transcript — scalar fields + JSON text for nested data |
| `cds_intervals` | 194,927 | Normalized CDS exon intervals — replaces `gr_genes`, enables range queries |

### `refcds` columns
| Column | Type | Description |
|---|---|---|
| `gene_name` | TEXT PK | Gene symbol |
| `gene_id` | TEXT | Ensembl gene ID |
| `protein_id` | TEXT | Ensembl protein ID |
| `CDS_length` | INTEGER | Total CDS length in bp |
| `chr` | TEXT | Chromosome |
| `strand` | INTEGER | `+1` or `-1` |
| `intervals_splice` | TEXT (JSON) | Splice site positions `[int, ...]` |
| `seq_cds` | TEXT | CDS DNA sequence |
| `seq_cds1up` | TEXT | 1bp upstream context per CDS base |
| `seq_cds1down` | TEXT | 1bp downstream context per CDS base |
| `seq_splice` | TEXT | Splice site sequences |
| `seq_splice1up` | TEXT | Upstream context at splice sites |
| `seq_splice1down` | TEXT | Downstream context at splice sites |
| `L` | TEXT (JSON) | Mutation opportunity matrix `[[int×4]×192]` (syn/mis/non/spl) |

### `cds_intervals` columns
| Column | Type | Description |
|---|---|---|
| `rowid` | INTEGER PK | Auto |
| `gene_name` | TEXT FK | References `refcds(gene_name)` |
| `chr` | TEXT | Chromosome |
| `start` | INTEGER | Exon start (1-based) |
| `end` | INTEGER | Exon end |

### Indexes
| Name | On | Purpose |
|---|---|---|
| `idx_refcds_gene` | `refcds(gene_name)` | O(log n) gene detail fetch |
| `idx_cds_chr_start` | `cds_intervals(chr, start)` | Fast mutation→gene range scan |
| `idx_cds_gene` | `cds_intervals(gene_name)` | Fast gene→exon lookup |
| `cds_rtree` | virtual rtree on `(start, end)` | Interval overlap queries |

### Example queries
```sql
-- Fetch gene details
SELECT gene_name, chr, strand, CDS_length FROM refcds WHERE gene_name = 'TP53';

-- Map a mutation at chr17:7577120 to a gene
SELECT r.gene_name FROM cds_intervals ci
JOIN refcds r USING (gene_name)
WHERE ci.chr = '17' AND ci.start <= 7577120 AND ci.end >= 7577120;

-- All exons for a gene
SELECT chr, start, end FROM cds_intervals WHERE gene_name = 'TP53';
```

---

## simple_breast.csv
Somatic mutation calls from a breast cancer cohort (18,637 mutations).

| Field | Type | Description |
|---|---|---|
| `sampleID` | string | Sample identifier |
| `chr` | string | Chromosome |
| `pos` | int | Genomic position (1-based) |
| `ref` | string | Reference allele |
| `mut` | string | Mutant allele |

## hg19.json

Serialized dNdScv reference data for human GRCh37/hg19, containing two top-level keys.

### `RefCDS` array of 20,091 genes

| Field | Type | Description |
|---|---|---|
| `gene_name` | string | Gene symbol |
| `gene_id` | string | Ensembl gene ID |
| `protein_id` | string | Ensembl protein ID |
| `CDS_length` | int | Total CDS length in bp |
| `chr` | string | Chromosome |
| `strand` | int | `+1` or `-1` |
| `intervals_cds` | `[[start, end], ...]` | Exon genomic coordinates |
| `intervals_splice` | `int[]` | Splice site genomic positions |
| `seq_cds` | string | CDS DNA sequence |
| `seq_cds1up` | string | 1bp upstream context per CDS base |
| `seq_cds1down` | string | 1bp downstream context per CDS base |
| `seq_splice` | string | Splice site sequences |
| `seq_splice1up` | string | Upstream context at splice sites |
| `seq_splice1down` | string | Downstream context at splice sites |
| `L` | `[[int×4], ...]×192` | Mutation opportunity counts (192 trinucleotide contexts × 4 types: syn/mis/non/spl) |

### `gr_genes` array of 1,068,978 genomic ranges

| Field | Type | Description |
|---|---|---|
| `seqnames` | string | Chromosome |
| `start` | int | Range start (1-based) |
| `end` | int | Range end |
| `width` | int | Range width in bp |
| `strand` | string | `+` or `-` |
| `name` | string | Gene symbol |

### Schema

```json
{
  "RefCDS": [
    {
      "gene_name": "A1BG",
      "gene_id": "ENSG00000121410",
      "protein_id": "ENSP00000263100",
      "CDS_length": 1488,
      "chr": "19",
      "strand": -1,
      "intervals_cds": [[58858388, 58858395], [58858719, 58859006]],
      "intervals_splice": [58858396, 58858397, 58858714],
      "seq_cds": "ATGCGT...",
      "seq_cds1up": "NATGCG...",
      "seq_cds1down": "TGCGTA...",
      "seq_splice": "ATGC...",
      "seq_splice1up": "NATG...",
      "seq_splice1down": "TGCA...",
      "L": [[0, 7, 0, 0], [3, 4, 0, 0]]
    }
  ],
  "gr_genes": [
    {
      "seqnames": "19",
      "start": 58858387,
      "end": 58864805,
      "width": 6419,
      "strand": "-",
      "name": "A1BG"
    }
  ]
}
```
