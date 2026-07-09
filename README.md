
## Investigations

- Has site-based
- Has overall dn/ds

### [dndscv.R](archive/dndscv/R/dndscv.R)

- Get first 5 cols of the vcf
    - sampleID, chr, pos, ref, alt
- Remove pos where ref/alt are same
- Remove N/A
- If gene_list is CDKN2A replace it with the isoforms "CDKN2A.p14arf","CDKN2A.p16INK4a"
- Load reference
- Filter by gene_list are all genes if not specified
