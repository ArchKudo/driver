
- Answers the question are specific genes accumulating more protein-altering mutations?
- Consequently we only care about the cds region which can encode into proteins and non-coding regions are ignored.

## Usage

### Input

Sample sheet (Required)
 - Should be a valid CSV
 - Should contain 5 columns, will fail otherwise
 - First row should be header
 - CSV columns are position based even though we take a header
 - Follows the order
  - sampleID, chr, pos, ref, alt


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

### [buildRef.R](archive/dndscv/R/buildref.R)

- First roadblock
    - Cannot directly convert provided hg38.rda to csv
    - Provides genomic ranges, refCDS

- Input parameters
    - cdsfile: transcript table?
    - genomefile: indexed reference?
    - ~~out~~
    - ~~numcode: Translate genes to protein the standard way always can ignore for now~~
    - excludechrs/onlychrs
    - ~~useids: Combine gene id/name also ignore~~


- Covariates business can be avoided without much loss of accuracy

- investing in rebuilding buildRef
    - Would require downloading from biomart
    - Would require translating gene to protein to check type of variation
    - can be a project on it's own

- just support hg19 reference for now, save it as go binary object instead of csv
- always use long names

- Get the table
    - "gene.id",
    - "gene.name",
    - "cds.id",
    - "chr",
    - "chr.coding.start",
    - "chr.coding.end",
    - "cds.start",
    - "cds.end",
    - "length",
    - "strand"

- Why is cds.start / end scalar and not array representing multiple exon start stops?



- Convert chr.coding.start..strand to numeric

- Check for empty gene names

- Replace by geneid:genename if longid

- Check for all unique genes

- build fasta index
- check for chromsomes in fasta file

- filter on chrs present in both fasta and cds table

- throw genes away not in chrs/contigs


## [driver.go](driver.go)


### [sample_sheet_test.go](sample_sheet_test.go)

- Contains unittests for loading samplesheet
