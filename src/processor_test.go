package main

import (
	"testing"
)

// Helper function to create test CDSRow
func makeCDSRow(geneID, geneName, cdsID, chr string, chrStart, chrEnd, cdsStart, cdsEnd, length, strand int) CDSRow {
	return CDSRow{
		GeneID:         geneID,
		GeneName:       geneName,
		CDSID:          cdsID,
		Chr:            chr,
		ChrCodingStart: chrStart,
		ChrCodingEnd:   chrEnd,
		CDSStart:       cdsStart,
		CDSEnd:         cdsEnd,
		Length:         length,
		Strand:         strand,
	}
}

func TestNewCDSProcessor(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("GENE1", "GENE1", "CDS1", "chr1", 1000, 2000, 1, 1001, 1001, 1),
	}
	processor := NewCDSProcessor(rows)

	if len(processor.Rows()) != 1 {
		t.Errorf("Expected 1 row, got %d", len(processor.Rows()))
	}
	if processor.Rows()[0].GeneName != "GENE1" {
		t.Errorf("Expected GENE1, got %s", processor.Rows()[0].GeneName)
	}
}

func TestFilterExcludedChromosomes(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G3", "GENE3", "CDS3", "chr3", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).FilterExcludedChromosomes([]string{"chr2"})

	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows after filtering, got %d", len(processor.Rows()))
	}

	for _, r := range processor.Rows() {
		if r.Chr == "chr2" {
			t.Errorf("chr2 should have been filtered out")
		}
	}
}

func TestFilterExcludedChromosomesEmpty(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).FilterExcludedChromosomes([]string{})

	if len(processor.Rows()) != 1 {
		t.Errorf("Empty exclude list should not filter, got %d rows", len(processor.Rows()))
	}
}

func TestFilterOnlyChromosomes(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G3", "GENE3", "CDS3", "chr3", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).FilterOnlyChromosomes([]string{"chr1", "chr3"})

	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows, got %d", len(processor.Rows()))
	}

	for _, r := range processor.Rows() {
		if r.Chr == "chr2" {
			t.Errorf("chr2 should have been filtered out")
		}
	}
}

func TestFilterOnlyChromosomesEmpty(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).FilterOnlyChromosomes([]string{})

	if len(processor.Rows()) != 1 {
		t.Errorf("Empty only list should not filter, got %d rows", len(processor.Rows()))
	}
}

func TestNormalizeChromosomeNamesWithOverlap(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1),
	}
	genomeChrs := []string{"chr1", "chr2", "chr3"}
	processor := NewCDSProcessor(rows).NormalizeChromosomeNames(genomeChrs)

	// Should keep both since chr1 and chr2 overlap
	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows, got %d", len(processor.Rows()))
	}
}

func TestNormalizeChromosomeNamesWithoutOverlapAddPrefix(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "2", 100, 200, 1, 101, 101, 1),
	}
	genomeChrs := []string{"chr1", "chr2", "chr3"}
	processor := NewCDSProcessor(rows).NormalizeChromosomeNames(genomeChrs)

	// Should add "chr" prefix and keep rows
	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows, got %d", len(processor.Rows()))
	}
	if processor.Rows()[0].Chr != "chr1" {
		t.Errorf("Expected chr1, got %s", processor.Rows()[0].Chr)
	}
}

func TestRemoveIncompleteRecords(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1), // Missing GeneID
		makeCDSRow("G3", "", "CDS3", "chr3", 100, 200, 1, 101, 101, 1),    // Missing GeneName
		makeCDSRow("G4", "GENE4", "", "chr4", 100, 200, 1, 101, 101, 1),   // Missing CDSID
	}
	validChrs := []string{"chr1"}
	processor := NewCDSProcessor(rows).WithValidChromosomes(validChrs).RemoveIncompleteRecords()

	if len(processor.Rows()) != 1 {
		t.Errorf("Expected 1 complete record, got %d", len(processor.Rows()))
	}
	if processor.Rows()[0].GeneName != "GENE1" {
		t.Errorf("Expected GENE1, got %s", processor.Rows()[0].GeneName)
	}
}

func TestRemoveIncompleteRecordsInvalidChromosomes(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1),
	}
	validChrs := []string{"chr1"}
	processor := NewCDSProcessor(rows).WithValidChromosomes(validChrs).RemoveIncompleteRecords()

	if len(processor.Rows()) != 1 {
		t.Errorf("Expected 1 record with valid chromosome, got %d", len(processor.Rows()))
	}
}

func TestValidateCoordinates(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 0, 500, 1, 251, 251, 1), // ChrCodingStart is 0 (< 1) and ChrCodingEnd > length
	}
	contigs := []Contig{{Name: "chr1", Length: 250}}
	processor := NewCDSProcessor(rows).WithGenome(contigs)

	processor, err := processor.ValidateCoordinates()
	if err == nil {
		t.Errorf("Expected error for out-of-bounds coordinates")
	}
	if len(processor.Rows()) != 1 {
		t.Errorf("Expected 1 valid row after filtering out-of-bounds rows, got %d", len(processor.Rows()))
	}
}

func TestValidateCoordinatesNoGenome(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows)

	_, err := processor.ValidateCoordinates()
	if err == nil {
		t.Errorf("Expected error when genome lengths not set")
	}
}

func TestFindFullCDS(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1), // Full CDS: starts at 1, ends at 101
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 100, 200, 5, 101, 97, 1),  // Not full: doesn't start at 1
		makeCDSRow("G3", "GENE3", "CDS3", "chr1", 100, 200, 1, 99, 100, 1),  // Not full: ends at 99 but length is 100
	}
	processor := NewCDSProcessor(rows)
	fullCDS := processor.FindFullCDS()

	if len(fullCDS) != 1 {
		t.Errorf("Expected 1 full CDS, got %d", len(fullCDS))
	}
	if _, ok := fullCDS["CDS1"]; !ok {
		t.Errorf("Expected CDS1 to be marked as full")
	}
}

func TestTrimBoundaryBases(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 1, 200, 1, 200, 200, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 100, 200, 1, 101, 101, 1),
	}
	contigs := []Contig{{Name: "chr1", Length: 200}}
	processor := NewCDSProcessor(rows).WithGenome(contigs).TrimBoundaryBases()

	result := processor.Rows()
	// First row should have +3 to start (was 1) and -3 to end (was 200, chr length is 200)
	if result[0].ChrCodingStart != 4 {
		t.Errorf("Expected ChrCodingStart 4, got %d", result[0].ChrCodingStart)
	}
	if result[0].ChrCodingEnd != 197 {
		t.Errorf("Expected ChrCodingEnd 197, got %d", result[0].ChrCodingEnd)
	}
}

func TestDeduplicateByCDSKey(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1), // Duplicate
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 300, 400, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).DeduplicateByCDSKey()

	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 unique CDS after deduplication, got %d", len(processor.Rows()))
	}
}

func TestSortByGeneAndLength(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "ZEBRA", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "APPLE", "CDS2", "chr1", 100, 200, 1, 201, 201, 1),
		makeCDSRow("G3", "APPLE", "CDS3", "chr1", 100, 200, 1, 151, 151, 1),
	}
	processor := NewCDSProcessor(rows).SortByGeneAndLength()

	result := processor.Rows()
	// Should be sorted by gene name first (APPLE, APPLE, ZEBRA)
	// Then by length descending (APPLE 201, then 151)
	if result[0].GeneName != "APPLE" || result[0].Length != 201 {
		t.Errorf("Expected APPLE with length 201 first")
	}
	if result[1].GeneName != "APPLE" || result[1].Length != 151 {
		t.Errorf("Expected APPLE with length 151 second")
	}
	if result[2].GeneName != "ZEBRA" {
		t.Errorf("Expected ZEBRA last")
	}
}

func TestFilterByFullCDS(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G3", "GENE3", "CDS3", "chr1", 100, 200, 1, 101, 101, 1),
	}
	fullCDS := map[string]struct{}{"CDS1": {}, "CDS3": {}}
	processor := NewCDSProcessor(rows).FilterByFullCDS(fullCDS)

	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows in fullCDS, got %d", len(processor.Rows()))
	}
	for _, r := range processor.Rows() {
		if r.CDSID == "CDS2" {
			t.Errorf("CDS2 should have been filtered out")
		}
	}
}

func TestFilterValidCDSLength(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 99, 1),  // 99 % 3 = 0, valid
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 100, 200, 1, 102, 100, 1), // 100 % 3 != 0, invalid
		makeCDSRow("G3", "GENE3", "CDS3", "chr1", 100, 200, 1, 104, 102, 1), // 102 % 3 = 0, valid
	}
	processor := NewCDSProcessor(rows).FilterValidCDSLength()

	if len(processor.Rows()) != 2 {
		t.Errorf("Expected 2 rows with valid length, got %d", len(processor.Rows()))
	}
	for _, r := range processor.Rows() {
		if r.Length%3 != 0 {
			t.Errorf("Length %d is not divisible by 3", r.Length)
		}
	}
}

func TestSortByCoordinates(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr2", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 200, 300, 1, 101, 101, 1),
		makeCDSRow("G3", "GENE3", "CDS3", "chr1", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).SortByCoordinates()

	result := processor.Rows()
	// Should be sorted by chr, then by ChrCodingStart, then by CDSID
	if result[0].Chr != "chr1" || result[0].ChrCodingStart != 100 {
		t.Errorf("Expected chr1:100 first")
	}
	if result[1].Chr != "chr1" || result[1].ChrCodingStart != 200 {
		t.Errorf("Expected chr1:200 second")
	}
	if result[2].Chr != "chr2" {
		t.Errorf("Expected chr2 last")
	}
}

func TestChainedOperations(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 99, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 99, 1),
		makeCDSRow("", "GENE3", "CDS3", "chr1", 100, 200, 1, 101, 99, 1), // Missing GeneID
		makeCDSRow("G4", "GENE4", "CDS4", "chr3", 100, 200, 1, 101, 100, 1),
	}

	validChrs := []string{"chr1", "chr2"}
	processor := NewCDSProcessor(rows).
		FilterOnlyChromosomes(validChrs).
		WithValidChromosomes(validChrs).
		RemoveIncompleteRecords().
		FilterValidCDSLength().
		SortByCoordinates()

	result := processor.Rows()
	// Should have 2 records (CDS1 and CDS2 - CDS3 missing GeneID, CDS4 filtered out)
	if len(result) != 2 {
		t.Errorf("Expected 2 rows after chaining, got %d", len(result))
	}
}

func TestProcessorWithGenome(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
	}
	contigs := []Contig{{Name: "chr1", Seq: "ACGTACGTACGT", Length: 12}}

	processor := NewCDSProcessor(rows).WithGenome(contigs)

	if processor.genome["chr1"] != "ACGTACGTACGT" {
		t.Errorf("Genome not set correctly")
	}
	if processor.lengths["chr1"] != 12 {
		t.Errorf("Lengths not set correctly")
	}
}

func TestRows(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr1", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows)

	result := processor.Rows()
	if len(result) != 2 {
		t.Errorf("Expected 2 rows, got %d", len(result))
	}
}

func TestEdgeCaseEmptyInput(t *testing.T) {
	processor := NewCDSProcessor([]CDSRow{}).
		FilterExcludedChromosomes([]string{"chr1"}).
		FilterOnlyChromosomes([]string{"chr2"}).
		RemoveIncompleteRecords()

	if len(processor.Rows()) != 0 {
		t.Errorf("Expected 0 rows for empty input, got %d", len(processor.Rows()))
	}
}

func TestEdgeCaseAllFiltered(t *testing.T) {
	rows := []CDSRow{
		makeCDSRow("G1", "GENE1", "CDS1", "chr1", 100, 200, 1, 101, 101, 1),
		makeCDSRow("G2", "GENE2", "CDS2", "chr2", 100, 200, 1, 101, 101, 1),
	}
	processor := NewCDSProcessor(rows).
		FilterExcludedChromosomes([]string{"chr1", "chr2"})

	if len(processor.Rows()) != 0 {
		t.Errorf("Expected 0 rows after filtering all, got %d", len(processor.Rows()))
	}
}
