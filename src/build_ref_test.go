package main

import (
	"bufio"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"testing"
)

type expectedGeneSummary struct {
	GeneID          string
	GeneName        string
	ProteinID       string
	CDSLength       int
	Chr             string
	Strand          int
	IntervalsCDS    []Interval
	IntervalsSplice []int
	SeqCDS          string
	SeqCDS1Up       string
	SeqCDS1Down     string
	SeqSplice       string
	SeqSplice1Up    string
	SeqSplice1Down  string
	L               [192][4]int
}

type expectedBuildRef struct {
	Order  []string
	ByGene map[string]expectedGeneSummary
}

func loadExpectedFromRDA(t *testing.T, rdaPath string) expectedBuildRef {
	t.Helper()

	if _, err := exec.LookPath("docker"); err != nil {
		t.Skip("docker not found; skipping RDA-backed buildref test")
	}
	if err := exec.Command("docker", "exec", "r-dev", "true").Run(); err != nil {
		t.Skip("r-dev container not available; skipping RDA-backed buildref test")
	}

	rScript := strings.Join([]string{
		"args <- commandArgs(trailingOnly=TRUE)",
		"load(args[1])",
		"cat(paste('N', length(RefCDS), sep='\t'), '\n', sep='')",
		"for (i in seq_along(RefCDS)) {",
		"  x <- RefCDS[[i]]",
		"  iv <- paste(apply(x$intervals_cds, 1, function(r) sprintf('%d-%d', r[1], r[2])), collapse=',')",
		"  sp <- paste(x$intervals_splice, collapse=',')",
		"  fields <- c(",
		"    'G',",
		"    x$gene_name,",
		"    x$gene_id,",
		"    x$protein_id,",
		"    as.character(x$CDS_length),",
		"    x$chr,",
		"    as.character(x$strand),",
		"    iv,",
		"    sp,",
		"    Biostrings::toString(x$seq_cds),",
		"    Biostrings::toString(x$seq_cds1up),",
		"    Biostrings::toString(x$seq_cds1down),",
		"    if (!is.null(x$seq_splice)) Biostrings::toString(x$seq_splice) else '',",
		"    if (!is.null(x$seq_splice1up)) Biostrings::toString(x$seq_splice1up) else '',",
		"    if (!is.null(x$seq_splice1down)) Biostrings::toString(x$seq_splice1down) else '',",
		"    paste(as.integer(x$L), collapse=',')",
		"  )",
		"  cat(paste(fields, collapse='\t'), '\n', sep='')",
		"}",
	}, "\n")

	containerPath := filepath.Join("/workspace", filepath.ToSlash(rdaPath))
	cmd := exec.Command("docker", "exec", "-i", "r-dev", "Rscript", "-", containerPath)
	cmd.Stdin = strings.NewReader(rScript)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("failed to load expected values from RDA: %v: %s", err, strings.TrimSpace(string(out)))
	}

	expected := expectedBuildRef{ByGene: map[string]expectedGeneSummary{}}
	scanner := bufio.NewScanner(strings.NewReader(string(out)))
	scanner.Buffer(make([]byte, 0, 1024*1024), 32*1024*1024)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "Loading required namespace:") {
			continue
		}
		parts := strings.Split(line, "\t")
		switch parts[0] {
		case "N":
			if len(parts) != 2 {
				t.Fatalf("invalid N line: %q", line)
			}
			n, err := strconv.Atoi(parts[1])
			if err != nil {
				t.Fatalf("invalid gene count in line %q: %v", line, err)
			}
			expected.Order = make([]string, 0, n)
		case "G":
			if len(parts) != 16 {
				t.Fatalf("invalid G line: %q", line)
			}
			gene := expectedGeneSummary{
				GeneName:        parts[1],
				GeneID:          parts[2],
				ProteinID:       parts[3],
				CDSLength:       mustAtoi(t, parts[4], line),
				Chr:             parts[5],
				Strand:          mustAtoi(t, parts[6], line),
				IntervalsCDS:    parseIntervals(t, parts[7]),
				IntervalsSplice: parseIntList(t, parts[8]),
				SeqCDS:          parts[9],
				SeqCDS1Up:       parts[10],
				SeqCDS1Down:     parts[11],
				SeqSplice:       parts[12],
				SeqSplice1Up:    parts[13],
				SeqSplice1Down:  parts[14],
				L:               parseLMatrix(t, parts[15]),
			}
			expected.Order = append(expected.Order, gene.GeneName)
			expected.ByGene[gene.GeneName] = gene
		default:
			t.Fatalf("unexpected oracle line: %q", line)
		}
	}
	if err := scanner.Err(); err != nil {
		t.Fatalf("failed parsing oracle output: %v", err)
	}
	if len(expected.Order) == 0 {
		t.Fatalf("oracle did not return any genes")
	}
	return expected
}

func mustAtoi(t *testing.T, raw, line string) int {
	t.Helper()
	n, err := strconv.Atoi(raw)
	if err != nil {
		t.Fatalf("invalid int %q in line %q: %v", raw, line, err)
	}
	return n
}

func parseIntervals(t *testing.T, raw string) []Interval {
	t.Helper()
	if raw == "" {
		return nil
	}
	parts := strings.Split(raw, ",")
	out := make([]Interval, 0, len(parts))
	for _, part := range parts {
		bounds := strings.Split(part, "-")
		if len(bounds) != 2 {
			t.Fatalf("invalid interval %q", part)
		}
		out = append(out, Interval{
			Start: mustAtoi(t, bounds[0], part),
			End:   mustAtoi(t, bounds[1], part),
		})
	}
	return out
}

func parseIntList(t *testing.T, raw string) []int {
	t.Helper()
	if raw == "" {
		return nil
	}
	parts := strings.Split(raw, ",")
	out := make([]int, 0, len(parts))
	for _, part := range parts {
		out = append(out, mustAtoi(t, part, raw))
	}
	return out
}

func parseLMatrix(t *testing.T, raw string) [192][4]int {
	t.Helper()
	var L [192][4]int
	parts := strings.Split(raw, ",")
	if len(parts) != 192*4 {
		t.Fatalf("invalid flattened L length: got %d, want %d", len(parts), 192*4)
	}
	for col := 0; col < 4; col++ {
		for row := 0; row < 192; row++ {
			L[row][col] = mustAtoi(t, parts[col*192+row], raw)
		}
	}
	return L
}

func lTotals(L [192][4]int) (syn, mis, non, splice int) {
	for i := 0; i < 192; i++ {
		syn += L[i][0]
		mis += L[i][1]
		non += L[i][2]
		splice += L[i][3]
	}
	return syn, mis, non, splice
}

func TestBuildRef_ExtDataMatchesRDAOracle(t *testing.T) {
	base := filepath.Join("..", "archive", "dndscv", "inst", "extdata")
	cds := filepath.Join(base, "BioMart_human_GRCh37_chr3_segment.txt")
	fasta := filepath.Join(base, "chr3_segment.fa")
	rda := filepath.Join(base, "refcds_example_chr3_segment.rda")

	expected := loadExpectedFromRDA(t, rda)

	got, err := buildRef(cds, fasta, nil, nil, false)
	if err != nil {
		t.Fatalf("buildRef returned error: %v", err)
	}

	if len(got.RefCDS) != len(expected.Order) {
		t.Fatalf("unexpected gene count: got %d, want %d", len(got.RefCDS), len(expected.Order))
	}

	gotOrder := make([]string, 0, len(got.RefCDS))
	for _, gene := range got.RefCDS {
		gotOrder = append(gotOrder, gene.GeneName)
	}
	if !reflect.DeepEqual(gotOrder, expected.Order) {
		t.Fatalf("gene order mismatch: got %v, want %v", gotOrder, expected.Order)
	}

	wantGRGenes := make([]GeneRange, 0)
	for idx, geneName := range expected.Order {
		gene := expected.ByGene[geneName]
		for _, iv := range gene.IntervalsCDS {
			wantGRGenes = append(wantGRGenes, GeneRange{Index: idx + 1, Chr: gene.Chr, Start: iv.Start, End: iv.End, Gene: gene.GeneName})
		}
		for _, pos := range gene.IntervalsSplice {
			wantGRGenes = append(wantGRGenes, GeneRange{Index: idx + 1, Chr: gene.Chr, Start: pos, End: pos, Gene: gene.GeneName})
		}
	}
	if !reflect.DeepEqual(got.GRGenes, wantGRGenes) {
		t.Fatalf("gr_genes mismatch: got %d entries, want %d", len(got.GRGenes), len(wantGRGenes))
	}

	for _, gene := range got.RefCDS {
		want := expected.ByGene[gene.GeneName]
		if gene.GeneID != want.GeneID {
			t.Fatalf("gene_id mismatch for %s: got %s, want %s", gene.GeneName, gene.GeneID, want.GeneID)
		}
		if gene.ProteinID != want.ProteinID {
			t.Fatalf("protein_id mismatch for %s: got %s, want %s", gene.GeneName, gene.ProteinID, want.ProteinID)
		}
		if gene.CDSLength != want.CDSLength {
			t.Fatalf("CDS_length mismatch for %s: got %d, want %d", gene.GeneName, gene.CDSLength, want.CDSLength)
		}
		if gene.Chr != want.Chr {
			t.Fatalf("chr mismatch for %s: got %s, want %s", gene.GeneName, gene.Chr, want.Chr)
		}
		if gene.Strand != want.Strand {
			t.Fatalf("strand mismatch for %s: got %d, want %d", gene.GeneName, gene.Strand, want.Strand)
		}
		if !reflect.DeepEqual(gene.IntervalsCDS, want.IntervalsCDS) {
			t.Fatalf("intervals_cds mismatch for %s: got %v, want %v", gene.GeneName, gene.IntervalsCDS, want.IntervalsCDS)
		}
		if !reflect.DeepEqual(gene.IntervalsSplice, want.IntervalsSplice) {
			t.Fatalf("intervals_splice mismatch for %s: got %v, want %v", gene.GeneName, gene.IntervalsSplice, want.IntervalsSplice)
		}
		if gene.SeqCDS != want.SeqCDS {
			t.Fatalf("seq_cds mismatch for %s", gene.GeneName)
		}
		if gene.SeqCDS1Up != want.SeqCDS1Up {
			t.Fatalf("seq_cds1up mismatch for %s", gene.GeneName)
		}
		if gene.SeqCDS1Down != want.SeqCDS1Down {
			t.Fatalf("seq_cds1down mismatch for %s", gene.GeneName)
		}
		if gene.SeqSplice != want.SeqSplice {
			t.Fatalf("seq_splice mismatch for %s", gene.GeneName)
		}
		if gene.SeqSplice1Up != want.SeqSplice1Up {
			t.Fatalf("seq_splice1up mismatch for %s", gene.GeneName)
		}
		if gene.SeqSplice1Down != want.SeqSplice1Down {
			t.Fatalf("seq_splice1down mismatch for %s", gene.GeneName)
		}
		if gene.L != want.L {
			syn, mis, non, splice := lTotals(gene.L)
			wsyn, wmis, wnon, wsplice := lTotals(want.L)
			t.Fatalf("L matrix mismatch for %s: got totals [%d %d %d %d], want [%d %d %d %d]", gene.GeneName, syn, mis, non, splice, wsyn, wmis, wnon, wsplice)
		}
	}
}

func TestBuildRef_OnlyChrFilter(t *testing.T) {
	base := filepath.Join("..", "archive", "dndscv", "inst", "extdata")
	cds := filepath.Join(base, "BioMart_human_GRCh37_chr3_segment.txt")
	fasta := filepath.Join(base, "chr3_segment.fa")

	got, err := buildRef(cds, fasta, nil, []string{"3"}, false)
	if err != nil {
		t.Fatalf("buildRef returned error with onlychrs filter: %v", err)
	}
	if len(got.RefCDS) == 0 {
		t.Fatalf("expected at least one gene after onlychrs filter")
	}
	for _, gene := range got.RefCDS {
		if gene.Chr != "3" {
			t.Fatalf("onlychrs filter failed, got gene on chr %s", gene.Chr)
		}
	}
}

func TestParseInt(t *testing.T) {
	t.Run("valid", func(t *testing.T) {
		got, err := parseInt(" 42 ")
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if got != 42 {
			t.Fatalf("parseInt mismatch: got %d, want 42", got)
		}
	})

	t.Run("missing", func(t *testing.T) {
		if _, err := parseInt("-"); err == nil {
			t.Fatalf("expected error for missing value")
		}
	})

	t.Run("invalid", func(t *testing.T) {
		if _, err := parseInt("abc"); err == nil {
			t.Fatalf("expected error for invalid integer")
		}
	})
}

func TestUniqueAndSliceHelpers(t *testing.T) {
	if got := unique([]string{"a", "b", "a", "c", "b"}); !reflect.DeepEqual(got, []string{"a", "b", "c"}) {
		t.Fatalf("unique mismatch for strings: got %v", got)
	}

	if got := unique([]int{1, 1, 2, 2, 3, 1, 1}); !reflect.DeepEqual(got, []int{1, 2, 3}) {
		t.Fatalf("unique mismatch for ints: got %v", got)
	}

	if got := addToAll([]int{1, 5, 10}, -2); !reflect.DeepEqual(got, []int{-1, 3, 8}) {
		t.Fatalf("addToAll mismatch: got %v", got)
	}

	if !containsBase("ACGTN", 'N') {
		t.Fatalf("containsBase should find present base")
	}
	if containsBase("ACGT", 'N') {
		t.Fatalf("containsBase should not find absent base")
	}
}

func TestComplementFunctions(t *testing.T) {
	if got := compBase('A'); got != 'T' {
		t.Fatalf("compBase mismatch for A: got %c", got)
	}
	if got := compBase('X'); got != 'N' {
		t.Fatalf("compBase mismatch for unknown base: got %c", got)
	}

	if got := complement("AcgTn"); got != "TGCAN" {
		t.Fatalf("complement mismatch: got %s", got)
	}

	if got := reverseComplement("ACGTN"); got != "NACGT" {
		t.Fatalf("reverseComplement mismatch: got %s", got)
	}
}

func TestTranslateDNA(t *testing.T) {
	got, err := translateDNA("ATGGAA")
	if err != nil {
		t.Fatalf("unexpected translateDNA error: %v", err)
	}
	if got != "ME" {
		t.Fatalf("translateDNA mismatch: got %s, want ME", got)
	}

	if _, err := translateDNA("ATGGA"); err == nil {
		t.Fatalf("expected length error for non-multiple-of-3 sequence")
	}

	if _, err := translateDNA("ATGNNN"); err == nil {
		t.Fatalf("expected unsupported codon error")
	}
}

func TestFastaLoadAndIndex(t *testing.T) {
	dir := t.TempDir()
	fastaPath := filepath.Join(dir, "test.fa")
	fasta := ">chrB desc\nACGT\nTT\n>chrA\nGGCC\n"
	if err := os.WriteFile(fastaPath, []byte(fasta), 0o644); err != nil {
		t.Fatalf("failed to write fasta: %v", err)
	}

	contigs, err := readFasta(fastaPath)
	if err != nil {
		t.Fatalf("readFasta returned error: %v", err)
	}
	if len(contigs) != 2 {
		t.Fatalf("unexpected contig count: got %d, want 2", len(contigs))
	}
	if contigs[0].Name != "chrB" || contigs[0].Seq != "ACGTTT" || contigs[0].Length != 6 {
		t.Fatalf("unexpected first contig data: %#v", contigs[0])
	}
	if contigs[1].Name != "chrA" || contigs[1].Seq != "GGCC" || contigs[1].Length != 4 {
		t.Fatalf("unexpected second contig data: %#v", contigs[1])
	}

	// When .fai is absent readIndex falls back to FASTA and returns the
	// os.ReadFile error alongside the chromosome list so callers know the
	// index was missing. We verify the list is correct and ignore that error.
	chrs, _ := readIndex(fastaPath)
	if !reflect.DeepEqual(chrs, []string{"chrA", "chrB"}) {
		t.Fatalf("unexpected chromosome list without .fai: got %v", chrs)
	}

	faiPath := fastaPath + ".fai"
	fai := "chrB\t6\t0\t6\t7\nchrA\t4\t0\t4\t5\n"
	if err := os.WriteFile(faiPath, []byte(fai), 0o644); err != nil {
		t.Fatalf("failed to write fai: %v", err)
	}

	chrs, err = readIndex(fastaPath)
	if err != nil {
		t.Fatalf("readIndex with .fai returned error: %v", err)
	}
	if !reflect.DeepEqual(chrs, []string{"chrB", "chrA"}) {
		t.Fatalf("unexpected chromosome list with .fai: got %v", chrs)
	}
}

func TestCDSAndSpliceSequenceHelpers(t *testing.T) {
	genome := map[string]string{"chr1": "ACGTACGTACGT"}

	cdsSeq, err := getCDSSeq(genome, "chr1", []Interval{{Start: 2, End: 4}, {Start: 6, End: 8}}, 1)
	if err != nil {
		t.Fatalf("getCDSSeq returned error: %v", err)
	}
	if cdsSeq != "CGTCGT" {
		t.Fatalf("getCDSSeq mismatch: got %s", cdsSeq)
	}

	cdsSeqNeg, err := getCDSSeq(genome, "chr1", []Interval{{Start: 2, End: 4}, {Start: 6, End: 8}}, -1)
	if err != nil {
		t.Fatalf("getCDSSeq (negative strand) returned error: %v", err)
	}
	if cdsSeqNeg != "ACGACG" {
		t.Fatalf("getCDSSeq negative strand mismatch: got %s", cdsSeqNeg)
	}

	ctxUp, err := getCDSContext(genome, "chr1", []Interval{{Start: 2, End: 4}}, 1, true)
	if err != nil {
		t.Fatalf("getCDSContext upstream returned error: %v", err)
	}
	if ctxUp != "ACG" {
		t.Fatalf("getCDSContext upstream mismatch: got %s", ctxUp)
	}

	ctxDownNeg, err := getCDSContext(genome, "chr1", []Interval{{Start: 2, End: 4}}, -1, false)
	if err != nil {
		t.Fatalf("getCDSContext downstream on negative strand returned error: %v", err)
	}
	if ctxDownNeg != "CGT" {
		t.Fatalf("getCDSContext downstream negative mismatch: got %s", ctxDownNeg)
	}

	plusCDS := []CDSRow{{ChrCodingStart: 10, ChrCodingEnd: 20, Strand: 1}, {ChrCodingStart: 30, ChrCodingEnd: 40, Strand: 1}}
	if got := getSpliceSites(plusCDS); !reflect.DeepEqual(got, []int{21, 22, 25, 28, 29}) {
		t.Fatalf("getSpliceSites (+) mismatch: got %v", got)
	}

	minusCDS := []CDSRow{{ChrCodingStart: 10, ChrCodingEnd: 20, Strand: -1}, {ChrCodingStart: 30, ChrCodingEnd: 40, Strand: -1}}
	if got := getSpliceSites(minusCDS); !reflect.DeepEqual(got, []int{21, 22, 25, 28, 29}) {
		t.Fatalf("getSpliceSites (-) mismatch: got %v", got)
	}

	spSeq, err := getSpliceSeq(genome, "chr1", []int{1, 2, 3}, 1)
	if err != nil {
		t.Fatalf("getSpliceSeq returned error: %v", err)
	}
	if spSeq != "ACG" {
		t.Fatalf("getSpliceSeq mismatch: got %s", spSeq)
	}

	spSeqNeg, err := getSpliceSeq(genome, "chr1", []int{1, 2, 3}, -1)
	if err != nil {
		t.Fatalf("getSpliceSeq negative returned error: %v", err)
	}
	if spSeqNeg != "TGC" {
		t.Fatalf("getSpliceSeq negative mismatch: got %s", spSeqNeg)
	}
}

func TestSplitFilterAndKeyHelpers(t *testing.T) {
	rows := []CDSRow{
		{GeneID: "g1", GeneName: "A", CDSID: "c1", Length: 9, CDSStart: 1, CDSEnd: 9},
		{GeneID: "g1", GeneName: "A", CDSID: "c1", Length: 9, CDSStart: 1, CDSEnd: 9},
		{GeneID: "g1", GeneName: "A", CDSID: "c2", Length: 10, CDSStart: 1, CDSEnd: 10},
		{GeneID: "g2", GeneName: "B", CDSID: "c3", Length: 12, CDSStart: 2, CDSEnd: 12},
	}

	unique := uniqueCDSKeyRows(rows)
	if len(unique) != 3 {
		t.Fatalf("uniqueCDSKeyRows mismatch: got %d rows", len(unique))
	}

	full := findFullCDS(rows)
	if _, ok := full["c1"]; !ok {
		t.Fatalf("findFullCDS should contain c1")
	}
	if _, ok := full["c2"]; !ok {
		t.Fatalf("findFullCDS should contain c2")
	}
	if _, ok := full["c3"]; ok {
		t.Fatalf("findFullCDS should not contain c3")
	}

	filteredTable := filterCDSTable(unique, full)
	if len(filteredTable) != 1 || filteredTable[0].CDSID != "c1" {
		t.Fatalf("filterCDSTable mismatch: got %v", filteredTable)
	}

	filteredRows := filterRowsByCDS(rows, full)
	if len(filteredRows) != 3 {
		t.Fatalf("filterRowsByCDS mismatch: got %d rows", len(filteredRows))
	}

	byCDS := splitByCDSID(rows)
	if len(byCDS["c1"]) != 2 || len(byCDS["c2"]) != 1 {
		t.Fatalf("splitByCDSID mismatch: got %v", byCDS)
	}

	byGene := splitByGeneName(rows)
	if len(byGene["A"]) != 3 || len(byGene["B"]) != 1 {
		t.Fatalf("splitByGeneName mismatch: got %v", byGene)
	}

	keys := sortedKeys(map[string]int{"z": 1, "a": 2, "m": 3})
	if !reflect.DeepEqual(keys, []string{"a", "m", "z"}) {
		t.Fatalf("sortedKeys mismatch: got %v", keys)
	}
}

func TestReadCDSTable(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "cds.tsv")
	content := strings.Join([]string{
		"gene.id\tgene.name\tcds.id\tchr\tchr.coding.start\tchr.coding.end\tcds.start\tcds.end\tlength\tstrand",
		"g1\tA\tc1\t1\t10\t20\t1\t11\t11\t1",
		"g2\tB\tc2\t1\t30\t40\t1\t11\t11\t-1",
		"bad\trow\twith\tfew\tcols",
	}, "\n")
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("failed to write cds table: %v", err)
	}

	rows, err := readCDS(path)
	if err != nil {
		t.Fatalf("readCDS returned error: %v", err)
	}
	if len(rows) != 2 {
		t.Fatalf("readCDS mismatch: got %d rows", len(rows))
	}
	if rows[1].Strand != -1 {
		t.Fatalf("readCDS strand parsing mismatch: got %d", rows[1].Strand)
	}
}

func TestInitCodonMetaAndBuildLMatrix(t *testing.T) {
	meta := initCodonMeta()
	if len(meta.trinucs) != 64 {
		t.Fatalf("initCodonMeta trinuc size mismatch: got %d", len(meta.trinucs))
	}
	if len(meta.subsIndex) != 192 {
		t.Fatalf("initCodonMeta substitution index size mismatch: got %d", len(meta.subsIndex))
	}

	// This synthetic codon is chosen so the expected impact totals are known a priori:
	// AAA -> one synonymous, seven missense, one nonsense over all single-base coding changes.
	entry := RefCDSEntry{
		SeqCDS:          "AAA",
		SeqCDS1Up:       "AAA",
		SeqCDS1Down:     "AAA",
		IntervalsSplice: []int{1},
		SeqSplice:       "A",
		SeqSplice1Up:    "C",
		SeqSplice1Down:  "G",
	}

	L, err := buildLMatrix(entry, meta)
	if err != nil {
		t.Fatalf("buildLMatrix returned error: %v", err)
	}
	syn, mis, non, splice := lTotals(L)
	if syn != 1 || mis != 7 || non != 1 || splice != 3 {
		t.Fatalf("unexpected L totals: got [%d %d %d %d], want [1 7 1 3]", syn, mis, non, splice)
	}

	badCDS := entry
	badCDS.SeqCDS1Down = "AA"
	if _, err := buildLMatrix(badCDS, meta); err == nil {
		t.Fatalf("expected error for CDS context length mismatch")
	}

	badSplice := entry
	badSplice.SeqSplice1Down = ""
	if _, err := buildLMatrix(badSplice, meta); err == nil {
		t.Fatalf("expected error for splice context length mismatch")
	}
}

func TestBuildRefResultMarshalJSON(t *testing.T) {
	r := BuildRefResult{
		RefCDS: []RefCDSEntry{
			{GeneName: "Z"},
			{GeneName: "A"},
			{GeneName: "M"},
		},
	}

	b, err := r.MarshalJSON()
	if err != nil {
		t.Fatalf("MarshalJSON returned error: %v", err)
	}

	var out map[string][]string
	if err := json.Unmarshal(b, &out); err != nil {
		t.Fatalf("json.Unmarshal returned error: %v", err)
	}

	if !reflect.DeepEqual(out["genes"], []string{"A", "M", "Z"}) {
		t.Fatalf("MarshalJSON genes mismatch: got %v", out["genes"])
	}
}
