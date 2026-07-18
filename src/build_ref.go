package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"os"
	"slices"
	"sort"
	"strconv"
	"strings"
)

type CDSRow struct {
	GeneID         string
	GeneName       string
	CDSID          string
	Chr            string
	ChrCodingStart int
	ChrCodingEnd   int
	CDSStart       int
	CDSEnd         int
	Length         int
	Strand         int
}

// Fetched from BioMart with CDSRow fields
func readCDS(cdsfile string) ([]CDSRow, error) {
	f, err := os.Open(cdsfile)
	if err != nil {
		return nil, fmt.Errorf("failed to open cdsfile: %w", err)
	}
	defer f.Close()

	r := csv.NewReader(f)
	r.Comma = '\t'
	r.FieldsPerRecord = -1
	recs, err := r.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("failed to read cdsfile: %w", err)
	}
	if len(recs) < 2 {
		return nil, fmt.Errorf("cdsfile has no data rows")
	}

	rows := make([]CDSRow, 0, len(recs)-1)
	for i := 1; i < len(recs); i++ {
		rec := recs[i]
		if len(rec) < 10 {
			continue
		}
		row := CDSRow{
			GeneID:   rec[0],
			GeneName: rec[1],
			CDSID:    rec[2],
			Chr:      rec[3],
		}
		if row.ChrCodingStart, err = parseInt(rec[4]); err != nil {
			continue
		}
		if row.ChrCodingEnd, err = parseInt(rec[5]); err != nil {
			continue
		}
		if row.CDSStart, err = parseInt(rec[6]); err != nil {
			continue
		}
		if row.CDSEnd, err = parseInt(rec[7]); err != nil {
			continue
		}
		if row.Length, err = parseInt(rec[8]); err != nil {
			continue
		}
		if row.Strand, err = parseInt(rec[9]); err != nil {
			continue
		}
		rows = append(rows, row)
	}
	return rows, nil
}

func findFullCDS(rows []CDSRow) map[string]struct{} {
	startSet := map[string]struct{}{}
	endSet := map[string]struct{}{}
	for _, r := range rows {
		if r.CDSStart == 1 {
			startSet[r.CDSID] = struct{}{}
		}
		if r.CDSEnd == r.Length {
			endSet[r.CDSID] = struct{}{}
		}
	}
	full := make(map[string]struct{})
	for id := range startSet {
		if _, ok := endSet[id]; ok {
			full[id] = struct{}{}
		}
	}
	return full
}

// Struct for parsing fasta
type Contig struct {
	Name   string
	Seq    string
	Length int
}

// Get headers, seq, length of sequence from a fasta file
func readFasta(path string) ([]Contig, error) {
	const MB = 1024 * 1024

	file, err := os.Open(path)

	if err != nil {
		return nil, err
	}

	defer file.Close()

	contigs := make([]Contig, 0)

	scanner := bufio.NewScanner(file)
	initial := make([]byte, 0, MB)
	scanner.Buffer(initial, 32*MB)

	buffer := strings.Builder{}
	var contig Contig
	haveHeader := false

	flush := func() {
		if !haveHeader {
			return
		}
		contig.Seq = strings.ToUpper(buffer.String())
		contig.Length = len(contig.Seq)
		contigs = append(contigs, contig)
		buffer.Reset()
	}

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, ">") {
			flush()
			header := strings.TrimPrefix(line, ">")
			parts := strings.Fields(header)
			if len(parts) == 0 {
				return contigs, fmt.Errorf("couldn't parse fasta header: %s", line)
			}
			haveHeader = true
			contig = Contig{Name: parts[0]}
			buffer.Reset()
			contig.Name = parts[0]
			continue
		}

		if !haveHeader {
			return contigs, fmt.Errorf("fasta sequence found before header: %s", line)
		}

		buffer.WriteString(line)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	flush()

	if len(contigs) == 0 {
		return contigs, fmt.Errorf("no sequences found in fasta")
	}

	return contigs, nil
}

// readIndex returns the ordered list of chromosome names from the FASTA index
// (.fai). Falls back to parsing the FASTA directly when the index is absent;
// in that case names are sorted alphabetically.
func readIndex(fasta string) ([]string, error) {
	data, err := os.ReadFile(fasta + ".fai")
	if err != nil {
		// No index: fall back to parsing the FASTA file.
		contigs, ef := readFasta(fasta)
		if ef != nil {
			return nil, fmt.Errorf("failed to get chromosomes via both methods: fai: %w; fasta: %w", err, ef)
		}
		chrs := make([]string, 0, len(contigs))
		for _, c := range contigs {
			chrs = append(chrs, c.Name)
		}
		sort.Strings(chrs)
		return chrs, err
	}

	// Parse tab-separated .fai: first field of every non-empty line is the name.
	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	chrs := make([]string, 0, len(lines))
	for _, line := range lines {
		if line == "" {
			continue
		}
		fields := strings.SplitN(line, "\t", 2)
		if len(fields) == 0 || fields[0] == "" {
			return nil, fmt.Errorf("couldn't parse index line: %s", line)
		}
		chrs = append(chrs, fields[0])
	}
	return chrs, nil
}

type Interval struct {
	Start int
	End   int
}

type RefCDSEntry struct {
	GeneName        string
	GeneID          string
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
	CodonImpact     []int
	CodonRates      []int
	CodonExpectedNS []float64
	L               [192][4]int
}

type GeneRange struct {
	Index int
	Chr   string
	Start int
	End   int
	Gene  string
}

type BuildRefResult struct {
	RefCDS  []RefCDSEntry
	GRGenes []GeneRange
}

type codonMeta struct {
	trinucs     []string
	trinucIndex map[string]int
	subsIndex   map[string]int
	impact      [64][64]int
}

type CDSProcessor struct {
	rows      []CDSRow
	genome    map[string]string
	lengths   map[string]int
	validChrs []string
}

func NewCDSProcessor(rows []CDSRow) *CDSProcessor {
	return &CDSProcessor{
		rows: rows,
	}
}

func (p *CDSProcessor) WithGenome(genome map[string]string, lengths map[string]int) *CDSProcessor {
	p.genome = genome
	p.lengths = lengths
	return p
}

func (p *CDSProcessor) WithValidChromosomes(chrs []string) *CDSProcessor {
	p.validChrs = chrs
	return p
}

func (p *CDSProcessor) FilterExcludedChromosomes(excludeChrs []string) *CDSProcessor {
	if len(excludeChrs) == 0 {
		return p
	}
	ex := make(map[string]struct{}, len(excludeChrs))
	for _, c := range excludeChrs {
		ex[c] = struct{}{}
	}
	filtered := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if _, ok := ex[r.Chr]; !ok {
			filtered = append(filtered, r)
		}
	}
	p.rows = filtered
	return p
}

func (p *CDSProcessor) FilterOnlyChromosomes(onlyChrs []string) *CDSProcessor {
	if len(onlyChrs) == 0 {
		return p
	}
	only := make(map[string]struct{}, len(onlyChrs))
	for _, c := range onlyChrs {
		only[c] = struct{}{}
	}
	filtered := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if _, ok := only[r.Chr]; ok {
			filtered = append(filtered, r)
		}
	}
	p.rows = filtered
	return p
}

func (p *CDSProcessor) NormalizeChromosomeNames(genomeChrs []string) *CDSProcessor {
	reftableChrs := uniqueStr(func() []string {
		out := make([]string, 0, len(p.rows))
		for _, r := range p.rows {
			out = append(out, r.Chr)
		}
		return out
	}())

	chrSet := make(map[string]struct{}, len(reftableChrs))
	for _, c := range reftableChrs {
		chrSet[c] = struct{}{}
	}

	hasOverlap := false
	for _, c := range genomeChrs {
		if _, ok := chrSet[c]; ok {
			hasOverlap = true
			break
		}
	}

	if !hasOverlap {
		// Try adding "chr" prefix
		for i := range p.rows {
			p.rows[i].Chr = "chr" + p.rows[i].Chr
		}
		chrSet = make(map[string]struct{}, len(p.rows))
		for _, r := range p.rows {
			chrSet[r.Chr] = struct{}{}
		}
		filteredChrs := make([]string, 0, len(genomeChrs))
		for _, c := range genomeChrs {
			if _, ok := chrSet[c]; ok {
				filteredChrs = append(filteredChrs, c)
			}
		}
		genomeChrs = filteredChrs
		if len(genomeChrs) == 0 {
			// Keep rows but record that normalization failed
			return p
		}
	}

	// Keep only chromosomes in common
	chrSet = make(map[string]struct{}, len(genomeChrs))
	for _, c := range genomeChrs {
		chrSet[c] = struct{}{}
	}
	filtered := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if _, ok := chrSet[r.Chr]; ok {
			filtered = append(filtered, r)
		}
	}
	p.rows = filtered
	return p
}

func (p *CDSProcessor) RemoveIncompleteRecords() *CDSProcessor {
	validChrSet := make(map[string]struct{}, len(p.validChrs))
	for _, c := range p.validChrs {
		validChrSet[c] = struct{}{}
	}

	cleaned := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if r.GeneID == "" || r.GeneName == "" || r.CDSID == "" {
			continue
		}
		if _, ok := validChrSet[r.Chr]; !ok {
			continue
		}
		cleaned = append(cleaned, r)
	}
	p.rows = cleaned
	return p
}

func (p *CDSProcessor) ValidateCoordinates() (*CDSProcessor, error) {
	if len(p.lengths) == 0 {
		return p, fmt.Errorf("genome lengths not set")
	}

	outside := 0
	inside := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		chrLen, ok := p.lengths[r.Chr]
		if !ok || r.ChrCodingStart < 1 || r.ChrCodingEnd > chrLen || r.ChrCodingStart > r.ChrCodingEnd {
			outside++
			continue
		}
		inside = append(inside, r)
	}
	p.rows = inside
	if outside > 0 {
		return p, fmt.Errorf("aborting buildref. %d rows in cdsfile have coordinates that fall outside of the corresponding chromosome length. please ensure that you are using the same assembly for the cdsfile and genomefile", outside)
	}
	return p, nil
}

func (p *CDSProcessor) FindFullCDS() map[string]struct{} {
	return findFullCDS(p.rows)
}

// TrimBoundaryBases removes first 3 bases at chromosome start and last 3 bases at chromosome end
func (p *CDSProcessor) TrimBoundaryBases() *CDSProcessor {
	for i := range p.rows {
		if p.rows[i].ChrCodingStart == 1 {
			p.rows[i].ChrCodingStart += 3
			p.rows[i].CDSStart += 3
		}
		if len(p.lengths) > 0 && p.lengths[p.rows[i].Chr] == p.rows[i].ChrCodingEnd {
			p.rows[i].ChrCodingEnd -= 3
			p.rows[i].CDSEnd -= 3
		}
	}
	return p
}

// DeduplicateByCDSKey removes duplicate CDS records, keeping only the key fields
// TODO: Maybe add this to preprocessing, i.e CDSRow should be mapset?
func uniqueCDSKeyRows(rows []CDSRow) []CDSRow {
	seen := map[string]struct{}{}
	out := make([]CDSRow, 0, len(rows))
	for _, r := range rows {
		k := strings.Join([]string{r.GeneID, r.GeneName, r.CDSID, strconv.Itoa(r.Length)}, "\t")
		if _, ok := seen[k]; ok {
			continue
		}
		seen[k] = struct{}{}
		out = append(out, CDSRow{GeneID: r.GeneID, GeneName: r.GeneName, CDSID: r.CDSID, Length: r.Length})
	}
	return out
}
func (p *CDSProcessor) DeduplicateByCDSKey() *CDSProcessor {
	p.rows = uniqueCDSKeyRows(p.rows)
	return p
}

// SortByGeneAndLength sorts records by gene name and length (descending)
func (p *CDSProcessor) SortByGeneAndLength() *CDSProcessor {
	sort.SliceStable(p.rows, func(i, j int) bool {
		if p.rows[i].GeneName == p.rows[j].GeneName {
			return p.rows[i].Length > p.rows[j].Length
		}
		return p.rows[i].GeneName < p.rows[j].GeneName
	})
	return p
}

// FilterByFullCDS keeps only CDS entries that are in the full CDS set
func (p *CDSProcessor) FilterByFullCDS(fullCDS map[string]struct{}) *CDSProcessor {
	filtered := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if _, ok := fullCDS[r.CDSID]; ok {
			filtered = append(filtered, r)
		}
	}
	p.rows = filtered
	return p
}

// FilterValidCDSLength keeps only CDS with length divisible by 3
func (p *CDSProcessor) FilterValidCDSLength() *CDSProcessor {
	filtered := make([]CDSRow, 0, len(p.rows))
	for _, r := range p.rows {
		if r.Length%3 == 0 {
			filtered = append(filtered, r)
		}
	}
	p.rows = filtered
	return p
}

// SortByCoordinates sorts records by chromosome, coding start, and CDSID
func (p *CDSProcessor) SortByCoordinates() *CDSProcessor {
	sort.Slice(p.rows, func(i, j int) bool {
		if p.rows[i].Chr == p.rows[j].Chr {
			if p.rows[i].ChrCodingStart == p.rows[j].ChrCodingStart {
				return p.rows[i].CDSID < p.rows[j].CDSID
			}
			return p.rows[i].ChrCodingStart < p.rows[j].ChrCodingStart
		}
		return p.rows[i].Chr < p.rows[j].Chr
	})
	return p
}

// Rows returns the processed CDS rows
func (p *CDSProcessor) Rows() []CDSRow {
	return p.rows
}

func buildRef(cdsFile, genomeFile string, excludeChrs, onlyChrs []string, useIDs bool) (BuildRefResult, error) {

	reftable, err := readCDS(cdsFile)
	if err != nil {
		return BuildRefResult{}, err
	}

	for i := range reftable {
		if useIDs {
			reftable[i].GeneName = reftable[i].GeneID + ":" + reftable[i].GeneName
		}
	}

	validChrs, err := readIndex(genomeFile)
	if err != nil {
		return BuildRefResult{}, err
	}

	// Read genome and chromosome lengths
	contigs, err := readFasta(genomeFile)
	if err != nil {
		return BuildRefResult{}, err
	}
	genome := make(map[string]string, len(contigs))
	lengths := make(map[string]int, len(contigs))
	for _, c := range contigs {
		genome[c.Name] = c.Seq
		lengths[c.Name] = c.Length
	}

	processor := NewCDSProcessor(reftable).
		FilterExcludedChromosomes(excludeChrs).
		FilterOnlyChromosomes(onlyChrs).
		NormalizeChromosomeNames(validChrs).
		WithValidChromosomes(validChrs).
		WithGenome(genome, lengths).
		RemoveIncompleteRecords()

	// Validate coordinates
	var errVal error
	processor, errVal = processor.ValidateCoordinates()
	if errVal != nil {
		return BuildRefResult{}, errVal
	}

	// Find full CDS and apply boundary trimming
	fullCDS := processor.FindFullCDS()
	processor.TrimBoundaryBases()

	// Create separate processor for cdsTable (unique CDS records)
	cdsTableProcessor := NewCDSProcessor(processor.Rows()).
		DeduplicateByCDSKey().
		SortByGeneAndLength().
		FilterValidCDSLength().
		FilterByFullCDS(fullCDS)

	// Prepare reftable (all exon records)
	reftable = NewCDSProcessor(processor.Rows()).
		FilterByFullCDS(fullCDS).
		SortByCoordinates().
		Rows()

	cdsSplit := splitByCDSID(reftable)
	geneSplit := splitByGeneName(cdsTableProcessor.Rows())
	geneOrder := sortedKeys(geneSplit)

	ref := make([]RefCDSEntry, 0, len(geneOrder))
	meta := initCodonMeta()

	for _, geneName := range geneOrder {
		geneCDSs := geneSplit[geneName]
		valid := false
		for _, candidate := range geneCDSs {
			cds := cdsSplit[candidate.CDSID]
			if len(cds) == 0 {
				continue
			}
			strand := cds[0].Strand
			chr := cds[0].Chr

			intervals := make([]Interval, 0, len(cds))
			for _, exon := range cds {
				intervals = append(intervals, Interval{Start: exon.ChrCodingStart, End: exon.ChrCodingEnd})
			}

			cdsSeq, err := getCDSSeq(genome, chr, intervals, strand)
			if err != nil {
				continue
			}
			if containsBase(cdsSeq, 'N') {
				continue
			}

			pep, err := translateDNA(cdsSeq)
			if err != nil || len(pep) == 0 {
				continue
			}
			if strings.ContainsRune(pep[:len(pep)-1], '*') {
				continue
			}

			splPos := getSpliceSites(cds)
			splSeq := ""
			splSeq1Up := ""
			splSeq1Down := ""

			cdsSeq1Up, err := getCDSContext(genome, chr, intervals, strand, true)
			if err != nil {
				continue
			}
			cdsSeq1Down, err := getCDSContext(genome, chr, intervals, strand, false)
			if err != nil {
				continue
			}

			if len(splPos) > 0 {
				splSeq, err = getSpliceSeq(genome, chr, splPos, strand)
				if err != nil {
					continue
				}
				if strand == 1 {
					splSeq1Up, err = getSpliceSeq(genome, chr, addToAll(splPos, -1), strand)
					if err != nil {
						continue
					}
					splSeq1Down, err = getSpliceSeq(genome, chr, addToAll(splPos, 1), strand)
					if err != nil {
						continue
					}
				} else {
					splSeq1Up, err = getSpliceSeq(genome, chr, addToAll(splPos, 1), strand)
					if err != nil {
						continue
					}
					splSeq1Down, err = getSpliceSeq(genome, chr, addToAll(splPos, -1), strand)
					if err != nil {
						continue
					}
				}
			}

			entry := RefCDSEntry{
				GeneName:        candidate.GeneName,
				GeneID:          candidate.GeneID,
				ProteinID:       candidate.CDSID,
				CDSLength:       candidate.Length,
				Chr:             chr,
				Strand:          strand,
				IntervalsCDS:    intervals,
				IntervalsSplice: splPos,
				SeqCDS:          cdsSeq,
				SeqCDS1Up:       cdsSeq1Up,
				SeqCDS1Down:     cdsSeq1Down,
				SeqSplice:       splSeq,
				SeqSplice1Up:    splSeq1Up,
				SeqSplice1Down:  splSeq1Down,
			}
			entry.L, err = buildLMatrix(entry, meta)
			if err != nil {
				continue
			}

			ref = append(ref, entry)
			valid = true
			break
		}
		if !valid {
			continue
		}
	}

	geneRanges := make([]GeneRange, 0)
	for i := range ref {
		idx := i + 1
		for _, iv := range ref[i].IntervalsCDS {
			geneRanges = append(geneRanges, GeneRange{Index: idx, Chr: ref[i].Chr, Start: iv.Start, End: iv.End, Gene: ref[i].GeneName})
		}
		for _, p := range ref[i].IntervalsSplice {
			geneRanges = append(geneRanges, GeneRange{Index: idx, Chr: ref[i].Chr, Start: p, End: p, Gene: ref[i].GeneName})
		}
	}

	return BuildRefResult{RefCDS: ref, GRGenes: geneRanges}, nil
}

func parseInt(v string) (int, error) {
	v = strings.TrimSpace(v)
	if v == "" || v == "-" {
		return 0, fmt.Errorf("missing numeric")
	}
	n, err := strconv.Atoi(v)
	if err != nil {
		return 0, err
	}
	return n, nil
}

func filterCDSTable(rows []CDSRow, full map[string]struct{}) []CDSRow {
	out := make([]CDSRow, 0, len(rows))
	for _, r := range rows {
		if r.Length%3 != 0 {
			continue
		}
		if _, ok := full[r.CDSID]; !ok {
			continue
		}
		out = append(out, r)
	}
	return out
}

func filterRowsByCDS(rows []CDSRow, full map[string]struct{}) []CDSRow {
	out := make([]CDSRow, 0, len(rows))
	for _, r := range rows {
		if _, ok := full[r.CDSID]; ok {
			out = append(out, r)
		}
	}
	return out
}

func splitByCDSID(rows []CDSRow) map[string][]CDSRow {
	m := map[string][]CDSRow{}
	for _, r := range rows {
		m[r.CDSID] = append(m[r.CDSID], r)
	}
	return m
}

func splitByGeneName(rows []CDSRow) map[string][]CDSRow {
	m := map[string][]CDSRow{}
	for _, r := range rows {
		m[r.GeneName] = append(m[r.GeneName], r)
	}
	return m
}

func sortedKeys[T any](m map[string]T) []string {
	out := make([]string, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	sort.Strings(out)
	return out
}

func getCDSSeq(genome map[string]string, chr string, intervals []Interval, strand int) (string, error) {
	seq, ok := genome[chr]
	if !ok {
		return "", fmt.Errorf("chromosome not found: %s", chr)
	}
	var b strings.Builder
	for _, iv := range intervals {
		if iv.Start < 1 || iv.End > len(seq) || iv.Start > iv.End {
			return "", fmt.Errorf("interval out of range: %s:%d-%d", chr, iv.Start, iv.End)
		}
		b.WriteString(seq[iv.Start-1 : iv.End])
	}
	out := b.String()
	if strand == -1 {
		out = reverseComplement(out)
	}
	return out, nil
}

func getCDSContext(genome map[string]string, chr string, intervals []Interval, strand int, upstream bool) (string, error) {
	offset := -1
	if !upstream {
		offset = 1
	}
	if strand == -1 {
		offset = -offset
	}
	shifted := make([]Interval, len(intervals))
	for i, iv := range intervals {
		shifted[i] = Interval{Start: iv.Start + offset, End: iv.End + offset}
	}
	return getCDSSeq(genome, chr, shifted, strand)
}

func getSpliceSites(cds []CDSRow) []int {
	if len(cds) <= 1 {
		return nil
	}
	strand := cds[0].Strand
	positions := make([]int, 0, len(cds)*5)
	if strand == 1 {
		for i := 0; i < len(cds)-1; i++ {
			spl5 := cds[i].ChrCodingEnd
			positions = append(positions, spl5+1, spl5+2, spl5+5)
		}
		for i := 1; i < len(cds); i++ {
			spl3 := cds[i].ChrCodingStart
			positions = append(positions, spl3-1, spl3-2)
		}
	} else {
		for i := 1; i < len(cds); i++ {
			spl5 := cds[i].ChrCodingStart
			positions = append(positions, spl5-1, spl5-2, spl5-5)
		}
		for i := 0; i < len(cds)-1; i++ {
			spl3 := cds[i].ChrCodingEnd
			positions = append(positions, spl3+1, spl3+2)
		}
	}
	sort.Ints(positions)
	return uniqueInt(positions)
}

func getSpliceSeq(genome map[string]string, chr string, positions []int, strand int) (string, error) {
	seq, ok := genome[chr]
	if !ok {
		return "", fmt.Errorf("chromosome not found: %s", chr)
	}
	buf := make([]byte, 0, len(positions))
	for _, p := range positions {
		if p < 1 || p > len(seq) {
			return "", fmt.Errorf("splice position out of range: %s:%d", chr, p)
		}
		buf = append(buf, seq[p-1])
	}
	out := string(buf)
	if strand == -1 {
		out = complement(out)
	}
	return out, nil
}

func addToAll(nums []int, delta int) []int {
	out := make([]int, len(nums))
	for i, n := range nums {
		out[i] = n + delta
	}
	return out
}

func complement(seq string) string {
	b := []byte(strings.ToUpper(seq))
	for i := range b {
		b[i] = compBase(b[i])
	}
	return string(b)
}

func reverseComplement(seq string) string {
	b := []byte(strings.ToUpper(seq))
	for i, j := 0, len(b)-1; i <= j; i, j = i+1, j-1 {
		bi := compBase(b[j])
		bj := compBase(b[i])
		b[i], b[j] = bi, bj
	}
	return string(b)
}

func compBase(b byte) byte {
	switch b {
	case 'A':
		return 'T'
	case 'C':
		return 'G'
	case 'G':
		return 'C'
	case 'T':
		return 'A'
	default:
		return 'N'
	}
}

var standardCode = map[string]rune{
	"TTT": 'F', "TTC": 'F', "TTA": 'L', "TTG": 'L', "CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	"ATT": 'I', "ATC": 'I', "ATA": 'I', "ATG": 'M', "GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', "CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T', "GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"TAT": 'Y', "TAC": 'Y', "TAA": '*', "TAG": '*', "CAT": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q',
	"AAT": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K', "GAT": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E',
	"TGT": 'C', "TGC": 'C', "TGA": '*', "TGG": 'W', "CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
	"AGT": 'S', "AGC": 'S', "AGA": 'R', "AGG": 'R', "GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
}

func translateDNA(seq string) (string, error) {
	if len(seq)%3 != 0 {
		return "", fmt.Errorf("sequence length is not a multiple of 3")
	}
	seq = strings.ToUpper(seq)
	out := make([]rune, 0, len(seq)/3)
	for i := 0; i < len(seq); i += 3 {
		codon := seq[i : i+3]
		aa, ok := standardCode[codon]
		if !ok {
			return "", fmt.Errorf("unsupported codon: %s", codon)
		}
		out = append(out, aa)
	}
	return string(out), nil
}

func containsBase(seq string, base byte) bool {
	for i := 0; i < len(seq); i++ {
		if seq[i] == base {
			return true
		}
	}
	return false
}

func initCodonMeta() codonMeta {
	nt := []byte{'A', 'C', 'G', 'T'}
	trinucs := make([]string, 0, 64)
	for _, a := range nt {
		for _, b := range nt {
			for _, c := range nt {
				trinucs = append(trinucs, string([]byte{a, b, c}))
			}
		}
	}
	trinucIndex := make(map[string]int, 64)
	for i, t := range trinucs {
		trinucIndex[t] = i
	}

	subsIndex := make(map[string]int, 192)
	idx := 0
	for _, old := range trinucs {
		for _, newBase := range []byte{'A', 'C', 'G', 'T'} {
			if old[1] == newBase {
				continue
			}
			newTri := string([]byte{old[0], newBase, old[2]})
			subsIndex[old+">"+newTri] = idx
			idx++
		}
	}

	var impact [64][64]int
	for i, from := range trinucs {
		fromAA, _ := translateDNA(from)
		for j, to := range trinucs {
			toAA, _ := translateDNA(to)
			switch {
			case toAA == fromAA:
				impact[i][j] = 1
			case toAA == "*":
				impact[i][j] = 3
			case toAA != "*" && fromAA != "*" && toAA != fromAA:
				impact[i][j] = 2
			default:
				impact[i][j] = 0
			}
		}
	}

	return codonMeta{trinucs: trinucs, trinucIndex: trinucIndex, subsIndex: subsIndex, impact: impact}
}

func buildLMatrix(entry RefCDSEntry, meta codonMeta) ([192][4]int, error) {
	var L [192][4]int
	nt := []byte{'A', 'C', 'G', 'T'}

	cds := []byte(entry.SeqCDS)
	up := []byte(entry.SeqCDS1Up)
	down := []byte(entry.SeqCDS1Down)
	if len(cds) != len(up) || len(cds) != len(down) {
		return L, fmt.Errorf("cds context lengths mismatch")
	}
	if len(cds)%3 != 0 {
		return L, fmt.Errorf("cds length is not multiple of 3")
	}

	for pos := 0; pos < len(cds); pos++ {
		oldTri := string([]byte{up[pos], cds[pos], down[pos]})
		codonStart := (pos / 3) * 3
		oldCodon := string(cds[codonStart : codonStart+3])
		for _, nb := range nt {
			if nb == cds[pos] {
				continue
			}
			newTri := string([]byte{up[pos], nb, down[pos]})
			newCodonBytes := []byte(oldCodon)
			newCodonBytes[pos%3] = nb
			newCodon := string(newCodonBytes)

			fromIdx, ok1 := meta.trinucIndex[oldCodon]
			toIdx, ok2 := meta.trinucIndex[newCodon]
			matIdx, ok3 := meta.subsIndex[oldTri+">"+newTri]
			if !ok1 || !ok2 || !ok3 {
				continue
			}

			impact := meta.impact[fromIdx][toIdx]
			switch impact {
			case 1:
				L[matIdx][0]++
			case 2:
				L[matIdx][1]++
			case 3:
				L[matIdx][2]++
			}
		}
	}

	if len(entry.IntervalsSplice) > 0 {
		spl := []byte(entry.SeqSplice)
		sup := []byte(entry.SeqSplice1Up)
		sdown := []byte(entry.SeqSplice1Down)
		if len(spl) != len(sup) || len(spl) != len(sdown) {
			return L, fmt.Errorf("splice context lengths mismatch")
		}
		for i := 0; i < len(spl); i++ {
			oldTri := string([]byte{sup[i], spl[i], sdown[i]})
			for _, nb := range nt {
				if nb == spl[i] {
					continue
				}
				newTri := string([]byte{sup[i], nb, sdown[i]})
				idx, ok := meta.subsIndex[oldTri+">"+newTri]
				if !ok {
					continue
				}
				L[idx][3]++
			}
		}
	}

	return L, nil
}

func uniqueStr(items []string) []string {
	seen := map[string]struct{}{}
	out := make([]string, 0, len(items))
	for _, item := range items {
		if _, ok := seen[item]; ok {
			continue
		}
		seen[item] = struct{}{}
		out = append(out, item)
	}
	return out
}

func uniqueInt(items []int) []int {
	if len(items) == 0 {
		return nil
	}
	out := make([]int, 0, len(items))
	last := items[0] - 1
	for _, item := range items {
		if item != last {
			out = append(out, item)
			last = item
		}
	}
	return out
}

type jsonSummary struct {
	Genes []string `json:"genes"`
}

func (b BuildRefResult) MarshalJSON() ([]byte, error) {
	genes := make([]string, 0, len(b.RefCDS))
	for _, g := range b.RefCDS {
		genes = append(genes, g.GeneName)
	}
	slices.Sort(genes)
	return json.Marshal(jsonSummary{Genes: genes})
}
