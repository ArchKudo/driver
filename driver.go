package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
)

type Mutation struct {
	sampleID string
	chr      string
	pos      string
	ref      string
	alt      string
}

func readSampleSheet(filename string) []Mutation {
	file, err := os.Open(filename)

	if err != nil {
		panic(fmt.Errorf("Couldn't load samplesheet: %w", err))
	}

	defer file.Close()

	reader := csv.NewReader(file)

	// skip header
	_, err = reader.Read()

	if err != nil {
		panic(fmt.Errorf("Error at reading header: %w", err))
	}

	// read records
	records, err := reader.ReadAll()

	if err != nil {
		panic(fmt.Errorf("Error at reading records: %w", err))
	}

	if len(records) == 0 {
		panic(fmt.Errorf("No records found!"))
	}

	mutations := make([]Mutation, 0, len(records))

	for i, r := range records {

		if len(r) != 5 {
			panic(fmt.Errorf("Found %d columns expected 5 at pos %d", len(r), i))
		}

		mut := Mutation{
			sampleID: r[0],
			chr:      r[1],
			pos:      r[2],
			ref:      r[3],
			alt:      r[4],
		}

		mutations = append(mutations, mut)

	}

	return mutations

}

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

func readCDSTable(cdsfile string) ([]CDSRow, error) {
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

func buildRef(cdsFile, genomeFile string, excludeChrs, onlyChrs []string, useIDs bool) (BuildRefResult, error) {

	reftable, err := readCDSTable(cdsFile)
	if err != nil {
		return BuildRefResult{}, err
	}

	for i := range reftable {
		if useIDs {
			reftable[i].GeneName = reftable[i].GeneID + ":" + reftable[i].GeneName
		}
	}

	validChrs, err := fastaIndexChromosomes(genomeFile)
	if err != nil {
		return BuildRefResult{}, err
	}

	if len(excludeChrs) > 0 {
		ex := make(map[string]struct{}, len(excludeChrs))
		for _, c := range excludeChrs {
			ex[c] = struct{}{}
		}
		filtered := make([]string, 0, len(validChrs))
		for _, c := range validChrs {
			if _, ok := ex[c]; !ok {
				filtered = append(filtered, c)
			}
		}
		validChrs = filtered
	}

	if len(onlyChrs) > 0 {
		only := make(map[string]struct{}, len(onlyChrs))
		for _, c := range onlyChrs {
			only[c] = struct{}{}
		}
		filtered := make([]string, 0, len(validChrs))
		for _, c := range validChrs {
			if _, ok := only[c]; ok {
				filtered = append(filtered, c)
			}
		}
		validChrs = filtered
	}

	reftableChrs := uniqueStr(func() []string {
		out := make([]string, 0, len(reftable))
		for _, r := range reftable {
			out = append(out, r.Chr)
		}
		return out
	}())

	hasOverlap := false
	chrSet := make(map[string]struct{}, len(reftableChrs))
	for _, c := range reftableChrs {
		chrSet[c] = struct{}{}
	}
	for _, c := range validChrs {
		if _, ok := chrSet[c]; ok {
			hasOverlap = true
			break
		}
	}

	if !hasOverlap {
		for i := range reftable {
			reftable[i].Chr = "chr" + reftable[i].Chr
		}
		chrSet = make(map[string]struct{}, len(reftable))
		for _, r := range reftable {
			chrSet[r.Chr] = struct{}{}
		}
		filtered := make([]string, 0, len(validChrs))
		for _, c := range validChrs {
			if _, ok := chrSet[c]; ok {
				filtered = append(filtered, c)
			}
		}
		validChrs = filtered
		if len(validChrs) == 0 {
			return BuildRefResult{}, fmt.Errorf("no chromosome names in common between the genome file and the CDS table")
		}
	} else {
		filtered := make([]string, 0, len(validChrs))
		for _, c := range validChrs {
			if _, ok := chrSet[c]; ok {
				filtered = append(filtered, c)
			}
		}
		validChrs = filtered
	}

	validChrSet := make(map[string]struct{}, len(validChrs))
	for _, c := range validChrs {
		validChrSet[c] = struct{}{}
	}

	cleaned := make([]CDSRow, 0, len(reftable))
	for _, r := range reftable {
		if r.GeneID == "" || r.GeneName == "" || r.CDSID == "" {
			continue
		}
		if _, ok := validChrSet[r.Chr]; !ok {
			continue
		}
		cleaned = append(cleaned, r)
	}
	reftable = cleaned

	genome, lengths, err := loadFasta(genomeFile)
	if err != nil {
		return BuildRefResult{}, err
	}

	outside := 0
	inside := make([]CDSRow, 0, len(reftable))
	for _, r := range reftable {
		chrLen, ok := lengths[r.Chr]
		if !ok || r.ChrCodingStart < 1 || r.ChrCodingEnd > chrLen || r.ChrCodingStart > r.ChrCodingEnd {
			outside++
			continue
		}
		inside = append(inside, r)
	}
	if outside > 0 {
		return BuildRefResult{}, fmt.Errorf("aborting buildref. %d rows in cdsfile have coordinates that fall outside of the corresponding chromosome length. please ensure that you are using the same assembly for the cdsfile and genomefile", outside)
	}
	reftable = inside

	fullCDS := findFullCDS(reftable)

	for i := range reftable {
		if reftable[i].ChrCodingStart == 1 {
			reftable[i].ChrCodingStart += 3
			reftable[i].CDSStart += 3
		}
		if lengths[reftable[i].Chr] == reftable[i].ChrCodingEnd {
			reftable[i].ChrCodingEnd -= 3
			reftable[i].CDSEnd -= 3
		}
	}

	cdsTable := uniqueCDSKeyRows(reftable)
	sort.SliceStable(cdsTable, func(i, j int) bool {
		if cdsTable[i].GeneName == cdsTable[j].GeneName {
			return cdsTable[i].Length > cdsTable[j].Length
		}
		return cdsTable[i].GeneName < cdsTable[j].GeneName
	})
	cdsTable = filterCDSTable(cdsTable, fullCDS)
	reftable = filterRowsByCDS(reftable, fullCDS)

	sort.Slice(reftable, func(i, j int) bool {
		if reftable[i].Chr == reftable[j].Chr {
			if reftable[i].ChrCodingStart == reftable[j].ChrCodingStart {
				return reftable[i].CDSID < reftable[j].CDSID
			}
			return reftable[i].ChrCodingStart < reftable[j].ChrCodingStart
		}
		return reftable[i].Chr < reftable[j].Chr
	})

	cdsSplit := splitByCDSID(reftable)
	geneSplit := splitByGeneName(cdsTable)
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

func fastaIndexChromosomes(genomeFile string) ([]string, error) {
	faiPath := genomeFile + ".fai"
	if _, err := os.Stat(faiPath); err == nil {
		f, err := os.Open(faiPath)
		if err != nil {
			return nil, err
		}
		defer f.Close()
		s := bufio.NewScanner(f)
		chrs := make([]string, 0)
		for s.Scan() {
			line := strings.TrimSpace(s.Text())
			if line == "" {
				continue
			}
			fields := strings.Split(line, "\t")
			if len(fields) == 0 || fields[0] == "" {
				continue
			}
			chrs = append(chrs, fields[0])
		}
		if err := s.Err(); err != nil {
			return nil, err
		}
		return chrs, nil
	}

	_, lengths, err := loadFasta(genomeFile)
	if err != nil {
		return nil, err
	}
	chrs := make([]string, 0, len(lengths))
	for c := range lengths {
		chrs = append(chrs, c)
	}
	sort.Strings(chrs)
	return chrs, nil
}

func loadFasta(path string) (map[string]string, map[string]int, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()

	seqs := map[string]string{}
	lengths := map[string]int{}

	s := bufio.NewScanner(f)
	buf := make([]byte, 0, 1024*1024)
	s.Buffer(buf, 32*1024*1024)

	cur := ""
	b := strings.Builder{}
	flush := func() {
		if cur == "" {
			return
		}
		seq := strings.ToUpper(b.String())
		seqs[cur] = seq
		lengths[cur] = len(seq)
		b.Reset()
	}

	for s.Scan() {
		line := strings.TrimSpace(s.Text())
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			flush()
			header := strings.TrimPrefix(line, ">")
			parts := strings.Fields(header)
			if len(parts) == 0 {
				return nil, nil, fmt.Errorf("invalid fasta header")
			}
			cur = parts[0]
			continue
		}
		if cur == "" {
			return nil, nil, fmt.Errorf("fasta sequence before header")
		}
		b.WriteString(line)
	}
	if err := s.Err(); err != nil {
		return nil, nil, err
	}
	flush()

	if len(seqs) == 0 {
		return nil, nil, fmt.Errorf("no sequences found in fasta")
	}
	return seqs, lengths, nil
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

func main() {

	file := filepath.Join("data", "simple_breast.csv")

	sample_sheet := readSampleSheet(file)

	for i := 0; i < 5; i++ {
		fmt.Println(sample_sheet[i])
	}
}
