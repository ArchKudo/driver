package main

import (
	"encoding/csv"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// Struct for columns in sample sheet
type Mutation struct {
	sampleID string
	chr      string
	pos      string
	ref      string
	alt      string
}

// Parse the main input sample sheet
// Psuedo vcf file
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

// Struct for CDS table

type DndsOutput struct {
	RefCDS     []RefCDSEntry         `json:"refcds,omitempty"`
	Mutations  []AnnotatedMutation   `json:"annotmuts,omitempty"`
	GeneMuts   []GeneMutationSummary `json:"genemuts,omitempty"`
	GlobalDNDS map[string]float64    `json:"globaldnds,omitempty"`
	N          [192][4]int           `json:"N,omitempty"`
	L          [192][4]int           `json:"L,omitempty"`
}

type AnnotatedMutation struct {
	SampleID string `json:"sampleID"`
	Chr      string `json:"chr"`
	Pos      int    `json:"pos"`
	Ref      string `json:"ref"`
	Mut      string `json:"mut"`
	Gene     string `json:"gene"`
	Impact   string `json:"impact"`
}

type GeneMutationSummary struct {
	GeneName string `json:"gene_name"`
	NSyn     int    `json:"n_syn"`
	NMis     int    `json:"n_mis"`
	NNon     int    `json:"n_non"`
	NSpl     int    `json:"n_spl"`
	Total    int    `json:"total"`
}

type RecurrentSite struct {
	Chr     string  `json:"chr"`
	Pos     int     `json:"pos"`
	Ref     string  `json:"ref"`
	Mut     string  `json:"mut"`
	Gene    string  `json:"gene"`
	Count   int     `json:"count"`
	Mu      float64 `json:"mu"`
	DNDS    float64 `json:"dnds"`
	PValue  float64 `json:"pvalue"`
	QValue  float64 `json:"qvalue"`
	Context string  `json:"context"`
}

type SiteDNDSOutput struct {
	Sites []RecurrentSite `json:"recursites"`
}

type RecurrentCodon struct {
	Gene   string  `json:"gene"`
	Codon  int     `json:"codon"`
	Count  int     `json:"count"`
	Mu     float64 `json:"mu"`
	DNDS   float64 `json:"dnds"`
	PValue float64 `json:"pvalue"`
	QValue float64 `json:"qvalue"`
}

type CodonDNDSOutput struct {
	Codons []RecurrentCodon `json:"recurcodons"`
}

type GeneSetDNDS struct {
	GeneSet map[string]float64 `json:"globaldnds_geneset"`
	Rest    map[string]float64 `json:"globaldnds_rest"`
}

type GeneCI struct {
	Gene   string  `json:"gene"`
	Value  float64 `json:"value"`
	Low95  float64 `json:"low95"`
	High95 float64 `json:"high95"`
}

type GeneCIOutput struct {
	Rows []GeneCI `json:"geneci"`
}

type WithinGeneResult struct {
	Gene     string         `json:"gene"`
	Regions  []RegionResult `json:"regions"`
	Mutation int            `json:"mutation_count"`
}

type RegionResult struct {
	Label string  `json:"label"`
	Start int     `json:"start"`
	End   int     `json:"end"`
	Obs   int     `json:"obs"`
	Exp   float64 `json:"exp"`
}

type FitLNPResult struct {
	MLE       float64 `json:"mle"`
	CI95Low   float64 `json:"ci95_low"`
	CI95High  float64 `json:"ci95_high"`
	ThetaMode string  `json:"theta_mode"`
}

type serializableInterval struct {
	Start int `json:"start"`
	End   int `json:"end"`
}

type serializableRefCDSEntry struct {
	GeneName        string                 `json:"gene_name"`
	GeneID          string                 `json:"gene_id"`
	ProteinID       string                 `json:"protein_id"`
	CDSLength       int                    `json:"cds_length"`
	Chr             string                 `json:"chr"`
	Strand          int                    `json:"strand"`
	IntervalsCDS    []serializableInterval `json:"intervals_cds"`
	IntervalsSplice []int                  `json:"intervals_splice"`
	SeqCDS          string                 `json:"seq_cds"`
	SeqCDS1Up       string                 `json:"seq_cds1up"`
	SeqCDS1Down     string                 `json:"seq_cds1down"`
	SeqSplice       string                 `json:"seq_splice,omitempty"`
	SeqSplice1Up    string                 `json:"seq_splice1up,omitempty"`
	SeqSplice1Down  string                 `json:"seq_splice1down,omitempty"`
	CodonImpact     []int                  `json:"codon_impact,omitempty"`
	CodonRates      []int                  `json:"codon_rates,omitempty"`
	CodonExpectedNS []float64              `json:"codon_expected_ns,omitempty"`
	L               [192][4]int            `json:"L"`
}

type serializableBuildRef struct {
	RefCDS  []serializableRefCDSEntry `json:"refcds"`
	GRGenes []GeneRange               `json:"gr_genes"`
}

func dndscvEquivalent(mutations []Mutation, ref BuildRefResult, geneList []string, outmats bool) DndsOutput {
	allowed := make(map[string]struct{}, len(geneList))
	if len(geneList) > 0 {
		for _, g := range geneList {
			allowed[g] = struct{}{}
		}
	}

	geneByChr := map[string][]RefCDSEntry{}
	for _, g := range ref.RefCDS {
		if len(allowed) > 0 {
			if _, ok := allowed[g.GeneName]; !ok {
				continue
			}
		}
		geneByChr[g.Chr] = append(geneByChr[g.Chr], g)
	}

	annot := make([]AnnotatedMutation, 0, len(mutations))
	byGene := map[string]*GeneMutationSummary{}
	var N [192][4]int
	var L [192][4]int

	for _, gene := range ref.RefCDS {
		for i := 0; i < 192; i++ {
			for j := 0; j < 4; j++ {
				L[i][j] += gene.L[i][j]
			}
		}
	}

	for _, mut := range mutations {
		pos, err := strconv.Atoi(mut.pos)
		if err != nil {
			continue
		}
		cands := geneByChr[mut.chr]
		if len(cands) == 0 {
			continue
		}
		for _, gene := range cands {
			if !positionInIntervals(pos, gene.IntervalsCDS) {
				continue
			}
			impact := inferImpact(mut.ref, mut.alt)
			annot = append(annot, AnnotatedMutation{
				SampleID: mut.sampleID,
				Chr:      mut.chr,
				Pos:      pos,
				Ref:      mut.ref,
				Mut:      mut.alt,
				Gene:     gene.GeneName,
				Impact:   impact,
			})
			if _, ok := byGene[gene.GeneName]; !ok {
				byGene[gene.GeneName] = &GeneMutationSummary{GeneName: gene.GeneName}
			}
			row := byGene[gene.GeneName]
			switch impact {
			case "Synonymous":
				row.NSyn++
				N[0][0]++
			case "Nonsense":
				row.NNon++
				N[0][2]++
			default:
				row.NMis++
				N[0][1]++
			}
			row.Total++
			break
		}
	}

	geneRows := make([]GeneMutationSummary, 0, len(byGene))
	for _, v := range byGene {
		geneRows = append(geneRows, *v)
	}
	sort.Slice(geneRows, func(i, j int) bool { return geneRows[i].GeneName < geneRows[j].GeneName })

	nSyn, nMis, nNon, lSyn, lMis, lNon := 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
	for _, g := range geneRows {
		nSyn += float64(g.NSyn)
		nMis += float64(g.NMis)
		nNon += float64(g.NNon)
	}
	for i := 0; i < 192; i++ {
		lSyn += float64(L[i][0])
		lMis += float64(L[i][1])
		lNon += float64(L[i][2])
	}
	global := map[string]float64{
		"wmis": safeRatio(nMis, lMis, nSyn, lSyn),
		"wnon": safeRatio(nNon, lNon, nSyn, lSyn),
		"wspl": 1,
	}

	out := DndsOutput{RefCDS: ref.RefCDS, Mutations: annot, GeneMuts: geneRows, GlobalDNDS: global}
	if outmats {
		out.N = N
		out.L = L
	}
	return out
}

func buildcodonEquivalent(ref []RefCDSEntry) []RefCDSEntry {
	meta := initCodonMeta()
	out := make([]RefCDSEntry, len(ref))
	for i, entry := range ref {
		newEntry := entry
		cds := []byte(entry.SeqCDS)
		up := []byte(entry.SeqCDS1Up)
		down := []byte(entry.SeqCDS1Down)
		if len(cds) == 0 || len(cds)%3 != 0 || len(cds) != len(up) || len(cds) != len(down) {
			out[i] = newEntry
			continue
		}
		impact := make([]int, 0, len(cds)*3)
		rates := make([]int, 0, len(cds)*3)
		expectedNS := make([]float64, len(cds)/3)
		for pos := 0; pos < len(cds); pos++ {
			oldTri := string([]byte{up[pos], cds[pos], down[pos]})
			codonStart := (pos / 3) * 3
			oldCodon := string(cds[codonStart : codonStart+3])
			for _, nb := range []byte{'A', 'C', 'G', 'T'} {
				if nb == cds[pos] {
					continue
				}
				newTri := string([]byte{up[pos], nb, down[pos]})
				newCodon := []byte(oldCodon)
				newCodon[pos%3] = nb
				fromIdx, ok1 := meta.trinucIndex[oldCodon]
				toIdx, ok2 := meta.trinucIndex[string(newCodon)]
				rateIdx, ok3 := meta.subsIndex[oldTri+">"+newTri]
				if !ok1 || !ok2 || !ok3 {
					continue
				}
				imp := meta.impact[fromIdx][toIdx]
				impact = append(impact, imp)
				rates = append(rates, rateIdx+1)
				if imp == 2 || imp == 3 {
					expectedNS[codonStart/3]++
				}
			}
		}
		newEntry.CodonImpact = impact
		newEntry.CodonRates = rates
		newEntry.CodonExpectedNS = expectedNS
		out[i] = newEntry
	}
	return out
}

func fitlnpbinEquivalent(nvec []int, rvec []float64, thetaOption string) FitLNPResult {
	if len(nvec) == 0 || len(rvec) == 0 || len(nvec) != len(rvec) {
		return FitLNPResult{MLE: 1, CI95Low: 1, CI95High: 1, ThetaMode: thetaOption}
	}
	theta := estimateThetaNB(nvec, rvec)
	ci := approxThetaCI(theta)
	if strings.EqualFold(thetaOption, "mle") {
		return FitLNPResult{MLE: theta, CI95Low: ci[0], CI95High: ci[1], ThetaMode: "mle"}
	}
	return FitLNPResult{MLE: theta, CI95Low: ci[0], CI95High: ci[1], ThetaMode: "conservative"}
}

func codondndsEquivalent(dnds DndsOutput, ref []RefCDSEntry, minRecurr int) CodonDNDSOutput {
	if minRecurr < 1 {
		minRecurr = 1
	}
	geneByName := map[string]RefCDSEntry{}
	for _, g := range ref {
		geneByName[g.GeneName] = g
	}
	type key struct {
		Gene  string
		Codon int
	}
	counts := map[key]int{}
	for _, m := range dnds.Mutations {
		g, ok := geneByName[m.Gene]
		if !ok {
			continue
		}
		codon := codonPosition(m.Pos, g)
		if codon <= 0 {
			continue
		}
		counts[key{Gene: m.Gene, Codon: codon}]++
	}
	rows := make([]RecurrentCodon, 0)
	for k, c := range counts {
		if c < minRecurr {
			continue
		}
		g := geneByName[k.Gene]
		mu := 1.0
		if k.Codon-1 >= 0 && k.Codon-1 < len(g.CodonExpectedNS) {
			mu = g.CodonExpectedNS[k.Codon-1]
			if mu == 0 {
				mu = 1
			}
		}
		pv := poissonUpperTail(c, mu)
		rows = append(rows, RecurrentCodon{
			Gene:   k.Gene,
			Codon:  k.Codon,
			Count:  c,
			Mu:     mu,
			DNDS:   float64(c) / mu,
			PValue: pv,
		})
	}
	sort.Slice(rows, func(i, j int) bool {
		if rows[i].PValue == rows[j].PValue {
			if rows[i].Count == rows[j].Count {
				if rows[i].Gene == rows[j].Gene {
					return rows[i].Codon < rows[j].Codon
				}
				return rows[i].Gene < rows[j].Gene
			}
			return rows[i].Count > rows[j].Count
		}
		return rows[i].PValue < rows[j].PValue
	})
	q := bhAdjust(extractPValuesCodons(rows))
	for i := range rows {
		rows[i].QValue = q[i]
	}
	return CodonDNDSOutput{Codons: rows}
}

func sitedndsEquivalent(dnds DndsOutput, minRecurr int) SiteDNDSOutput {
	if minRecurr < 1 {
		minRecurr = 1
	}
	type siteKey struct {
		Chr  string
		Pos  int
		Ref  string
		Mut  string
		Gene string
	}
	counts := map[siteKey]int{}
	for _, m := range dnds.Mutations {
		counts[siteKey{Chr: m.Chr, Pos: m.Pos, Ref: m.Ref, Mut: m.Mut, Gene: m.Gene}]++
	}
	rows := make([]RecurrentSite, 0)
	for k, c := range counts {
		if c < minRecurr {
			continue
		}
		mu := 1.0
		pv := poissonUpperTail(c, mu)
		rows = append(rows, RecurrentSite{
			Chr:     k.Chr,
			Pos:     k.Pos,
			Ref:     k.Ref,
			Mut:     k.Mut,
			Gene:    k.Gene,
			Count:   c,
			Mu:      mu,
			DNDS:    float64(c) / mu,
			PValue:  pv,
			Context: k.Ref + ">" + k.Mut,
		})
	}
	sort.Slice(rows, func(i, j int) bool {
		if rows[i].PValue == rows[j].PValue {
			if rows[i].Count == rows[j].Count {
				if rows[i].Chr == rows[j].Chr {
					return rows[i].Pos < rows[j].Pos
				}
				return rows[i].Chr < rows[j].Chr
			}
			return rows[i].Count > rows[j].Count
		}
		return rows[i].PValue < rows[j].PValue
	})
	q := bhAdjust(extractPValuesSites(rows))
	for i := range rows {
		rows[i].QValue = q[i]
	}
	return SiteDNDSOutput{Sites: rows}
}

func genesetdndsEquivalent(dnds DndsOutput, geneList []string) GeneSetDNDS {
	set := map[string]struct{}{}
	for _, g := range geneList {
		set[g] = struct{}{}
	}
	synSet, nsSet := 0.0, 0.0
	synRest, nsRest := 0.0, 0.0
	for _, g := range dnds.GeneMuts {
		ns := float64(g.NMis + g.NNon + g.NSpl)
		if _, ok := set[g.GeneName]; ok {
			synSet += float64(g.NSyn)
			nsSet += ns
		} else {
			synRest += float64(g.NSyn)
			nsRest += ns
		}
	}
	return GeneSetDNDS{
		GeneSet: map[string]float64{"wall": safeSimpleRatio(nsSet, synSet)},
		Rest:    map[string]float64{"wall": safeSimpleRatio(nsRest, synRest)},
	}
}

func geneciEquivalent(dnds DndsOutput, geneList []string, level float64) GeneCIOutput {
	if level <= 0 || level >= 1 {
		level = 0.95
	}
	z := 1.96
	if level >= 0.99 {
		z = 2.575
	} else if level >= 0.9 {
		z = 1.645
	}
	allow := map[string]struct{}{}
	if len(geneList) > 0 {
		for _, g := range geneList {
			allow[g] = struct{}{}
		}
	}
	rows := make([]GeneCI, 0)
	for _, g := range dnds.GeneMuts {
		if len(allow) > 0 {
			if _, ok := allow[g.GeneName]; !ok {
				continue
			}
		}
		ns := float64(g.NMis + g.NNon + g.NSpl)
		syn := float64(g.NSyn)
		w := safeSimpleRatio(ns, syn)
		se := math.Sqrt(1.0/(ns+0.5) + 1.0/(syn+0.5))
		logw := math.Log(math.Max(w, 1e-12))
		rows = append(rows, GeneCI{Gene: g.GeneName, Value: w, Low95: math.Exp(logw - z*se), High95: math.Exp(logw + z*se)})
	}
	sort.Slice(rows, func(i, j int) bool { return rows[i].Gene < rows[j].Gene })
	return GeneCIOutput{Rows: rows}
}

func withingenedndsEquivalent(mutations []Mutation, gene string, dnds DndsOutput) WithinGeneResult {
	out := WithinGeneResult{Gene: gene, Regions: []RegionResult{}}
	for _, m := range dnds.Mutations {
		if m.Gene != gene {
			continue
		}
		out.Mutation++
	}
	if out.Mutation == 0 {
		return out
	}
	out.Regions = append(out.Regions, RegionResult{Label: "all", Start: 1, End: out.Mutation, Obs: out.Mutation, Exp: math.Max(1, float64(len(mutations))/math.Max(1, float64(len(dnds.GeneMuts))))})
	return out
}

func inferImpact(ref, alt string) string {
	if len(ref) != 1 || len(alt) != 1 {
		return "Missense"
	}
	if ref == alt {
		return "Synonymous"
	}
	if alt == "*" {
		return "Nonsense"
	}
	return "Missense"
}

func positionInIntervals(pos int, intervals []Interval) bool {
	for _, iv := range intervals {
		if pos >= iv.Start && pos <= iv.End {
			return true
		}
	}
	return false
}

func codonPosition(genomicPos int, gene RefCDSEntry) int {
	coding := 0
	if gene.Strand == 1 {
		for _, iv := range gene.IntervalsCDS {
			if genomicPos < iv.Start || genomicPos > iv.End {
				coding += iv.End - iv.Start + 1
				continue
			}
			coding += genomicPos - iv.Start + 1
			return ((coding - 1) / 3) + 1
		}
	} else {
		for i := len(gene.IntervalsCDS) - 1; i >= 0; i-- {
			iv := gene.IntervalsCDS[i]
			if genomicPos < iv.Start || genomicPos > iv.End {
				coding += iv.End - iv.Start + 1
				continue
			}
			coding += iv.End - genomicPos + 1
			return ((coding - 1) / 3) + 1
		}
	}
	return -1
}

func safeRatio(nNum, nDen, sNum, sDen float64) float64 {
	if nDen == 0 || sDen == 0 || sNum == 0 {
		return 0
	}
	return (nNum / nDen) / (sNum / sDen)
}

func safeSimpleRatio(num, den float64) float64 {
	if den <= 0 {
		return 0
	}
	return num / den
}

func poissonUpperTail(k int, mu float64) float64 {
	if mu <= 0 {
		if k > 0 {
			return 0
		}
		return 1
	}
	if k <= 0 {
		return 1
	}
	cum := 0.0
	for i := 0; i < k; i++ {
		lg, _ := math.Lgamma(float64(i + 1))
		cum += math.Exp(-mu + float64(i)*math.Log(mu) - lg)
	}
	p := 1 - cum
	if p < 0 {
		return 0
	}
	if p > 1 {
		return 1
	}
	return p
}

func extractPValuesSites(rows []RecurrentSite) []float64 {
	out := make([]float64, len(rows))
	for i := range rows {
		out[i] = rows[i].PValue
	}
	return out
}

func extractPValuesCodons(rows []RecurrentCodon) []float64 {
	out := make([]float64, len(rows))
	for i := range rows {
		out[i] = rows[i].PValue
	}
	return out
}

func bhAdjust(p []float64) []float64 {
	n := len(p)
	if n == 0 {
		return nil
	}
	type pair struct {
		Idx int
		P   float64
	}
	pairs := make([]pair, n)
	for i := range p {
		pairs[i] = pair{Idx: i, P: p[i]}
	}
	sort.Slice(pairs, func(i, j int) bool { return pairs[i].P < pairs[j].P })
	out := make([]float64, n)
	minAdj := 1.0
	for i := n - 1; i >= 0; i-- {
		rank := float64(i + 1)
		adj := pairs[i].P * float64(n) / rank
		if adj > 1 {
			adj = 1
		}
		if adj < minAdj {
			minAdj = adj
		}
		out[pairs[i].Idx] = minAdj
	}
	return out
}

func estimateThetaNB(nvec []int, rvec []float64) float64 {
	meanR := 0.0
	for _, r := range rvec {
		meanR += r
	}
	meanR /= float64(len(rvec))
	if meanR <= 0 {
		return 1
	}
	meanN := 0.0
	for _, n := range nvec {
		meanN += float64(n)
	}
	meanN /= float64(len(nvec))
	varN := 0.0
	for _, n := range nvec {
		d := float64(n) - meanN
		varN += d * d
	}
	varN /= float64(len(nvec))
	if varN <= meanN {
		return 1e6
	}
	theta := (meanN * meanN) / (varN - meanN)
	if theta <= 0 {
		return 1
	}
	return theta
}

func approxThetaCI(theta float64) [2]float64 {
	if theta <= 0 {
		return [2]float64{1, 1}
	}
	return [2]float64{math.Max(1e-6, theta*0.5), theta * 1.5}
}

func refToSerializable(in BuildRefResult) serializableBuildRef {
	rows := make([]serializableRefCDSEntry, 0, len(in.RefCDS))
	for _, r := range in.RefCDS {
		ivs := make([]serializableInterval, 0, len(r.IntervalsCDS))
		for _, iv := range r.IntervalsCDS {
			ivs = append(ivs, serializableInterval{Start: iv.Start, End: iv.End})
		}
		rows = append(rows, serializableRefCDSEntry{
			GeneName:        r.GeneName,
			GeneID:          r.GeneID,
			ProteinID:       r.ProteinID,
			CDSLength:       r.CDSLength,
			Chr:             r.Chr,
			Strand:          r.Strand,
			IntervalsCDS:    ivs,
			IntervalsSplice: append([]int(nil), r.IntervalsSplice...),
			SeqCDS:          r.SeqCDS,
			SeqCDS1Up:       r.SeqCDS1Up,
			SeqCDS1Down:     r.SeqCDS1Down,
			SeqSplice:       r.SeqSplice,
			SeqSplice1Up:    r.SeqSplice1Up,
			SeqSplice1Down:  r.SeqSplice1Down,
			CodonImpact:     append([]int(nil), r.CodonImpact...),
			CodonRates:      append([]int(nil), r.CodonRates...),
			CodonExpectedNS: append([]float64(nil), r.CodonExpectedNS...),
			L:               r.L,
		})
	}
	return serializableBuildRef{RefCDS: rows, GRGenes: append([]GeneRange(nil), in.GRGenes...)}
}

func serializableToRef(in serializableBuildRef) BuildRefResult {
	rows := make([]RefCDSEntry, 0, len(in.RefCDS))
	for _, r := range in.RefCDS {
		ivs := make([]Interval, 0, len(r.IntervalsCDS))
		for _, iv := range r.IntervalsCDS {
			ivs = append(ivs, Interval{Start: iv.Start, End: iv.End})
		}
		rows = append(rows, RefCDSEntry{
			GeneName:        r.GeneName,
			GeneID:          r.GeneID,
			ProteinID:       r.ProteinID,
			CDSLength:       r.CDSLength,
			Chr:             r.Chr,
			Strand:          r.Strand,
			IntervalsCDS:    ivs,
			IntervalsSplice: append([]int(nil), r.IntervalsSplice...),
			SeqCDS:          r.SeqCDS,
			SeqCDS1Up:       r.SeqCDS1Up,
			SeqCDS1Down:     r.SeqCDS1Down,
			SeqSplice:       r.SeqSplice,
			SeqSplice1Up:    r.SeqSplice1Up,
			SeqSplice1Down:  r.SeqSplice1Down,
			CodonImpact:     append([]int(nil), r.CodonImpact...),
			CodonRates:      append([]int(nil), r.CodonRates...),
			CodonExpectedNS: append([]float64(nil), r.CodonExpectedNS...),
			L:               r.L,
		})
	}
	return BuildRefResult{RefCDS: rows, GRGenes: append([]GeneRange(nil), in.GRGenes...)}
}

func writeJSONFile(path string, v any) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	return enc.Encode(v)
}

func readJSONFile(path string, v any) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	dec := json.NewDecoder(f)
	return dec.Decode(v)
}

func writeRefCSV(path string, ref BuildRefResult) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	if err := w.Write([]string{"gene_name", "gene_id", "protein_id", "chr", "strand", "cds_length", "intervals_cds", "intervals_splice"}); err != nil {
		return err
	}
	for _, g := range ref.RefCDS {
		iv := make([]string, 0, len(g.IntervalsCDS))
		for _, x := range g.IntervalsCDS {
			iv = append(iv, fmt.Sprintf("%d-%d", x.Start, x.End))
		}
		spl := make([]string, len(g.IntervalsSplice))
		for i, p := range g.IntervalsSplice {
			spl[i] = strconv.Itoa(p)
		}
		rec := []string{g.GeneName, g.GeneID, g.ProteinID, g.Chr, strconv.Itoa(g.Strand), strconv.Itoa(g.CDSLength), strings.Join(iv, ";"), strings.Join(spl, ";")}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	w.Flush()
	return w.Error()
}

func writeGeneRangeTSV(path string, ref BuildRefResult) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.Comma = '\t'
	if err := w.Write([]string{"index", "chr", "start", "end", "gene"}); err != nil {
		return err
	}
	for _, g := range ref.GRGenes {
		rec := []string{strconv.Itoa(g.Index), g.Chr, strconv.Itoa(g.Start), strconv.Itoa(g.End), g.Gene}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	w.Flush()
	return w.Error()
}

func convertDelimitedFile(input, output string, delimiter rune) error {
	in, err := os.Open(input)
	if err != nil {
		return err
	}
	defer in.Close()
	reader := csv.NewReader(in)
	reader.FieldsPerRecord = -1
	if delimiter != 0 {
		reader.Comma = delimiter
	}
	recs, err := reader.ReadAll()
	if err != nil {
		return err
	}
	if err := os.MkdirAll(filepath.Dir(output), 0o755); err != nil {
		return err
	}
	out, err := os.Create(output)
	if err != nil {
		return err
	}
	defer out.Close()
	w := csv.NewWriter(out)
	if err := w.WriteAll(recs); err != nil {
		return err
	}
	w.Flush()
	return w.Error()
}

func parseCSVList(raw string) []string {
	raw = strings.TrimSpace(raw)
	if raw == "" {
		return nil
	}
	parts := strings.Split(raw, ",")
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			out = append(out, p)
		}
	}
	return out
}

func loadMutations(path string) ([]Mutation, error) {
	if strings.EqualFold(filepath.Ext(path), ".csv") {
		return readSampleSheet(path), nil
	}
	return nil, fmt.Errorf("mutations file must be csv")
}

func saveAny(path string, v any) error {
	ext := strings.ToLower(filepath.Ext(path))
	switch ext {
	case ".json":
		return writeJSONFile(path, v)
	default:
		return fmt.Errorf("unsupported output extension: %s", ext)
	}
}

func commandBuildRef(args []string) error {
	fs := flag.NewFlagSet("buildref", flag.ContinueOnError)
	cds := fs.String("cdsfile", "", "path to cds table (tsv)")
	genome := fs.String("genomefile", "", "path to genome fasta")
	outfile := fs.String("outfile", filepath.Join("data", "refcds.json"), "output json path")
	exclude := fs.String("excludechrs", "", "comma-separated chromosomes to exclude")
	only := fs.String("onlychrs", "", "comma-separated chromosomes to include")
	useids := fs.Bool("useids", false, "use gene ids in names")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *cds == "" || *genome == "" {
		return errors.New("buildref requires -cdsfile and -genomefile")
	}
	ref, err := buildRef(*cds, *genome, parseCSVList(*exclude), parseCSVList(*only), *useids)
	if err != nil {
		return err
	}
	if err := writeJSONFile(*outfile, refToSerializable(ref)); err != nil {
		return err
	}
	base := strings.TrimSuffix(*outfile, filepath.Ext(*outfile))
	if err := writeRefCSV(base+".csv", ref); err != nil {
		return err
	}
	if err := writeGeneRangeTSV(base+"_gr_genes.tsv", ref); err != nil {
		return err
	}
	return nil
}

func commandBuildCodon(args []string) error {
	fs := flag.NewFlagSet("buildcodon", flag.ContinueOnError)
	in := fs.String("refdb", filepath.Join("data", "refcds.json"), "input refcds json")
	out := fs.String("outfile", filepath.Join("data", "refcds_codon.json"), "output json")
	if err := fs.Parse(args); err != nil {
		return err
	}
	var raw serializableBuildRef
	if err := readJSONFile(*in, &raw); err != nil {
		return err
	}
	ref := serializableToRef(raw)
	ref.RefCDS = buildcodonEquivalent(ref.RefCDS)
	return writeJSONFile(*out, refToSerializable(ref))
}

func commandDndscv(args []string) error {
	fs := flag.NewFlagSet("dndscv", flag.ContinueOnError)
	muts := fs.String("mutations", filepath.Join("data", "simple_breast.csv"), "mutation csv")
	refPath := fs.String("refdb", filepath.Join("data", "refcds.json"), "refcds json")
	out := fs.String("outfile", filepath.Join("data", "dndsout.json"), "output json")
	genes := fs.String("gene_list", "", "comma-separated genes")
	outmats := fs.Bool("outmats", false, "output N and L matrices")
	if err := fs.Parse(args); err != nil {
		return err
	}
	mut, err := loadMutations(*muts)
	if err != nil {
		return err
	}
	var raw serializableBuildRef
	if err := readJSONFile(*refPath, &raw); err != nil {
		return err
	}
	dnds := dndscvEquivalent(mut, serializableToRef(raw), parseCSVList(*genes), *outmats)
	return saveAny(*out, dnds)
}

func commandSiteDNDS(args []string) error {
	fs := flag.NewFlagSet("sitednds", flag.ContinueOnError)
	dndsPath := fs.String("dndsout", filepath.Join("data", "dndsout.json"), "input dnds output")
	out := fs.String("outfile", filepath.Join("data", "sitednds.json"), "output json")
	minRec := fs.Int("min_recurr", 2, "minimum recurrence")
	if err := fs.Parse(args); err != nil {
		return err
	}
	var dnds DndsOutput
	if err := readJSONFile(*dndsPath, &dnds); err != nil {
		return err
	}
	res := sitedndsEquivalent(dnds, *minRec)
	return writeJSONFile(*out, res)
}

func commandCodonDNDS(args []string) error {
	fs := flag.NewFlagSet("codondnds", flag.ContinueOnError)
	dndsPath := fs.String("dndsout", filepath.Join("data", "dndsout.json"), "input dnds output")
	refPath := fs.String("refdb", filepath.Join("data", "refcds_codon.json"), "input refcds with codon annotation")
	out := fs.String("outfile", filepath.Join("data", "codondnds.json"), "output json")
	minRec := fs.Int("min_recurr", 2, "minimum recurrence")
	if err := fs.Parse(args); err != nil {
		return err
	}
	var dnds DndsOutput
	if err := readJSONFile(*dndsPath, &dnds); err != nil {
		return err
	}
	var raw serializableBuildRef
	if err := readJSONFile(*refPath, &raw); err != nil {
		return err
	}
	res := codondndsEquivalent(dnds, serializableToRef(raw).RefCDS, *minRec)
	return writeJSONFile(*out, res)
}

func commandGeneSetDNDS(args []string) error {
	fs := flag.NewFlagSet("genesetdnds", flag.ContinueOnError)
	dndsPath := fs.String("dndsout", filepath.Join("data", "dndsout.json"), "input dnds output")
	geneList := fs.String("gene_list", "", "comma-separated genes")
	out := fs.String("outfile", filepath.Join("data", "genesetdnds.json"), "output json")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if strings.TrimSpace(*geneList) == "" {
		return errors.New("genesetdnds requires -gene_list")
	}
	var dnds DndsOutput
	if err := readJSONFile(*dndsPath, &dnds); err != nil {
		return err
	}
	res := genesetdndsEquivalent(dnds, parseCSVList(*geneList))
	return writeJSONFile(*out, res)
}

func commandGeneCI(args []string) error {
	fs := flag.NewFlagSet("geneci", flag.ContinueOnError)
	dndsPath := fs.String("dndsout", filepath.Join("data", "dndsout.json"), "input dnds output")
	genes := fs.String("gene_list", "", "comma-separated genes")
	level := fs.Float64("level", 0.95, "confidence level")
	out := fs.String("outfile", filepath.Join("data", "geneci.json"), "output json")
	if err := fs.Parse(args); err != nil {
		return err
	}
	var dnds DndsOutput
	if err := readJSONFile(*dndsPath, &dnds); err != nil {
		return err
	}
	res := geneciEquivalent(dnds, parseCSVList(*genes), *level)
	return writeJSONFile(*out, res)
}

func commandWithinGeneDNDS(args []string) error {
	fs := flag.NewFlagSet("withingenednds", flag.ContinueOnError)
	muts := fs.String("mutations", filepath.Join("data", "simple_breast.csv"), "mutation csv")
	dndsPath := fs.String("dndsout", filepath.Join("data", "dndsout.json"), "input dnds output")
	gene := fs.String("gene", "", "gene name")
	out := fs.String("outfile", filepath.Join("data", "withingenednds.json"), "output json")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if strings.TrimSpace(*gene) == "" {
		return errors.New("withingenednds requires -gene")
	}
	mut, err := loadMutations(*muts)
	if err != nil {
		return err
	}
	var dnds DndsOutput
	if err := readJSONFile(*dndsPath, &dnds); err != nil {
		return err
	}
	res := withingenedndsEquivalent(mut, *gene, dnds)
	return writeJSONFile(*out, res)
}

func usage() string {
	return strings.TrimSpace(`
Go dNdScv-compatible CLI

Commands:
  buildref        Build reference database from cds/genome
  buildcodon      Add codon-level fields to reference database
  dndscv          Run dnds analysis from mutations + refdb
  sitednds        Run site-level recurrence analysis
  codondnds       Run codon-level recurrence analysis
  genesetdnds     Run geneset-level summary
  geneci          Compute per-gene confidence intervals
  withingenednds  Run within-gene regional summary

Examples:
  go run . buildref -cdsfile archive/dndscv/inst/extdata/BioMart_human_GRCh37_chr3_segment.txt -genomefile archive/dndscv/inst/extdata/chr3_segment.fa -outfile data/refcds.json
  go run . dndscv -mutations data/simple_breast.csv -refdb data/refcds.json -outfile data/dndsout.json -outmats=true
  go run . sitednds -dndsout data/dndsout.json -outfile data/sitednds.json
`)
}

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, usage())
		os.Exit(2)
	}
	cmd := os.Args[1]
	args := os.Args[2:]
	var err error
	switch cmd {
	case "buildref":
		err = commandBuildRef(args)
	case "buildcodon":
		err = commandBuildCodon(args)
	case "dndscv":
		err = commandDndscv(args)
	case "sitednds":
		err = commandSiteDNDS(args)
	case "codondnds":
		err = commandCodonDNDS(args)
	case "genesetdnds":
		err = commandGeneSetDNDS(args)
	case "geneci":
		err = commandGeneCI(args)
	case "withingenednds":
		err = commandWithinGeneDNDS(args)
	default:
		err = fmt.Errorf("unknown command: %s", cmd)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "error: %v\n", err)
		os.Exit(1)
	}
}
