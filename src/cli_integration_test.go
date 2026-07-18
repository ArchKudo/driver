package main

import (
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

func TestBuildcodonEquivalentAddsFields(t *testing.T) {
	entry := RefCDSEntry{
		GeneName:    "G",
		GeneID:      "id",
		ProteinID:   "p",
		CDSLength:   6,
		SeqCDS:      "ATGGAA",
		SeqCDS1Up:   "TATGGA",
		SeqCDS1Down: "TGGAAT",
		L:           [192][4]int{},
	}
	out := buildcodonEquivalent([]RefCDSEntry{entry})
	if len(out) != 1 {
		t.Fatalf("unexpected output length: %d", len(out))
	}
	if len(out[0].CodonImpact) == 0 {
		t.Fatalf("expected codon impact to be populated")
	}
	if len(out[0].CodonRates) == 0 {
		t.Fatalf("expected codon rates to be populated")
	}
	if len(out[0].CodonExpectedNS) != 2 {
		t.Fatalf("expected codon expected NS length 2, got %d", len(out[0].CodonExpectedNS))
	}
}

func TestDndscvSiteAndCodonFlow(t *testing.T) {
	ref := BuildRefResult{RefCDS: []RefCDSEntry{{
		GeneName:     "GENE1",
		GeneID:       "G1",
		ProteinID:    "P1",
		Chr:          "1",
		Strand:       1,
		CDSLength:    6,
		IntervalsCDS: []Interval{{Start: 100, End: 105}},
		SeqCDS:       "ATGGAA",
		SeqCDS1Up:    "TATGGA",
		SeqCDS1Down:  "TGGAAT",
	}}}
	ref.RefCDS = buildcodonEquivalent(ref.RefCDS)
	muts := []Mutation{
		{sampleID: "S1", chr: "1", pos: "100", ref: "A", alt: "C"},
		{sampleID: "S2", chr: "1", pos: "100", ref: "A", alt: "C"},
		{sampleID: "S3", chr: "1", pos: "103", ref: "G", alt: "A"},
	}
	dnds := dndscvEquivalent(muts, ref, nil, true)
	if len(dnds.Mutations) != 3 {
		t.Fatalf("expected 3 annotated mutations, got %d", len(dnds.Mutations))
	}
	sites := sitedndsEquivalent(dnds, 2)
	if len(sites.Sites) != 1 {
		t.Fatalf("expected one recurrent site, got %d", len(sites.Sites))
	}
	codons := codondndsEquivalent(dnds, ref.RefCDS, 2)
	if len(codons.Codons) != 1 {
		t.Fatalf("expected one recurrent codon, got %d", len(codons.Codons))
	}
}

func TestCLICommandsEndToEnd(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping cli integration in short mode")
	}

	tmp := t.TempDir()
	root := filepath.Clean("..")

	cdsSrc := filepath.Join(root, "archive", "dndscv", "inst", "extdata", "BioMart_human_GRCh37_chr3_segment.txt")
	faSrc := filepath.Join(root, "archive", "dndscv", "inst", "extdata", "chr3_segment.fa")
	faIdxSrc := filepath.Join(root, "archive", "dndscv", "inst", "extdata", "chr3_segment.fa.fai")

	cdsDst := filepath.Join(tmp, "cds.tsv")
	faDst := filepath.Join(tmp, "genome.fa")
	mutDst := filepath.Join(tmp, "muts.csv")
	copyMust(t, cdsSrc, cdsDst)
	copyMust(t, faSrc, faDst)
	copyMust(t, faIdxSrc, faDst+".fai")

	refJSON := filepath.Join(tmp, "refcds.json")
	refCodonJSON := filepath.Join(tmp, "refcds_codon.json")
	dndsJSON := filepath.Join(tmp, "dndsout.json")
	siteJSON := filepath.Join(tmp, "sitednds.json")
	codonJSON := filepath.Join(tmp, "codondnds.json")
	geneSetJSON := filepath.Join(tmp, "geneset.json")
	geneCIJSON := filepath.Join(tmp, "geneci.json")
	withinJSON := filepath.Join(tmp, "within.json")

	runMainCmd(t, root, "buildref", "-cdsfile", cdsDst, "-genomefile", faDst, "-outfile", refJSON)
	writeMutationCSVFromRef(t, refJSON, mutDst)
	runMainCmd(t, root, "buildcodon", "-refdb", refJSON, "-outfile", refCodonJSON)
	runMainCmd(t, root, "dndscv", "-mutations", mutDst, "-refdb", refJSON, "-outfile", dndsJSON, "-outmats=true")
	runMainCmd(t, root, "sitednds", "-dndsout", dndsJSON, "-outfile", siteJSON)
	runMainCmd(t, root, "codondnds", "-dndsout", dndsJSON, "-refdb", refCodonJSON, "-outfile", codonJSON)

	var dnds DndsOutput
	readJSONMust(t, dndsJSON, &dnds)
	if len(dnds.GeneMuts) == 0 {
		t.Fatalf("expected dndscv gene summaries")
	}
	gene := dnds.GeneMuts[0].GeneName
	runMainCmd(t, root, "genesetdnds", "-dndsout", dndsJSON, "-gene_list", gene, "-outfile", geneSetJSON)
	runMainCmd(t, root, "geneci", "-dndsout", dndsJSON, "-gene_list", gene, "-outfile", geneCIJSON)
	runMainCmd(t, root, "withingenednds", "-mutations", mutDst, "-dndsout", dndsJSON, "-gene", gene, "-outfile", withinJSON)

	mustExist(t, refJSON)
	mustExist(t, strings.TrimSuffix(refJSON, ".json")+".csv")
	mustExist(t, strings.TrimSuffix(refJSON, ".json")+"_gr_genes.tsv")
	mustExist(t, refCodonJSON)
	mustExist(t, dndsJSON)
	mustExist(t, siteJSON)
	mustExist(t, codonJSON)
	mustExist(t, geneSetJSON)
	mustExist(t, geneCIJSON)
	mustExist(t, withinJSON)
}

func writeMutationCSVFromRef(t *testing.T, refPath, mutPath string) {
	t.Helper()
	var raw serializableBuildRef
	readJSONMust(t, refPath, &raw)
	ref := serializableToRef(raw)
	if len(ref.RefCDS) == 0 || len(ref.RefCDS[0].IntervalsCDS) == 0 {
		t.Fatalf("reference database is empty")
	}

	entry := ref.RefCDS[0]
	chr := entry.Chr
	pos1 := entry.IntervalsCDS[0].Start
	pos2 := pos1
	if entry.IntervalsCDS[0].End > entry.IntervalsCDS[0].Start {
		pos2 = pos1 + 1
	}

	refBase := "A"
	if len(entry.SeqCDS) > 0 {
		refBase = strings.ToUpper(string(entry.SeqCDS[0]))
	}
	altBase := "C"
	if refBase == "C" {
		altBase = "A"
	}

	content := strings.Join([]string{
		"sampleID,chr,pos,ref,alt",
		"S1," + chr + "," + strconv.Itoa(pos1) + "," + refBase + "," + altBase,
		"S2," + chr + "," + strconv.Itoa(pos1) + "," + refBase + "," + altBase,
		"S3," + chr + "," + strconv.Itoa(pos2) + "," + refBase + "," + altBase,
	}, "\n") + "\n"

	if err := os.WriteFile(mutPath, []byte(content), 0o644); err != nil {
		t.Fatalf("failed to write mutations csv %s: %v", mutPath, err)
	}
}

func runMainCmd(t *testing.T, root string, args ...string) {
	t.Helper()
	cmd := exec.Command("go", append([]string{"run", "./src"}, args...)...)
	cmd.Dir = root
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("command failed (%v): %v\n%s", args, err, string(out))
	}
}

func copyMust(t *testing.T, src, dst string) {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("failed to read source %s: %v", src, err)
	}
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatalf("failed to write destination %s: %v", dst, err)
	}
}

func mustExist(t *testing.T, p string) {
	t.Helper()
	if _, err := os.Stat(p); err != nil {
		t.Fatalf("expected file %s to exist: %v", p, err)
	}
}

func readJSONMust(t *testing.T, p string, v any) {
	t.Helper()
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("failed to read %s: %v", p, err)
	}
	if err := json.Unmarshal(b, v); err != nil {
		t.Fatalf("failed to decode %s: %v", p, err)
	}
}
