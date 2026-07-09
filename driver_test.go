package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeTestCSV(t *testing.T, content string) string {
	t.Helper()

	dir := t.TempDir()
	path := filepath.Join(dir, "samples.csv")

	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatalf("failed to create test csv: %v", err)
	}

	return path
}

func assertPanicsWith(t *testing.T, expectedSubstring string, fn func()) {
	t.Helper()

	defer func() {
		r := recover()
		if r == nil {
			t.Fatalf("expected panic containing %q but no panic occurred", expectedSubstring)
		}

		msg := ""
		switch v := r.(type) {
		case error:
			msg = v.Error()
		default:
			msg = "<non-error panic>"
		}

		if !strings.Contains(msg, expectedSubstring) {
			t.Fatalf("panic message mismatch: got %q, want substring %q", msg, expectedSubstring)
		}
	}()

	fn()
}

func TestReadSampleSheet_Success(t *testing.T) {
	path := writeTestCSV(t, "sample_id,chr,pos,ref,alt\ns1,1,100,A,T\ns2,X,200,G,C\n")

	got := readSampleSheet(path)

	if len(got) != 2 {
		t.Fatalf("unexpected number of mutations: got %d, want 2", len(got))
	}

	wantFirst := Mutation{sampleID: "s1", chr: "1", pos: "100", ref: "A", alt: "T"}
	if got[0] != wantFirst {
		t.Fatalf("unexpected first mutation: got %#v, want %#v", got[0], wantFirst)
	}

	wantSecond := Mutation{sampleID: "s2", chr: "X", pos: "200", ref: "G", alt: "C"}
	if got[1] != wantSecond {
		t.Fatalf("unexpected second mutation: got %#v, want %#v", got[1], wantSecond)
	}
}

func TestReadSampleSheet_PanicsWhenFileMissing(t *testing.T) {
	assertPanicsWith(t, "Couldn't load samplesheet", func() {
		_ = readSampleSheet(filepath.Join(t.TempDir(), "missing.csv"))
	})
}

func TestReadSampleSheet_PanicsOnNoRecords(t *testing.T) {
	path := writeTestCSV(t, "sample_id,chr,pos,ref,alt\n")

	assertPanicsWith(t, "No records found", func() {
		_ = readSampleSheet(path)
	})
}

// Go already checks when number of header fields don't match fields in records
func TestReadSampleSheet_PanicsOnWrongColumnCount(t *testing.T) {
	path := writeTestCSV(t, "sample_id,chr,pos,ref,alt,type\ns1,1,100,A,C,tumour\n")

	assertPanicsWith(t, "expected 5", func() {
		_ = readSampleSheet(path)
	})
}
