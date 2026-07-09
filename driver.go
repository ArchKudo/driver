package main

import (
	"fmt"
	"os"
	// "errors"
	"encoding/csv"
	// "io"
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

func main() {

	file := "data/simple_breast.csv"

	sample_sheet := readSampleSheet(file)

	for i := 0; i < 5; i++ {
		fmt.Println(sample_sheet[i])
	}
}
