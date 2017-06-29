package main

// Program to calculate similarity matrix from proportional (compositional) data matrix using logratio and Kendall's 1971 Theorem I

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

func main() {
	var (
		mtxp = flag.String("m", "Rogers2011Environmental-counts.csv", "data matrix to be turned to similarity one")
		colp = flag.Int("c", -1, "column to be used as the logratio denominator")
	)

	flag.Parse()

	matrix := *mtxp
	col := *colp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)

	// recalc data to logratios
	Logratio(a, col)

	// compute similarity matrix using Kendall's 1971 Theorem I
	sim := Adj2Sim(a, "rows")
	sim.WriteCSV()
}
