package main

// Program to calculate similarity matrix from data matrix using Kendall's 1971 Theorem I

import (
	. "code.google.com/p/seriation"
	"flag"
	"fmt"
	"os"
)

func main() {
	var (
		mtxp  = flag.String("m", "Rogers2011Environmental-counts.csv", "data matrix to be turned to similarity one")
	)

	flag.Parse()

	matrix := *mtxp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)
//	nSamp := a.Rows()

	// compute similarity matrix using Kendall's 1971 Theorem I
	sim := Adj2Sim(a, "rows")
	sim.WriteCSV()
fmt.Println()
sim.Print()
}
