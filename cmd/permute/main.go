package main

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

// Permute similarity/distance matrix using supplied permutation vector.
func main() {
	var (
		mtxp = flag.String("m", "R50.csv", "square matrix to be permuted")
		vecp = flag.String("p", "perm.csv", "permutation vector")
		fmtp = flag.String("f", "csv", "output format, either csv or dat")
	)

	flag.Parse()

	matrix := *mtxp
	vector := *vecp

	outfmt := *fmtp

	file, err := os.Open(matrix) // similarity matrix
	if err != nil {
		panic("file with matrix does not exist")
	}

	a := ReadCsvMatrix64(file)
	file2, err := os.Open(vector) // permutation vector
	if err != nil {
		panic("file with permutation does not exist")
	}
	p := ReadCsvIntVector(file2)
	if a.Rows() != p.Len() {
		panic("dimensions do not match")
	}

	a.Permute(p, p)
	switch outfmt {
	case "dat":
		a.Print()
	case "csv":
		a.WriteCSV()
	default:
		a.WriteCSV()
	}
}
