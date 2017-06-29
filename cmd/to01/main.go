package main

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

// Transforms matrix to 0 - 1 according to given threshold.
func main() {
	var (
		mtxp = flag.String("f", "R50.csv", "matrix to be transformed")
		tp   = flag.Float64("t", 0.5, "threshold")
	)

	flag.Parse()

	matrix := *mtxp
	t := *tp

	file, err := os.Open(matrix) // similarity matrix
	if err != nil {
		panic("file with matrix does not exist")
	}

	a := ReadCsvMatrix64(file)
	a.ToZeroOne(t)
	a.WriteCSV()
}
