package main

// Program to compute a distance matrix for compositional data (proportions).

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

func main() {
	mtxp := flag.String("f", "dat.csv", "data matrix, rows = samples, cols = species")

	flag.Parse()

	matrix := *mtxp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	dat := ReadCsvMatrix64(file)

	// convert data to proportions
	Dat2Prop(dat)

	dat.WriteCSV()

}
