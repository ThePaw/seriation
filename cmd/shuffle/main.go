package main

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"math/rand"
	"os"
)

// Permute similarity/distance matrix at random.
func main() {
	var (
		mtxp  = flag.String("f", "R50.csv", "square matrix to be permuted")
		fmtp  = flag.String("F", "csv", "output format, either csv or dat")
		seedp = flag.Int64("s", 0, "seed to random number generator")
	)

	flag.Parse()

	matrix := *mtxp
	outfmt := *fmtp
	seed := *seedp
	rand.Seed(seed)

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)
	rows, _ := a.Dims()

	p := NewIntVector(rows)
	//	q := NewIntVector(cols)
	p.Perm()
	//	q.Order()
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
