package main

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"math/rand"
)

// Create similarity/distance matrix from random numbers.
func main() {
	var (
		fmtp  = flag.String("f", "csv", "output format, either csv or dat")
		rowsp = flag.Int("r", 100, "number of rows/columns")
		seedp = flag.Int64("s", 0, "seed to random number generator")
	)

	flag.Parse()

	rows := *rowsp
	outfmt := *fmtp
	seed := *seedp
	rand.Seed(seed)

	a := NewMatrix64(rows, rows)

	for i := 0; i < rows-1; i++ {
		for j := i + 1; j < rows; j++ {
			a[i][j] = rand.Float64()
			a[j][i] = a[i][j]
		}
	}

	for i := 0; i < rows; i++ {
		a[i][i] = 1.0
	}

	switch outfmt {
	case "dat":
		a.Print()
	case "csv":
		a.WriteCSV()
	default:
		a.WriteCSV()
	}
}
