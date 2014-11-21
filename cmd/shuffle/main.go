package main

import (
	. "code.google.com/p/seriation"
	"flag"
	"math/rand"
	"os"
)

// Permute similarity/distance matrix at random.
func main() {
	var (
		mtxp = flag.String("m", "R50.csv", "square matrix to be permuted")
		fmtp = flag.String("f", "csv", "output format, either csv or dat")
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
	nSamp := a.Rows()

	p := NewIntVector(nSamp)
	p.Order()
	//	fmt.Println("===============Permuted mtx========================================")
	for i := 0; i < 1; i++ {
		p.Perm()
		a.Permute(p, p)
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
