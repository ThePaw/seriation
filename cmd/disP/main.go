package main

// Program to compute a distance matrix for compositional data (proportions).

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

func main() {
	var (
		dis Matrix64
	)
	mtxp := flag.String("f", "dat.csv", "data matrix, rows = samples, cols = species")
	disfp := flag.String("i", "ManlyOverlap", "distance function")
	forcep := flag.Bool("F", false, "force distances to [0..1]?")

	flag.Parse()

	matrix := *mtxp
	disf := *disfp
	force := *forcep

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	dat := ReadCsvMatrix64(file)

	// convert data to proportions
	Dat2Prop(dat)

	// distance functions
	switch disf {
	case "Edwards":
		dis = Edwards(dat)
	case "Manly":
		dis = Manly(dat)
	case "ManlyOverlap":
		dis = ManlyOverlap(dat)
	case "Nei":
		dis = Nei(dat)
	case "Rogers72":
		dis = Rogers72(dat)
	case "Reynolds":
		dis = Reynolds(dat)
	case "Sforza":
		dis = Sforza(dat)
	case "Nei83":
		dis = Nei83(dat)
	case "Euclid":
		dis = Euclid(dat)
	case "Percentage":
		dis = Percentage(dat)
	case "Nei73":
		dis = Nei73(dat)
	case "PercentageS":
		dis = PercentageS(dat)
		dis.SimToDist()
	case "Kendall71S":
		dis = Kendall71S(dat)
		dis.SimToDist()
	default:
		panic("index not available")
	}
if force {
dis.ForceTo01()
}
	dis.WriteCSV()

}
