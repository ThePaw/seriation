package main

// Program to compute a similarity matrix for binary data.

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"os"
)

func main() {
	var (
		sim Matrix64
	)
	mtxp := flag.String("f", "dat.csv", "data matrix, rows = samples, cols = species")
	simfp := flag.String("i", "Soerensen", "similarity index")
	normp := flag.Bool("n", false, "normalize a, b, c, when applicable")

	flag.Parse()

	matrix := *mtxp
	simf := *simfp
	norm := *normp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	dat := ReadCsvMatrix64(file)
	//	nSamp := a.Rows()

	// indices
	switch simf {
	// asymmetric
	case "CoCoGaston":
		sim = CoCoGaston(dat, norm)
	case "Cody":
		sim = Cody(dat, norm)
	case "Dice":
		sim = Dice(dat, norm)
	case "Fager":
		sim = Fager(dat, norm)
	case "Harrison":
		sim = Harrison(dat, norm)
	case "Harte":
		sim = Harte(dat, norm)
	case "Jaccard":
		sim = Jaccard(dat, norm)
	case "Johnson1":
		sim = Johnson1(dat, norm)
	case "Johnson2":
		sim = Johnson2(dat, norm)
	case "Kulczinsky1":
		sim = Kulczinsky1(dat, norm)
	case "Kulczinsky2":
		sim = Kulczinsky2(dat, norm)
	case "Lamont":
		sim = Lamont(dat, norm)
	case "Lande":
		sim = Lande(dat, norm)
	case "Legendre2":
		sim = Legendre2(dat, norm)
	case "Lennon1":
		sim = Lennon1(dat, norm)
	case "Lennon2":
		sim = Lennon2(dat, norm)
	case "Maarel":
		sim = Maarel(dat, norm)
	case "Magurran":
		sim = Magurran(dat, norm)
	case "McConnagh":
		sim = McConnagh(dat, norm)
	case "Mountford":
		sim = Mountford(dat, norm)
	case "Ochiai":
		sim = Ochiai(dat, norm)
	case "Routledge1":
		sim = Routledge1(dat, norm)
	case "Routledge2":
		sim = Routledge2(dat, norm)
	case "Routledge3":
		sim = Routledge3(dat, norm)
	case "Ruggiero":
		sim = Ruggiero(dat, norm)
	case "Simpson1":
		sim = Simpson1(dat, norm)
	case "Simpson2":
		sim = Simpson2(dat, norm)
	case "Soerensen":
		sim = Soerensen(dat, norm)
	case "Sokal1":
		sim = Sokal1(dat, norm)
	case "Sorgenfrei":
		sim = Sorgenfrei(dat, norm)
	case "Weiher":
		sim = Weiher(dat, norm)
	case "Whittaker":
		sim = Whittaker(dat, norm)
	case "Williams1":
		sim = Williams1(dat, norm)
	case "Williams2":
		sim = Williams2(dat, norm)
	case "WilsonShmida":
		sim = WilsonShmida(dat, norm)
		// symmetric
	case "Baroni":
		sim = Baroni(dat)
	case "BinEuclidean":
		sim = BinEuclidean(dat)
	case "ChiSquare":
		sim = ChiSquare(dat)
	case "Dennis":
		sim = Dennis(dat)
	case "Ellis":
		sim = Ellis(dat)
	case "Eyraud":
		sim = Eyraud(dat)
	case "Forbes":
		sim = Forbes(dat)
	case "Fossum":
		sim = Fossum(dat)
	case "Gower":
		sim = Gower(dat)
	case "Hamann":
		sim = Hamann(dat)
	case "Legendre":
		sim = Legendre(dat)
	case "Manhattan":
		sim = Manhattan(dat)
	case "Margaleff":
		sim = Margaleff(dat)
	case "Peirce":
		sim = Peirce(dat)
	case "Rogers":
		sim = Rogers(dat)
	case "SimpleMatching":
		sim = SimpleMatching(dat)
	case "Sokal2":
		sim = Sokal2(dat)
	case "Sokal3":
		sim = Sokal3(dat)
	case "Sokal4":
		sim = Sokal4(dat)
	case "Stiles":
		sim = Stiles(dat)
	case "Yule1":
		sim = Yule1(dat)
	case "Yule2":
		sim = Yule2(dat)
	default:
		panic("index not available")
	}

	sim.WriteCSV()

}
