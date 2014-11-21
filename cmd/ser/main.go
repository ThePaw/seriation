package main

// Program to seriate a square symmetric similarity matrix.

import (
	. "code.google.com/p/seriation"
	"flag"
	"fmt"
	"math/rand"
	"os"
)

func main() {
	var (
		optMethod, optMethodForImpro OptMethod3
		objFn                        ObjFn
		isLoss, isDistFn             bool
		mtxp  = flag.String("m", "R50.csv", "similarity matrix to be seriated")
		objfp = flag.String("o", "H", "objective function")
		iterp = flag.Int("i", 1, "number of iterations")
		seedp = flag.Int64("s", 33, "seed")
	)

	// ============= Set parameters here: =======================

	method := "RobSA3"
	impro := 0
	window := 0
	improMethod := "RobSA3"

	// =======================================================

	flag.Parse()

	matrix := *mtxp
	objf := *objfp
	nIter := *iterp
	seed := *seedp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)
	nSamp := a.Rows()

	p := NewIntVector(nSamp)     
	p.Order()
	p.Perm()

	switch method {
	case "RobSA3":
		optMethod = RobSA3
	case "RobFA3":
		optMethod = RobFA3
	default:
		panic("method: unknown")
	}

	switch improMethod {
	case "RobSA3":
		optMethodForImpro = RobSA3
	case "RobFA3":
		optMethodForImpro = RobFA3
	default:
		panic("method: unknown")
	}

	switch objf {
	// distance gain
	case "Wrug":
		objFn = Wrug
		isLoss = false
		isDistFn = true
	case "Wrcug":
		objFn = Wrcug
		isLoss = false
		isDistFn = true
	case "Wrwg":
		objFn = Wrwg
		isLoss = false
		isDistFn = true
	case "Wrcwg":
		objFn = Wrcwg
		isLoss = false
		isDistFn = true
	case "H":
		objFn = H
		isLoss = false
		isDistFn = true
	case "Ine":
		objFn = Ine
		isLoss = false
		isDistFn = true

		// distance loss
	case "Lsq":
		objFn = Lsq
		isLoss = true
		isDistFn = true
	case "Ham":
		objFn = Ham
		isLoss = true
		isDistFn = true
	case "Are":
		objFn = Are
		isLoss = true
		isDistFn = true
	case "Ware":
		objFn = Ware
		isLoss = true
		isDistFn = true
	case "Dware":
		objFn = Dware
		isLoss = true
		isDistFn = true
	case "Nsd":
		objFn = Nsd
		isLoss = true
		isDistFn = true
	case "Msd":
		objFn = Msd
		isLoss = true
		isDistFn = true
		// similarity loss
	case "Par":
		objFn = Par
		isLoss = true
		isDistFn = false
	case "Psis":
		objFn = Psis
		isLoss = true
		isDistFn = false
	case "Bers":
		objFn = Bers
		isLoss = true
		isDistFn = false
	case "Gar12":
		objFn = Gar12
		isLoss = true
		isDistFn = true
	case "Gar25":
		objFn = Gar25
		isLoss = true
		isDistFn = true
	case "Gar37":
		objFn = Gar37
		isLoss = true
		isDistFn = true
	case "Gar50":
		objFn = Gar50
		isLoss = true
		isDistFn = true
	case "Gar75":
		objFn = Gar75
		isLoss = true
		isDistFn = true
	case "Gar112":
		objFn = Gar112
		isLoss = true
		isDistFn = true
	case "Gar125":
		objFn = Gar125
		isLoss = true
		isDistFn = true
	case "Gar187":
		objFn = Gar187
		isLoss = true
		isDistFn = true
	case "Gar375":
		objFn = Gar375
		isLoss = true
		isDistFn = true
	case "Rgar12":
		objFn = Rgar12
		isLoss = true
		isDistFn = true
	case "Rgar25":
		objFn = Rgar25
		isLoss = true
		isDistFn = true
	case "Rgar37":
		objFn = Rgar37
		isLoss = true
		isDistFn = true
	case "Rgar50":
		objFn = Rgar50
		isLoss = true
		isDistFn = true
	case "Rgar75":
		objFn = Rgar75
		isLoss = true
		isDistFn = true
	case "Rgar112":
		objFn = Rgar112
		isLoss = true
		isDistFn = true
	case "Rgar125":
		objFn = Rgar125
		isLoss = true
		isDistFn = true
	case "Rgar187":
		objFn = Rgar187
		isLoss = true
		isDistFn = true
	case "Rgar250":
		objFn = Rgar250
		isLoss = true
		isDistFn = true
	case "Rgar375":
		objFn = Rgar375
		isLoss = true
		isDistFn = true

	}
	// print header
	fmt.Println("=======================================================")
	fmt.Println("Objective function: ", objf)
	fmt.Println("Matrix: ", matrix)
	fmt.Println("Iterations: ", nIter)
	fmt.Println("Heuristic for search: ", method)
	fmt.Println("Improvement function #", impro)
	fmt.Println("Window width for improvement: ", window)
	fmt.Println("Heuristic for improvement: ", improMethod)
	fmt.Println("Seed: ", seed)
	fmt.Println("=======================================================")
	fmt.Println()


	// start with the same conditions
	rand.Seed(seed)

	bestPerm, success := Seriate(a, objFn, isLoss, isDistFn, optMethod, optMethodForImpro, impro, window, nIter)

	fmt.Println("Best permutation found:")
	bestPerm.Print()
	if success {
		fmt.Println("Sorted matrix is Robinson.")
	} else {
		fmt.Println("Sorted matrix is NOT Robinson.")
	}
// Print out seriated matrix
	fmt.Println()
	fmt.Println("Seriated matrix:")
a.Permute(bestPerm, bestPerm)
a.Print()
// a.WriteCSV()
}
