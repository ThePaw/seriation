package main

// Program to evaluate the performance of different objective functions for seriation of a square symmetric similarity matrix.

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"time"
)

func main() {
	var (
		optMethod, optMethodForImpro OptMethod3
		objFn                        ObjFn
		isLoss, isDistFn             bool
		seed                         int64
	)

	// ============= Set parameters here: =======================

	nIter := 100
	method := "RobSA3"
	impro := 0
	window := 0
	improMethod := "RobSA3"
	seed = 33

	// =======================================================

	var (
		mtxp  = flag.String("m", "R50.csv", "similarity matrix to be seriated")
		objfp = flag.String("o", "Lsq", "objective function")
	)

	flag.Parse()

	matrix := *mtxp
	objf := *objfp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)

	nSamp := a.Rows()

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
	//	fmt.Println("Sorted matrix: ")
	//	a.Print()
	fmt.Println("Best permutations: ")

	t0 := time.Now()

	// start with the same conditions
	rand.Seed(seed)

	rhoH, rankH, pOH, rhoMean, rhoStDev, rProp, hitsProp := ObjFnPerformance(a, objFn, isLoss, isDistFn, optMethod, optMethodForImpro, impro, window, nIter)

	fmt.Println()
	fmt.Println("Rank matrix")
	rankH.Print()
	fmt.Println("Pair-order violations matrix")
	pOH.Print()

	fmt.Println("Rho histogram: ")
	for i := 0; i < rhoH.Len(); i++ {
		fmt.Printf("%d ", rhoH[i])
	}
	fmt.Println()
	// mean rank hits
	sumDiag := 0
	for i := 0; i < nSamp; i++ {
		sumDiag += rankH[i][i]
	}
	rankHitsProp := float64(sumDiag) / float64(nIter*nSamp)
	fmt.Println()
	dur := time.Since(t0)
	sec := dur.Seconds()
	fmt.Println("Summary:")
	fmt.Println("objF nSamp nIter rhoMean rhoStDev rProp hitsProp rankHitsProp seconds")
	fmt.Println(objf, nSamp, nIter, rhoMean, rhoStDev, rProp, hitsProp, rankHitsProp, sec)
}
