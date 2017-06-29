package main

// Program to seriate a square symmetric similarity matrix.

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
	)
	mtxp := flag.String("f", "sim.csv", "similarity matrix to be seriated")
	objfp := flag.String("o", "H", "objective function")
	iterp := flag.Int("i", 1, "number of iterations")
	methodp := flag.String("m", "RobSA3", "method")
	improp := flag.Int("p", 0, "improvement type")
	windowp := flag.Int("w", 0, "window width for improvement")
	triesp := flag.Int("t", 0, "number of tries for improvement")
	improMethodp := flag.String("M", "RobSA3", "method   for improvement")
	seedp := flag.Int64("s", 33, "seed")
	isDistP := flag.Bool("d", false, "is input dissimilarity matrix?")
	verboseP := flag.Bool("v", false, "verbose output?")

	flag.Parse()

	matrix := *mtxp
	objf := *objfp
	nIter := *iterp
	method := *methodp
	impro := *improp
	window := *windowp
	tries := *triesp
	improMethod := *improMethodp
	seed := *seedp
	isDist := *isDistP
	verbose := *verboseP

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file)
	if isDist {
		a.SimToDist()
	}

	nSamp := a.Rows()

	b := a.Clone()

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

	objFn, isLoss, isDistFn = SelectObjFn(objf)

	if verbose {
		// print header
		fmt.Println("=======================================================")
		fmt.Println("Objective function: ", objf)
		fmt.Println("Matrix: ", matrix)
		fmt.Println("Iterations: ", nIter)
		fmt.Println("Heuristic for search: ", method)
		fmt.Println("Improvement function #", impro)
		fmt.Println("Window width for improvement: ", window)
		fmt.Println("Without SegmentOpt")
		fmt.Println("Tries for improvement: ", tries)
		fmt.Println("Heuristic for improvement: ", improMethod)
		fmt.Println("Seed: ", seed)
		fmt.Println("=======================================================")
		fmt.Println()
	}

	// start with the same conditions
	rand.Seed(seed)

	t0 := time.Now()
	bestPerm, success := Seriate(a, objFn, isLoss, isDistFn, optMethod, optMethodForImpro, impro, window, tries, nIter)
	dur := time.Since(t0)
	sec := dur.Seconds()

	if verbose {
		if success {
			fmt.Println("Sorted matrix IS Robinson.")
		} else {
			fmt.Println("Sorted matrix is NOT Robinson.")
		}
		fmt.Println("seconds: ", sec)

		b.SimToDist()
		hValue := H(b, bestPerm)
		fmt.Println("hValue: ", -hValue)
		fmt.Println("Best permutation that was found:")
	}
	bestPerm.WriteCSV()
}
