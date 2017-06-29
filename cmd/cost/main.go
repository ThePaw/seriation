package main

// Program to calculate the "cost" of a (permuted) square symmetric similarity matrix.

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"fmt"
	"os"
)

func main() {
	var (
		p                IntVector
		objFn            ObjFn
		isLoss, isDistFn bool
		cost             float64
	)
	mtxp := flag.String("f", "sim.csv", "similarity matrix")
	pp := flag.String("p", "none", "permutation vector")
	objfp := flag.String("o", "H", "objective function")
	forcep := flag.Bool("F", false, "force to 0 - 1?")
	verbosep := flag.Bool("v", false, "verbose output?")

	flag.Parse()

	matrix := *mtxp
	pvector := *pp
	objf := *objfp
	force := *forcep
	verbose := *verbosep

	file1, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}
	a := ReadCsvMatrix64(file1)

	if pvector != "none" {
		file2, err := os.Open(pvector) // For read access.
		if err != nil {
			panic("file does not exist")
		}
		p = ReadCsvIntVector(file2)
	}

	objFn, isLoss, isDistFn = SelectObjFn(objf)

	if verbose {
		// print header
		fmt.Println("=======================================================")
		fmt.Println("Objective function: ", objf)
		fmt.Println("Matrix: ", matrix)
		if pvector != "none" {
			fmt.Println("Permutation vector: ", pvector)
		}
		fmt.Println("=======================================================")
		fmt.Println()
	}

	// if objFn is distance-based, convert similarities to distances
	if isDistFn {
		a.SimToDist()
	}

	if force {
		a.ForceTo01()
	}

	if pvector != "none" {

		cost = Cost(a, p, objFn, isLoss, isDistFn)
	} else {
		cost = CostRaw(a, objFn, isLoss, isDistFn)
	}

	if verbose {
		p.Print()
	}

	fmt.Println(objf, cost)
	/*	p.Invert()
		p.Print()
		cost = Cost(a, p, objFn, isLoss, isDistFn)
		fmt.Println(objf, cost)
		p.Invert()
		p.Print()
		cost = Cost(a, p, objFn, isLoss, isDistFn)
		fmt.Println(objf, cost)
	*/
	if verbose {
		fmt.Println("=======================================================")
	}

}
