package main

// Program to calculate the "cost" of a square symmetric similarity matrix.

import (
	. "github.com/ThePaw/seriation"
	"flag"
	"fmt"
	"os"
)

func main() {
	mtxp := flag.String("f", "sim.csv", "similarity matrix")

	flag.Parse()

	matrix := *mtxp

	file1, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	a := ReadCsvMatrix64(file1)

	if a.IsR() {
		fmt.Println("Matrix is Robinsonian.")
	} else {
		a.SimToDist()
		if a.IsR() {
			fmt.Println("Matrix is anti-Robinsonian.")

		} else {
			fmt.Println("Matrix is NOT (anti)Robinsonian.")

		}

	}
}
