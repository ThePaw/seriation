package main

// Program to seriate data matrix using QAP.

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
		mtxp = flag.String("m", "dat/cc20x100sp.csv", "data matrix to be seriated")
		seed int64
	)

	// ============= Set parameters here: =======================
	seed = 33
	exp := 5.34
	exp2 := 1.34
	exp3 := 3.5
	di := 1
	fl := 2
	r := 5
	trials := 1
	improLagMax := 500
	nIter := 1
	// =======================================================

	flag.Parse()

	matrix := *mtxp

	file, err := os.Open(matrix) // For read access.
	if err != nil {
		panic("file does not exist")
	}

	mtx := ReadCsvMatrix64(file)
	nSamp := mtx.Rows()
	t0 := time.Now()

	// start with the same conditions
	rand.Seed(seed)

	rhoH, rankH, pOH, rhoMean, rhoStDev, qProp, hitsProp, empty := QAPserPerformance2(mtx, nIter, trials, improLagMax, r, exp, exp2, exp3, di, fl, seed)
	//rhoH, rankH, pOH, rhoMean, rhoStDev, rProp, hitsProp
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
	fmt.Println("nSamp, nIter, rhoMean, rhoStDev, qProp, hitsProp, rankHitsProp, empty, seconds")
	fmt.Println(nSamp, nIter, rhoMean, rhoStDev, qProp, hitsProp, rankHitsProp, empty, sec)
}
