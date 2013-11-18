package ser

import (
	"math"
	"os"
	"strings"
	"strconv"
	"fmt"
)

func ObjFnPerformanceMultiMtx(path string, matDim int, objFn ObjFn, isLoss, isDistFn bool, optMethod, optMethodForImpro OptMethod3, impro, window, nIter int) (rhoH IntVector, rankH, pOH IntMatrix, rhoMean, rhoStDev, rProp, hitsProp float64) {

	// init
	hitsSum := 0.0
	rhoSum := 0.0
	rhoM := 0.0
	rhoStDev = 0.0
	rSum := 0.0
	nSamp := matDim

	// alloc slices
	rhoH = NewIntVector(20)            // rho histogram
	pOH = NewIntMatrix(nSamp, nSamp)   // pair-order histogram
	rankH = NewIntMatrix(nSamp, nSamp) // ranks histogram
	aa := NewMatrix64(nSamp, nSamp)    // sorted similarity/distance matrix
	pKnown := NewIntVector(nSamp)      // known permutation
	pKnown.Order()
	p := pKnown.Clone()

	for it := 0; it < nIter; it++ {
		//read the matrix
		dir := strconv.Itoa(matDim)
		filename := strconv.Itoa(it)
		s := []string{path, dir, filename}
		foo := strings.Join(s, "/")

		file, err := os.Open(foo) // For read access.
		if err != nil {
			panic("file does not exist")
		}
		dist := ReadCsvMatrix64(file)
		nSamp := dist.Rows()
		if nSamp != matDim || !dist.IsSymmetric() {
			panic("bad matrix")
		}

		a := dist.Clone() // essential, because input matrix may be converted to similarity!

		p.Perm()
		//		a.ForceTo01()

			// if objFn is similarity-based, convert matrix to distances
			if !isDistFn {
				a.SimToDist()

				fmt.Println("CONVERTED")
			}

a.WriteCSV3()
		// solve for best permutation
		optMethod(a, p, objFn, isLoss)

		// try to improve the solution
		switch impro {
		case 0: // no improvement
		case 1:
			SegmentOpt(a, p, window, objFn, isLoss)
		case 2:
			SubMatOpt(a, p, window, objFn, isLoss, optMethodForImpro)
		case 3:
			SwapOpt(a, p, objFn, isLoss)
		case 4:
			RobSA3(a, p, objFn, isLoss)
		case 5:
			RobFA3(a, p, objFn, isLoss)
		case 6:
			// SegmentImpro + SwapOpt
			SegmentImpro(a, p, window, objFn, isLoss)
			SwapOpt(a, p, objFn, isLoss)
		default:
			// no improvement
		}

		// reverse, if needed
		rho := reverseIfNeeded2(p)

		// rank correlation
		rr := math.Abs(rho)

		// rho sample mean and unbiased (Bessel correction) variance estimates
		rhoSum += rr
		rhoDelta := rr - rhoM
		rhoM += rhoDelta / float64(it+1)
		rhoStDev += rhoDelta * (rr - rhoM)

		// add pair-orders to pair-order histogram
		addToPairOrderHistogram(p, pOH)

		// add ranks to histogram
		addToRankHistogram(p, rankH)

		// add rho to histogram
		addToRhoHistogram(rr, rhoH)

		// update perfect hits
		if p.Equals(pKnown) {
			hitsSum++
		}

		// is sorted similarity/distance matrix (A)R-matrix?
		for i := 0; i < nSamp; i++ {
			for j := 0; j < nSamp; j++ {
				aa[i][j] = a[p[i]][p[j]]
			}
		}

		if isDistFn {
			if aa.IsAR() {
				rSum++
			}
		} else {
			if aa.IsR() {
				rSum++
			}
		}
		p.Print()

	}

	// calc mean and st. deviation
	rhoMean = rhoSum / float64(nIter)

	rhoStDev /= float64(nIter - 1)
	rhoStDev = math.Sqrt(rhoStDev)

	// calc proportions
	hitsProp = hitsSum / float64(nIter)
	rProp = rSum / float64(nIter)
	return
}


func ObjFnPerformanceMultiMtxPerm(path string, matDim int, p IntVector, objFn ObjFn, isLoss, isDistFn bool, optMethod, optMethodForImpro OptMethod3, impro, window, nIter int) (rhoH IntVector, rankH, pOH IntMatrix, rhoMean, rhoStDev, rProp, hitsProp float64) {

	// init
	hitsSum := 0.0
	rhoSum := 0.0
	rhoM := 0.0
	rhoStDev = 0.0
	rSum := 0.0
	nSamp := matDim

	// alloc slices
	rhoH = NewIntVector(20)            // rho histogram
	pOH = NewIntMatrix(nSamp, nSamp)   // pair-order histogram
	rankH = NewIntMatrix(nSamp, nSamp) // ranks histogram
	aa := NewMatrix64(nSamp, nSamp)    // sorted similarity/distance matrix
	pKnown := NewIntVector(nSamp)      // known permutation
	pKnown.Order()

	for it := 0; it < nIter; it++ {
		//read the matrix
		dir := strconv.Itoa(matDim)
		filename := strconv.Itoa(it)
		s := []string{path, dir, filename}
		foo := strings.Join(s, "/")

		file, err := os.Open(foo) // For read access.
		if err != nil {
			panic("file does not exist")
		}
		dist := ReadCsvMatrix64(file)
		nSamp := dist.Rows()
		if nSamp != matDim {
			panic("bad matrix")
		}

		a := dist.Clone() // essential, because input matrix may be converted to distances!
		p.Perm()
		//		a.ForceTo01()

			// if objFn is similarity-based, convert matrix to distances
			if !isDistFn {
				a.SimToDist()
			}


		// solve for best permutation
		optMethod(a, p, objFn, isLoss)

		// try to improve the solution
		switch impro {
		case 0: // no improvement
		case 1:
			SegmentOpt(a, p, window, objFn, isLoss)
		case 2:
			SubMatOpt(a, p, window, objFn, isLoss, optMethodForImpro)
		case 3:
			SwapOpt(a, p, objFn, isLoss)
		case 4:
			RobSA3(a, p, objFn, isLoss)
		case 5:
			RobFA3(a, p, objFn, isLoss)
		case 6:
			// SegmentImpro + SwapOpt
			SegmentImpro(a, p, window, objFn, isLoss)
			SwapOpt(a, p, objFn, isLoss)
		default:
			// no improvement
		}

		// reverse, if needed
		rho := reverseIfNeeded2(p)

		// rank correlation
		rr := math.Abs(rho)

		// rho sample mean and unbiased (Bessel correction) variance estimates
		rhoSum += rr
		rhoDelta := rr - rhoM
		rhoM += rhoDelta / float64(it+1)
		rhoStDev += rhoDelta * (rr - rhoM)

		// add pair-orders to pair-order histogram
		addToPairOrderHistogram(p, pOH)

		// add ranks to histogram
		addToRankHistogram(p, rankH)

		// add rho to histogram
		addToRhoHistogram(rr, rhoH)

		// update perfect hits
		if p.Equals(pKnown) {
			hitsSum++
		}

		// is sorted similarity/distance matrix (A)R-matrix?
		for i := 0; i < nSamp; i++ {
			for j := 0; j < nSamp; j++ {
				aa[i][j] = a[p[i]][p[j]]
			}
		}

		if isDistFn {
			if aa.IsAR() {
				rSum++
			}
		} else {
			if aa.IsR() {
				rSum++
			}
		}
		p.Print()

	}

	// calc mean and st. deviation
	rhoMean = rhoSum / float64(nIter)

	rhoStDev /= float64(nIter - 1)
	rhoStDev = math.Sqrt(rhoStDev)

	// calc proportions
	hitsProp = hitsSum / float64(nIter)
	rProp = rSum / float64(nIter)
	return
}
