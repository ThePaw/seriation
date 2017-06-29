package ser

import (
	//	"fmt"
	"math"
	"math/rand"
)

// TruncSampler just truncates the matrix of species abundances to integer values.
func TruncSampler(mtx Matrix64) IntMatrix {
	nSamp := mtx.Rows()
	nSpec := mtx.Cols()
	out := NewIntMatrix(nSamp, nSpec)
	for i := 0; i < nSamp; i++ {
		for j := 0; j < nSpec; j++ {
			x := math.Floor(mtx[i][j])
			if x < 0 {
				x = 0
			}
			out[i][j] = int(x)
		}
	}
	return out
}

func QAPserPerformance2(mtx Matrix64, nIter, trials, improLagMax, r int, exp, exp2, exp3 float64, di, fl int, seed int64) (rhoH IntVector, rankH, pOH IntMatrix, rhoMean, rhoStDev, qProp, hitsProp float64, empty int) {
	var (
		dist, flow, imtx IntMatrix
	)

	// init
	rand.Seed(seed)
	totalCost := 0
	hitsSum := 0.0
	rhoSum := 0.0
	qSum := 0.0
	empty = 0

	nSamp, nSpec := mtx.Dims()

	// alloc slices
	rhoH = NewIntVector(20)            // rho histogram
	pOH = NewIntMatrix(nSamp, nSamp)   // pair-order histogram
	rankH = NewIntMatrix(nSamp, nSamp) // ranks histogram
	//pKnown := NewIntVector(nSamp)      // known permutation
	outCC := NewIntMatrix(nSamp, nSpec) // seriated data matrix

	//pKnown.Order()
	//p := pKnown.Clone()

	// apply sampler
	imtx = TruncSampler(mtx)

	for it := 0; it < nIter; it++ {
		isEmpty := true
		for isEmpty {

			//check whether there are no empty samples
			//isEmpty = false
			for i := 0; i < nSamp; i++ {
				isEmpty = true
				for j := 0; j < nSpec; j++ {
					if imtx[i][j] > 0 {
						isEmpty = false
						break
					}
				}
				if isEmpty {
					empty++
					break
				}
			}
		}

		// prepare distance and flow matrices
		switch di {
		case 0:
			dist = distances(nSamp)
		case 1:
			dist = distances2(nSamp, exp3)
		case 2:
			dist = distances3(nSamp, exp3)
		case 3:
			dist = distances4(nSamp, exp3)
		default:
			dist = distances3(nSamp, exp3)
		}

		switch fl {
		case 0:
			flow = flows(imtx, exp)
		case 1:
			flow = flows2(imtx, exp, exp2)
		case 2:
			flow = flows3(imtx, exp, exp2)
		default:
			flow = flows3(imtx, exp, exp2)
		}

		// solve QAP
		cc, p := QAP_fant(dist, flow, trials, improLagMax, r)
		totalCost += cc

		// reverse, if needed
		rho := reverseIfNeeded2(p)

		// rank correlation
		rr := math.Abs(rho)
		rhoSum += rr

		// add pair-orders to pair-order histogram
		addToPairOrderHistogram(p, pOH)

		// add ranks to histogram
		addToRankHistogram(p, rankH)

		// add rho to histogram
		addToRhoHistogram(rr, rhoH)

		// update perfect hits
		if rr == 1 {
			hitsSum++
		}

		// is sorted Coenocline a Q-matrix?
		for i := 0; i < nSamp; i++ {
			for j := 0; j < nSpec; j++ {
				outCC[i][j] = imtx[p[i]][j]
			}
		}
		if outCC.IsQ() {
			qSum++
		}

	}

	// calc mean and st. deviation
	rhoMean = rhoSum / float64(nIter)
	rhoStDev /= float64(nIter - 1)
	rhoStDev = math.Sqrt(rhoStDev)

	// turn sums to proportions
	hitsProp = hitsSum / float64(nIter)
	qProp = qSum / float64(nIter)

	return
}
