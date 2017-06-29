package ser

import (
	"fmt"
)

// Seriate does seriation of a similarity matrix.
func Seriate(sim Matrix64, objFn ObjFn, isLoss, isDistFn bool, optMethod, optMethodForImpro OptMethod3, impro, window, tries, nIter int) (bestPerm IntVector, success bool) {
	var (
		itBestCost, bestCost, cost float64
	)

	// init
	nSamp := sim.Rows()

	// alloc slices
	p := NewIntVector(nSamp) // random permutation
	p.Perm()
	itBestPerm := p.Clone()
	bestPerm = p.Clone()

	a := sim.Clone() // essential, because input matrix may be converted to distances!
	//		a.ForceTo01()

	// if objFn is distance-based, convert similarities to distances
	if isDistFn {
		a.SimToDist()
	}

	// iterate
	bestCost = +inf
	bestIt := 0
	for it := 0; it < nIter; it++ {
		itBestCost = bestCost
		p.Perm()

		// solve for best permutation
		cost = optMethod(a, p, objFn, isLoss) // cost inverted inside if !isLoss

		// update if improved
		if cost < itBestCost {
			itBestCost = cost
			itBestPerm.CopyFrom(p)
			fmt.Println("Iteration ", it+1, "SA:", itBestCost)
			itBestPerm.Print()
			//hValue := H(a, itBestPerm)
			//fmt.Println("hValue(SA): ", -hValue)

		}

		// try to improve the solution
		p.CopyFrom(itBestPerm)

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
			// SegmentImpro
			cost = SegmentImpro(a, p, window, tries, objFn, isLoss)
			// update if improved
			if cost < itBestCost {
				itBestCost = cost
				itBestPerm.CopyFrom(p)
				fmt.Println("Iteration ", it+1, "SegmentImpro: ", itBestCost)
				itBestPerm.Print()
				//hValue := H(a, itBestPerm)
				// fmt.Println("hValue(SegmentImpro): ", -hValue)
			}
			// SwapOpt
			p.CopyFrom(itBestPerm)

			cost = SwapOpt(a, p, objFn, isLoss)
			// update if improved
			if cost < itBestCost {
				itBestCost = cost
				itBestPerm.CopyFrom(p)
				fmt.Println("Iteration ", it+1, "SwapOpt: ", itBestCost)
				itBestPerm.Print()
				//hValue := H(a, itBestPerm)
				// fmt.Println("hValue(SwapOpt): ", -hValue)
			}
		default:
			// no improvement
		}

		// reverse, if needed
		reverseIfNeeded(itBestPerm)

		if itBestCost < bestCost {
			bestIt = it
			bestPerm.CopyFrom(itBestPerm)
			bestCost = itBestCost

			//			bestPerm.Print()

		}

		if a.IsR() {
			success = true
			break
		}

	} // end of iterations
	hValue := H(a, bestPerm)
	fmt.Println("Last improvement in iteration ", bestIt+1, "cost: ", bestCost, "hValue: ", -hValue)
	return
}

// Cost calculates objective function value for a similarity matrix given its permutation vector.
func Cost(sim Matrix64, p IntVector, objFn ObjFn, isLoss, isDistFn bool) (cost float64) {
	a := sim.Clone()
	if isDistFn {
		// if objFn is distance-based, convert similarities to distances

		a.SimToDist()
	}

	cost = objFn(a, p)

	if !isLoss {
		cost = -cost
	}

	return

}

// CostRaw calculates objective function value for a similarity matrix.
func CostRaw(sim Matrix64, objFn ObjFn, isLoss, isDistFn bool) (cost float64) {
	a := sim.Clone()
	nSamp := sim.Rows()
	p := NewIntVector(nSamp)
	p.Order()

	if isDistFn {
		// if objFn is distance-based, convert similarities to distances

		a.SimToDist()
	}

	cost = objFn(a, p)

	if !isLoss {
		cost = -cost
	}

	return

}

