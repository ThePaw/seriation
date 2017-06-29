package ser

import (
	"fmt"
)

func Seriate2(sim Matrix64, objFn ObjFn, isLoss, isDistFn bool, optMethod, optMethodForImpro OptMethod3, impro, window, tries, nIter int, pre, again bool) (bestPerm IntVector, success bool) {
	var (
		bestValue, cost float64
	)

	// init
	nSamp := sim.Rows()

	// alloc slices
	p := NewIntVector(nSamp) // known permutation
	p.Order()
	p.Perm()
	bestPerm = p.Clone()

	bestValue = +inf
	a := sim.Clone() // essential, because input matrix may be converted to distances!

	// if objFn is distance-based, convert similarities to distances
	if isDistFn {
		a.SimToDist()
	}

	if pre {
		// use Ham to pre-process the permutation
		cost = optMethod(a, p, Ham, true) // cost inverted inside if !isLoss
	}

	for it := 0; it < nIter; it++ {
		fmt.Println("%%% Iteration", it+1)

		// solve for best permutation
		cost = optMethod(a, p, objFn, isLoss) // cost inverted inside if !isLoss

		// update if improved
		if cost < bestValue {
			bestValue = cost
			bestPerm = p
			fmt.Println(bestValue)
		}

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
			// SegmentImpro + SwapOpt  + SegmentOpt
			cost = SegmentImpro(a, p, window, tries, objFn, isLoss)
			// update if improved
			if cost < bestValue {
				bestValue = cost
				bestPerm = p
				fmt.Println("SegmentImpro: ", bestValue)
			}

			cost = SwapOpt(a, p, objFn, isLoss)
			if cost < bestValue {
				bestValue = cost
				bestPerm = p
				fmt.Println("SwapOpt: ", bestValue)
			}
			/*
			   					cost=SegmentOpt(a, p, 6, objFn, isLoss)
			   		if cost < bestValue {
			   			bestValue = cost
			   			bestPerm = p
			   fmt.Println("SegmentOpt6: ",bestValue)
			   }
			*/

		case 7:
			// use mutator.go functions
		default:
			// no improvement
		}
		// update if improved
		if cost < bestValue {
			bestValue = cost
			bestPerm = p
			fmt.Println("%%", bestValue, "%%")
		}

		// reverse, if needed
		reverseIfNeeded(p)

		// is sorted similarity/distance matrix an (A)R-matrix?
		aa := a.Clone()
		aa.Permute(p, p)
		if isDistFn {
			aa.SimToDist()
		}
		if aa.IsR() {
			success = true
			break
		}
	} // end iterations

	// solve again ?
	if again && !success {

		// solve for the best permutation
		cost = optMethod(a, p, objFn, isLoss) // cost inverted inside if !isLoss

		// update if improved
		if cost < bestValue {
			bestValue = cost
			bestPerm = p
			fmt.Println(bestValue)
		}
	}

	return
}
