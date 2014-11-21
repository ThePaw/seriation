package ser

import (
)

func Seriate(sim Matrix64, objFn ObjFn, isLoss, isDistFn bool, optMethod, optMethodForImpro OptMethod3, impro, window, nIter int) (bestPerm IntVector, success bool) {

	// init
	nSamp := sim.Rows()

	// alloc slices
	pKnown := NewIntVector(nSamp) // known permutation
	pKnown.Order()
	bestPerm = pKnown.Clone()
	bestValue := +inf

	for it := 0; it < nIter; it++ {
		if success {
			break
		}
		p := pKnown.Clone()
		a := sim.Clone() // essential, because input matrix may be converted to distances!
		p.Perm()
		//		a.ForceTo01()

		// if objFn is distance-based, convert similarities to distances
		if isDistFn {
			a.SimToDist()
		}

		// solve for best permutation
		value := optMethod(a, p, objFn, isLoss)

		// update if improved
		if value < bestValue {
			bestValue = value
			bestPerm = p
		}

		/*		// try to improve the solution
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

		if value < bestValue {
				bestValue = value
				bestPerm = p
		}

		*/
		// reverse, if needed
		reverseIfNeeded(p)

		// is sorted similarity/distance matrix (A)R-matrix?
		aa := a.Clone()
		aa.Permute(p, p)
		if isDistFn {
			aa.SimToDist()
		}
		if aa.IsR() {
			success = true
		}
	}
	return
}
