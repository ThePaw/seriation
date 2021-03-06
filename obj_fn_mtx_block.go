// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective (loss and gain) functions for block seriation of m x n data matrices.

import (
	//	"fmt"
	"github.com/ThePaw/go-fn/fn"
	"math"
)

// Gain functions.

// Me returns the Measure of Effectiveness  of a permuted data matrix(McCormick 1972). It is a gain function.
func Me(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	var x0, x1, x2, x3, x4 float64
	rows, cols := mtx.Dims()

	if !(rowPerm.Len() == rows && colPerm.Len() == cols) {
		panic("bad dimensions")
	}
	gain := 0.0
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			x0 = mtx[rowPerm[i]][colPerm[j]]
			if j-1 < 0 {
				x1 = 0
			} else {
				x1 = mtx[rowPerm[i]][colPerm[j-1]]
			}
			if j+1 > cols-1 {
				x2 = 0
			} else {
				x2 = mtx[rowPerm[i]][colPerm[j+1]]
			}
			if i-1 < 0 {
				x3 = 0
			} else {
				x3 = mtx[rowPerm[i-1]][colPerm[j]]
			}

			if i+1 > rows-1 {
				x4 = 0
			} else {

				x4 = mtx[rowPerm[i+1]][colPerm[j]]
			}
			gain += x0 * (x1 + x2 + x3 + x4)
		}
	}
	return gain / 2
}


// Loss functions.

// Mir computes energy E(p) of a permuted data matrix according to Miklos (2005:3400), Eqs. 2, 3.  It is a loss function.
func Mir(mtx Matrix64, rowPerm, colPerm IntVector, normalize bool) float64 {
	// normalize: Eq 3 of Miklos (2005:3400)

	var av float64
	rows, cols := mtx.Dims()
	if !(rowPerm.Len() == rows && colPerm.Len() == cols) {
		panic("bad dimensions")
	}
	if normalize {
		av = 0.0
		for i := 0; i < rows; i++ {
			for k := 0; k < cols; k++ {
				for l := 0; l < cols; l++ {
					if k < l {
						av += math.Abs(mtx[rowPerm[i]][colPerm[k]] - mtx[rowPerm[i]][colPerm[l]])
					}
				}
			}
		}

		for k := 0; k < cols; k++ {
			for i := 0; i < rows; i++ {
				for j := 0; j < rows; j++ {
					if i < j {
						av += math.Abs(mtx[rowPerm[i]][colPerm[k]] - mtx[rowPerm[j]][colPerm[k]])
					}
				}
			}
		}

		denom := float64(rows)*fn.BinomCoeff(int64(cols), 2) + float64(cols)*fn.BinomCoeff(int64(rows), 2)
		av /= denom
	}
	loss := 0.0

	//sum #1
	for i := 0; i < rows; i++ {
		for j := 0; j < cols-1; j++ {
			loss += math.Abs(mtx[rowPerm[i]][colPerm[j]] - mtx[rowPerm[i]][colPerm[j+1]])
		}
	}

	//sum #2
	for i := 0; i < rows-1; i++ {
		for j := 0; j < cols; j++ {
			loss += math.Abs(mtx[rowPerm[i]][colPerm[j]] - mtx[rowPerm[i+1]][colPerm[j]])
		}
	}

	//sum #3
	for j := 0; j < cols; j++ {
		loss += (math.Abs(mtx[rowPerm[0]][colPerm[j]]-mtx[rowPerm[1]][colPerm[j]]) + (math.Abs(mtx[rowPerm[rows-2]][colPerm[j]] - mtx[rowPerm[rows-1]][colPerm[j]])))
	}

	//sum #4
	for i := 0; i < rows; i++ {
		loss += (math.Abs(mtx[rowPerm[i]][colPerm[0]]-mtx[rowPerm[i]][colPerm[1]]) + (math.Abs(mtx[rowPerm[i]][colPerm[cols-2]] - mtx[rowPerm[i]][colPerm[cols-1]])))
	}
	if normalize {
		loss /= av
	}
	return loss
}

