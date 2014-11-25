// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective (loss and gain) functions for block seriation of similarity matrices.

import (
	//	"fmt"
	"code.google.com/p/go-fn/fn"
	"math"
)

// Gain functions.

// Mes returns the Measure of Effectiveness  of a permuted similarity matrix(McCormick 1972). It is a gain function.
func Mes(mtx Matrix64, p IntVector) float64 {
	var x0, x1, x2, x3, x4 float64
	rows, cols := mtx.Dims()

	if !(p.Len() == rows && p.Len() == cols) {
		panic("bad dimensions")
	}
	gain := 0.0
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			x0 = mtx[p[i]][p[j]]
			if j-1 < 0 {
				x1 = 0
			} else {
				x1 = mtx[p[i]][p[j-1]]
			}
			if j+1 > cols-1 {
				x2 = 0
			} else {
				x2 = mtx[p[i]][p[j+1]]
			}
			if i-1 < 0 {
				x3 = 0
			} else {
				x3 = mtx[p[i-1]][p[j]]
			}

			if i+1 > rows-1 {
				x4 = 0
			} else {

				x4 = mtx[p[i+1]][p[j]]
			}
			gain += x0 * (x1 + x2 + x3 + x4)
		}
	}
	return gain / 2
}


// Loss functions.

// Mirs computes energy E(p) of a permuted similarity matrix according to Miklos (2005:3400), Eqs. 2, 3.  It is a loss function.
func Mirs(mtx Matrix64, p IntVector, normalize bool) float64 {
	// normalize: Eq 3 of Miklos (2005:3400)

	var av float64
	rows, cols := mtx.Dims()
	if !(p.Len() == rows && p.Len() == cols) {
		panic("bad dimensions")
	}
	if normalize {
		av = 0.0
		for i := 0; i < rows; i++ {
			for k := 0; k < cols; k++ {
				for l := 0; l < cols; l++ {
					if k < l {
						av += math.Abs(mtx[p[i]][p[k]] - mtx[p[i]][p[l]])
					}
				}
			}
		}

		for k := 0; k < cols; k++ {
			for i := 0; i < rows; i++ {
				for j := 0; j < rows; j++ {
					if i < j {
						av += math.Abs(mtx[p[i]][p[k]] - mtx[p[j]][p[k]])
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
			loss += math.Abs(mtx[p[i]][p[j]] - mtx[p[i]][p[j+1]])
		}
	}

	//sum #2
	for i := 0; i < rows-1; i++ {
		for j := 0; j < cols; j++ {
			loss += math.Abs(mtx[p[i]][p[j]] - mtx[p[i+1]][p[j]])
		}
	}

	//sum #3
	for j := 0; j < cols; j++ {
		loss += (math.Abs(mtx[p[0]][p[j]]-mtx[p[1]][p[j]]) + (math.Abs(mtx[p[rows-2]][p[j]] - mtx[p[rows-1]][p[j]])))
	}

	//sum #4
	for i := 0; i < rows; i++ {
		loss += (math.Abs(mtx[p[i]][p[0]]-mtx[p[i]][p[1]]) + (math.Abs(mtx[p[i]][p[cols-2]] - mtx[p[i]][p[cols-1]])))
	}
	if normalize {
		loss /= av
	}
	return loss
}

