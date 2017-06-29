// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Similarity matrices of compositional data (proportions, frequencies).

import (
	"math"
)

// PercentageS computes a "Percentage similarity" between the rows of a data matrix.
func PercentageS(dat Matrix64) (sim Matrix64) {
	// Ref.: http://ordination.okstate.edu/distsim.htm
	// Gauch (1982).
	// Not tested yet. It is asymmetric, should it be??

	rows := dat.Rows()
	sim = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += math.Min(p, q)
				sum2 += p + q
			}
			sim[i][j] = 200*(sum1/sum2)
			if sim[i][j] == math.NaN() {
				sim[i][j] = 0.0
			}
		}
	}
	return
}

