// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Similarity matrices, float.

import (
	"math"
)

// Kendall71S computes similarity matrix using Kendall's 1971 Theorem I.
func Kendall71S(dat Matrix64) Matrix64 {
	// Ref.:
	return Adj2Sim(dat, "rows")
}

// DegOverlap computes a "degree of overlap" between the rows of a data matrix.
func DegOverlap(dat Matrix64) Matrix64 {
	// Ref.: Paleo3 280: 469
	// Not tested yet. Looks like 1- ManlyOverlap().

	rows := dat.Rows()
	s := NewMatrix64(rows, rows)
	for i, row := range dat {
		for h, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			sum3 := 0.0

			for j, _ := range row { // columns of "dat"
				sum1 += dat[i][j] * dat[i][h]
				sum2 += dat[i][j] * dat[i][j]
				sum3 += dat[h][j] * dat[h][j]
			}
			s[i][h] = sum1 / (math.Sqrt(sum2) * math.Sqrt(sum3)) //////////// ??????????????

		}
	}
	return s
}
