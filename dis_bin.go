// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Dissimilarity matrices, binary.
// Thanks to (R vegan).
// Not tested yet.
// a = J, b = A-J, c = B-J, d = P-A-B+J
// thus:
// J=a, A=b+a, B=c+a, P=d+b+a+c+a-a = a+b+c+d

import (
	"math"
)

func ABJPN(dat Matrix64, i, j int) (A, B, J, P, N float64) {
	a, b, c, d, _, _, _ := abcd(dat, i, j, false)
	J = a
	A = b + a
	B = c + a
	P = a + b + c + d
	N = float64(dat.Cols())
	// M=dat.Cols() excl.double 0
	return
}

func EuclideanBinD(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat { // rows of "dat"
		for j, _ := range dat { // rows of "dat"
			A, B, J, _, _ := ABJPN(dat, i, j)
			s[i][j] = math.Sqrt(A + B - 2*J)
		}
	}
	return
}
