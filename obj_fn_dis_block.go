// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective functions for block seriation of permuted distance matrices.
// Gain functions.

// Moed returns the Measure of Effectiveness gain of a permuted distance matrix (McCormick 1972) according to Kostopoulos & Goulermas (2014, in preparation), for block seriation.
func Moed(dis Matrix64, p IntVector) float64 {
	var a, b, c, d, e float64
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if j > n-2 {
				a = 0
			} else {
				a = dis[p[i]][p[j+1]]
			}
			if j == 0 {
				b = 0
			} else {
				b = dis[p[i]][p[j-1]]
			}

			if i > n-2 {
				c = 0
			} else {
				c = dis[p[i+1]][p[j]]
			}
			if i == 0 {
				d = 0
			} else {

				d = dis[p[i-1]][p[j]]
			}
			e = dis[p[i]][p[j]]
			sum += e * (a + b + c + d)
		}
	}
	return sum / 2
}

// Moed2 returns the Measure of Effectiveness gain of a permuted distance matrix (McCormick 1972), for block seriation.
func Moed2(a Matrix64, p IntVector) float64 {
	var x0, x1, x2, x3, x4 float64
	rows := a.Rows()

	if !(p.Len() == rows) {
		panic("bad permutation vector length")
	}
	gain := 0.0
	for i := 0; i < rows; i++ {
		for j := 0; j < rows; j++ {
			x0 = a[p[i]][p[j]]
			if j-1 < 0 {
				x1 = 0
			} else {
				x1 = a[p[i]][p[j-1]]
			}
			if j+1 > rows-1 {
				x2 = 0
			} else {
				x2 = a[p[i]][p[j+1]]
			}
			if i-1 < 0 {
				x3 = 0
			} else {
				x3 = a[p[i-1]][p[j]]
			}

			if i+1 > rows-1 {
				x4 = 0
			} else {

				x4 = a[p[i+1]][p[j]]
			}
			gain += x0 * (x1 + x2 + x3 + x4)
		}
	}
	return gain / 2
}
