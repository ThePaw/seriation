// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective (loss and gain) functions for m x n data matrices.

import (
	"math"
)

// Ms returns the Moore Stress criterion  of a permuted data matrix (Niermann 2005:42, Eq. 1, 2).
func Ms(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	r, c := mtx.Dims()
	if !(rowPerm.Len() == r && colPerm.Len() == c) {
		panic("bad dimensions")
	}
	stress := 0.0
	for i := 1; i <= r; i++ {
		for j := 1; j <= c; j++ {
			for l := imax(1, i-1); l <= imin(r, i+1); l++ {
				for m := imax(1, j-1); m <= imin(c, j+1); m++ {
					val := mtx[rowPerm[i-1]][colPerm[j-1]] - mtx[rowPerm[l-1]][colPerm[m-1]]
					val *= val
					stress += val
				}
			}
		}
	}
	return stress
}

// Ns returns the  VonNeumann Stress criterion  of a permuted data matrix (Niermann 2005:42).
func Ns(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	r, c := mtx.Dims()
	if !(rowPerm.Len() == r && colPerm.Len() == c) {
		panic("bad dimensions")
	}
	stress := 0.0
	for i := 1; i <= r; i++ {
		for j := 1; j <= c; j++ {
			for l := imax(1, i-1); l <= imin(r, i+1); l++ {
				for m := imax(1, j-1); m <= imin(c, j+1); m++ {
					if l == i || m == j {
						val := mtx[rowPerm[i-1]][colPerm[j-1]] - mtx[rowPerm[l-1]][colPerm[m-1]]
						val *= val
						stress += val
					}
				}
			}
		}
	}
	return stress
}

// Psi computes energy ψ(p) of a permuted data matrix according to Podani (1994); see  Miklos (2005), Eq. 4.
func Psi(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	rows, cols := mtx.Dims()
	if !(rowPerm.Len() == rows && colPerm.Len() == cols) {
		panic("bad dimensions")
	}
	loss := 0.0
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			x := mtx[rowPerm[i]][colPerm[j]]
			a := math.Abs(float64(cols*(i+1))/float64(rows) - float64(j+1))
			b := math.Abs(float64(rows*(j+1))/float64(cols) - float64(i+1))
			loss += x*a + b
		}
	}
	return loss
}

// Ber2 returns loss of a permuted data matrix according to Kostopoulos & Goulermas (in preparation) = Bertin Classification Criterion of Pilhofer 2012: 2509, Eq. 1.
func Ber2(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 1; i < n; i++ {
		for j := 0; j < m-1; j++ {
			tmp := float64(0)
			for k := 0; k <= i-1; k++ {
				for l := j + 1; l < m; l++ {
					tmp += mtx[rowPerm[k]][colPerm[l]]
				}
			}
			sum += tmp * mtx[rowPerm[i]][colPerm[j]]
		}
	}
	return sum
}

///////////////////// Untested functions

// Bg returns  B(A) gain a permuted data matrix of Pilhofer 2012: 2509, Eq. 1.
func Bg(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 1; i < n; i++ {
		for j := 0; j < m-1; j++ {
			tmp := float64(0)
			for k := 0; k <= i-1; k++ {
				for l := 0; l < j; l++ {
					tmp += mtx[rowPerm[k]][colPerm[l]]
				}
			}
			sum += tmp * mtx[rowPerm[i]][colPerm[j]]
		}
	}
	return sum
}

// Ber returns Bertin loss of a permuted data matrix according to Kostopoulos & Goulermas MATLAB code.
func Ber(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 1; i < n; i++ {
		for j := 0; j < m-1; j++ {
			for k := 0; k <= i-1; k++ {
				for l := j + 1; l < m; l++ {
					sum += mtx[rowPerm[k]][colPerm[l]] * mtx[rowPerm[i]][colPerm[j]]
				}
			}
		}
	}
	return sum
}

// X(A) of Pilhofer 2012: 2509
func xA(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			for l := 0; l < m; l++ {
				if j != l {
					sum += mtx[rowPerm[i]][colPerm[j]] * mtx[rowPerm[i]][colPerm[l]]
				}
			}
		}
	}
	return sum
}

// yA returns Y(A) of Pilhofer 2012: 2509.
func yA(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			for k := 0; k < n; k++ {
				if i != k {
					sum += mtx[rowPerm[i]][colPerm[j]] * mtx[rowPerm[i]][colPerm[j]]
				}
			}
		}
	}
	return sum
}

// biA returns BI(A) of Pilhofer 2012: 2509.
func biA(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	b1 := Bg(mtx, rowPerm, colPerm)
	b2 := Ber(mtx, rowPerm, colPerm)
	x := xA(mtx, rowPerm, colPerm)
	y := yA(mtx, rowPerm, colPerm)
	nm2 := float64(n * n * m * m)

	return (b1 + b2 + x) * (b1 + b2 + y) / nm2

}

// Bci returns the Bertin Classification Index (BCI) of a permuted data matrix of Pilhofer 2012: 2509, Eq. 2.
func Bci(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	b1 := Ber(mtx, rowPerm, colPerm)
	b2 := biA(mtx, rowPerm, colPerm)
	return b1 / b2
}

// Bwh returns the Weighted Bertin Classification Criterion (WBCC)  of a permuted data matrix of Pilhofer 2012: 2509  using Hamming distance.
func Bwh(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			for k := 0; k <= i; k++ {
				for l := 0; l <= j; l++ {
					w1 := math.Abs(float64(i - k))
					w2 := math.Abs(float64(j - l))
					w := w1 + w2
					sum += w * mtx[rowPerm[i]][colPerm[j]] * mtx[rowPerm[k]][colPerm[l]]

				}
			}
		}
	}
	return sum
}

// Bwe returns the Weighted Bertin Classification Criterion (WBCC)  of a permuted data matrix of Pilhofer 2012: 2509  using Euclidean distance.
func Bwe(mtx Matrix64, rowPerm, colPerm IntVector) float64 {
	n, m := mtx.Dims()
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			for k := 0; k <= i; k++ {
				for l := 0; l <= j; l++ {
					w1 := (i - k) * (i - k)
					w2 := (j - l) * (j - l)
					w := math.Sqrt(float64(w1 + w2))
					sum += w * mtx[rowPerm[i]][colPerm[j]] * mtx[rowPerm[k]][colPerm[l]]

				}
			}
		}
	}
	return sum
}
