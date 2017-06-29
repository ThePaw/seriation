// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective (loss and gain) functions for distance (dissimilarity) matrices.

import (
	"fmt"
	"math"
	"os"
)

func f(x, y float64) float64 {
	if x < y {
		return 1
	}
	if x > y {
		return -1
	}
	return 0
}

func g(x, y float64) float64 {
	if x > y {
		return 1
	}
	return 0
}

// Gain functions

// Wrug within row unweighted gradient (WRUG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 6, g_{1}(\Psi). It is a gain function.
func Wrug(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0
	for k := 0; k < n-2; k++ {
		for l := k + 1; l < n-1; l++ {
			for m := l + 1; m < n; m++ {
				x := dis[p[k]][p[m]]
				y := dis[p[k]][p[l]]
				c += sign(x - y)
			}
		}
	}
	return float64(c)
}

// Wrcug returns within row and column unweighted gradient (WRCUG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 7, g_{2}(\Psi)). It is a gain function.
func Wrcug(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0
	for k := 0; k < n-2; k++ {
		for l := k + 1; l < n-1; l++ {
			for m := l + 1; m < n; m++ {
				x := dis[p[k]][p[m]]
				y := dis[p[k]][p[l]]
				c += sign(x - y)
				y = dis[p[l]][p[m]]
				c += sign(x - y)
			}
		}
	}
	return float64(c)
}

// Wrwg returns within row weighted gradient (WRWG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 8, g_{3}(\Psi)). It is a gain function.
func Wrwg(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for k := 0; k < n-2; k++ {
		for l := k + 1; l < n-1; l++ {
			for m := l + 1; m < n; m++ {
				x := dis[p[k]][p[m]]
				y := dis[p[k]][p[l]]
				c += x - y
			}
		}
	}
	return c
}

// Wrcwg returns within row and column weighted gradient (WRCWG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 9, g_{4}(\Psi)).  (? approx. -StrengLoss2).  It is a gain function.
func Wrcwg(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for k := 0; k < n-2; k++ {
		for l := k + 1; l < n-1; l++ {
			for m := l + 1; m < n; m++ {
				x := dis[p[k]][p[m]]
				y := dis[p[k]][p[l]]
				z := dis[p[l]][p[m]]
				c += 2*x - y - z
			}
		}
	}
	return c
}

// H returns Szczotka's gain criterion of a permuted distance matrix (Szczotka 1972; Hubert and Schultz 1976; Brusco and Stahl 2000: 201, Eq. 5, Z_{5} ; Brusco et al. 2008: 507, Eq. 7, h(\psi)). It is a gain function.
func H(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
	fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for i := 0; i < n-1; i++ {
		for j := i + 1; j < n; j++ {
			d := math.Abs(float64(i - j))
			x := dis[p[i]][p[j]]
			c += d * x
		}
	}
	return c
}

// HNormGain returns gain of the permuted matrix according to Szczotka 1972; see Brusco et al. 2008: 507-508, Eq. 7.
// TO BE IMPLEMENTED

// Ine returns the inertia gain criterion of a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 5, Eq. 11). It is a gain function.
func Ine(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			sum += dis[p[i]][p[j]] * math.Abs(float64((i-j)*(i-j)))
		}
	}
	return sum
}

// Loss functions.

// Lsq returns the least squares loss criterion of a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 5, Eq. 12). It is a loss function.
func Lsq(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			incr := dis[p[i]][p[j]] - math.Abs(float64(i-j))
			incr *= incr
			sum += incr
		}
	}
	return sum
}

// Msd returns the Moore stress loss criterion of the matrix, here applied to a permuted distance matrix, so that m=n (Niermann 2005: 42, Eq. 1, 2; Hahsler et al. 2008: 6, Eq. 15). It is a loss function.
func Msd(dis Matrix64, p IntVector) float64 {
	return Ms(dis, p, p)
}

// Nsd returns the von Neumann stress loss criterion of the matrix, here applied to a permuted distance matrix, so that m=n (Niermann 2005: 42). It is a loss function.
func Nsd(dis Matrix64, p IntVector) float64 {
	return Ns(dis, p, p)
}

// Gar returns generalised anti-Robinson violation loss criterion of a permuted distance matrix (Chen 2002; Tien et al. 2008; Wu et al. 2010: 773, GAR(w)). It is a loss function.
func Gar(dis Matrix64, p IntVector, w int) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	sum := 0.0
	for j := 0; j < n; j++ {
		for k := 0; k < j; k++ {
			for i := j - w; i < k; i++ {
				if i >= 0 {
					dik := dis[p[i]][p[k]]
					dij := dis[p[i]][p[j]]
					sum += g(dik, dij)
				}
			}
		}
	}
	for j := 0; j < n; j++ {
		for k := 0; k < j; k++ {
			for i := j - w; i < k; i++ {
				if i >= 0 {
					dkj := dis[p[k]][p[j]]
					dij := dis[p[i]][p[j]]
					sum += g(dkj, dij)
				}
			}
		}
	}
	return sum
}

// Gar5 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 5. It is a loss function.
func Gar5(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 5)
}

// Gar10 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 10. It is a loss function.
func Gar10(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 10)
}

// Gar12 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 12. It is a loss function.
func Gar12(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 12)
}

// Gar15 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 15. It is a loss function.
func Gar15(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 15)
}

// Gar25 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 25. It is a loss function.
func Gar25(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 25)
}

// Gar37 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 37. It is a loss function.
func Gar37(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 37)
}

// Gar50 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 50. It is a loss function.
func Gar50(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 50)
}

// Gar75 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 75. It is a loss function.
func Gar75(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 75)
}

// Gar112 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 112. It is a loss function.
func Gar112(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 112)
}

// Gar125 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 125. It is a loss function.
func Gar125(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 125)
}

// Gar187 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 187. It is a loss function.
func Gar187(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 187)
}

// Gar250 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 250. It is a loss function.
func Gar250(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 250)
}

// Gar375 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 375. It is a loss function.
func Gar375(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 375)
}

// Rgar returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix(Chen 2002; Tien et al. 2008). It is a loss function.
func Rgar(dis Matrix64, p IntVector, w int) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	gar := Gar(dis, p, w)
	return gar / (float64(n*w*(w-1)) - 2*float64(w)*float64(w*w-1)/3)
}

// Rgar5 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 5. It is a loss function.
func Rgar5(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 5)
}

// Rgar10 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 10. It is a loss function.
func Rgar10(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 10)
}

// Rgar12 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 12. It is a loss function.
func Rgar12(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 12)
}

// Rgar15 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 15. It is a loss function.
func Rgar15(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 15)
}

// Rgar25 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 25. It is a loss function.
func Rgar25(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 25)
}

// Rgar37 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 37. It is a loss function.
func Rgar37(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 37)
}

// Rgar50 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 50. It is a loss function.
func Rgar50(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 50)
}

// Rgar75 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 75. It is a loss function.
func Rgar75(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 75)
}

// Rgar112 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 112. It is a loss function.
func Rgar112(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 112)
}

// Rgar125 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 125. It is a loss function.
func Rgar125(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 125)
}

// Rgar187 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 187. It is a loss function.
func Rgar187(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 187)
}

// Rgar250 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 250. It is a loss function.
func Rgar250(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 250)
}

// Rgar375 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 375. It is a loss function.
func Rgar375(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 375)
}

// Ham returns the loss criterion as the length of the shortest Hamiltonian path (open traveling salesman problem, oTSP) through a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 4 , Eq. 10; Chen 2002: 2, MS). It is a loss function.
func Ham(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	sum := 0.0
	for i := 0; i < n-1; i++ {
		sum += dis[p[i]][p[i+1]]
	}
	return sum
}

// Are returns anti-Robinson events violation loss criterion of a permuted distance matrix (Hahsler et al. 2008: 4, Eq. 7, 8; Streng 1991; Streng and Schönfelder 1978; Chen 2002: 2, AR(i)  ; Wu et al. 2010: 773, AR_{n}). It is a loss function.
func Are(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				x := dis[p[i]][p[k]]
				y := dis[p[i]][p[j]]
				c += g(x, y)
			}
		}
	}
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				x := dis[p[k]][p[j]]
				y := dis[p[i]][p[j]]
				c += g(x, y)
			}
		}
	}
	return c
}

// Ware returns weighted anti-Robinson events violation loss criterion of a permuted distance matrix (Hahsler et al. 2008; Streng 1991; Tien et al. 2008; Streng and Schönfelder 1978; Chen 2002:21, AR(s); Wu et al. 2010: 773, AR_{s}). It is a loss function.
func Ware(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				x := dis[p[i]][p[k]]
				y := dis[p[i]][p[j]]
				d := math.Abs(x - y)
				c += d * g(x, y)
			}
		}
	}
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				x := dis[p[k]][p[j]]
				y := dis[p[i]][p[j]]
				d := math.Abs(x - y)
				c += d * g(x, y)
			}
		}
	}
	return c
}

// Dware returns doubly weighted anti-Robinson violation loss criterion of a permuted distance matrix (Hahsler et al. 2008; Streng 1991; Tien et al. 2008; Chen 2002: 21, AR(w); Wu et al. 2010: 773, AR_{w}). It is a loss function.
func Dware(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		fmt.Fprintln(os.Stderr, "warning: distance matrix is not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	c := 0.0
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				jk := math.Abs(float64(j - k))
				x := dis[p[i]][p[k]]
				y := dis[p[i]][p[j]]
				d := math.Abs(x - y)
				c += jk * d * g(x, y)
			}
		}
	}
	for i := 0; i < n-2; i++ {
		for j := i + 2; j < n; j++ {
			for k := i + 1; k < j; k++ {
				ik := math.Abs(float64(i - k))
				x := dis[p[k]][p[j]]
				y := dis[p[i]][p[j]]
				d := math.Abs(x - y)
				c += ik * d * g(x, y)
			}
		}
	}
	return c
}
