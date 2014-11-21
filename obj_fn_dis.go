// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective (loss and gain) functions for distance (dissimilarity) matrices.

import (
	"math"
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

// Wrug within row unweighted gradient (WRUG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 6, g_{1}(\Psi).
func Wrug(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Wrcug returns within row and column unweighted gradient (WRCUG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 7, g_{2}(\Psi)).
func Wrcug(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Wrwg returns within row weighted gradient (WRWG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 8, g_{3}(\Psi)).
func Wrwg(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Wrcwg returns within row and column weighted gradient (WRCWG) gain of a permuted distance matrix (Hubert et al. 2001, Chapter 4; Brusco 2002: 50, Eq. 9, g_{4}(\Psi)).  (? approx. -StrengLoss2)
func Wrcwg(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// H returns Szczotka's gain criterion of a permuted distance matrix (Szczotka 1972; Hubert and Schultz 1976; Brusco and Stahl 2000: 201, Eq. 5, Z_{5} ; Brusco et al. 2008: 507, Eq. 7, h(\psi)).
func H(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
	}
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

// Ine returns the inertia gain criterion of a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 5, Eq. 11).
func Ine(dis Matrix64, p IntVector) float64 {
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
			sum += dis[p[i]][p[j]] * math.Abs(float64((i-j)*(i-j)))
		}
	}
	return sum
}

// Lsq returns the least squares loss criterion of a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 5, Eq. 12).
func Lsq(dis Matrix64, p IntVector) float64 {
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
			incr := dis[p[i]][p[j]] - math.Abs(float64(i-j))
			incr *= incr
			sum += incr
		}
	}
	return sum
}

// Msd returns the Moore stress loss criterion of the matrix, here applied to a permuted distance matrix, so that m=n (Niermann 2005: 42, Eq. 1, 2; Hahsler et al. 2008: 6, Eq. 15).
func Msd(dis Matrix64, p IntVector) float64 {
	return MooreStressLoss(dis, p, p)
}

// Nsd returns the von Neumann stress loss criterion of the matrix, here applied to a permuted distance matrix, so that m=n (Niermann 2005: 42).
func Nsd(dis Matrix64, p IntVector) float64 {
	return VonNeumannStressLoss(dis, p, p)
}

// Gar returns generalised anti-Robinson violation loss criterion of a permuted distance matrix (Chen 2002; Tien et al. 2008; Wu et al. 2010: 773, GAR(w)).
func Gar(dis Matrix64, p IntVector, w int) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Gar5 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 5.
func Gar5(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 5)
}

// Gar10 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 10.
func Gar10(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 10)
}

// Gar12 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 12.
func Gar12(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 12)
}

// Gar15 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 15.
func Gar15(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 15)
}

// Gar25 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 25.
func Gar25(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 25)
}

// Gar37 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 37.
func Gar37(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 37)
}

// Gar50 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 50.
func Gar50(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 50)
}

// Gar75 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 75.
func Gar75(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 75)
}

// Gar112 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 112.
func Gar112(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 112)
}

// Gar125 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 125.
func Gar125(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 125)
}

// Gar187 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 187.
func Gar187(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 187)
}

// Gar250 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 250.
func Gar250(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 250)
}

// Gar375 returns generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 375.
func Gar375(dis Matrix64, p IntVector) float64 {
	return Gar(dis, p, 375)
}

// Rgar returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix(Chen 2002; Tien et al. 2008).
func Rgar(dis Matrix64, p IntVector, w int) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
	}
	n := p.Len()
	if dis.Rows() != n {
		panic("bad permutation vector length")
	}

	gar := Gar(dis, p, w)
	return gar / (float64(n*w*(w-1)) - 2*float64(w)*float64(w*w-1)/3)
}

// Rgar5 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 5.
func Rgar5(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 5)
}

// Rgar10 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 10.
func Rgar10(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 10)
}

// Rgar12 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 12.
func Rgar12(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 12)
}

// Rgar15 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 15.
func Rgar15(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 15)
}

// Rgar25 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 25.
func Rgar25(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 25)
}

// Rgar37 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 37.
func Rgar37(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 37)
}

// Rgar50 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 50.
func Rgar50(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 50)
}

// Rgar75 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 75.
func Rgar75(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 75)
}

// Rgar112 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 112.
func Rgar112(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 112)
}

// Rgar125 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 125.
func Rgar125(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 125)
}

// Rgar187 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 187.
func Rgar187(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 187)
}

// Rgar250 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 250.
func Rgar250(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 250)
}

// Rgar375 returns relative generalised anti-Robinson violation loss criterion of a permuted distance matrix with window width = 375.
func Rgar375(dis Matrix64, p IntVector) float64 {
	return Rgar(dis, p, 375)
}

// Ham returns the loss criterion as the length of the shortest Hamiltonian path (open traveling salesman problem, oTSP) through a permuted distance matrix (Caraux and Pinloche 2005; Hahsler et al. 2008: 4 , Eq. 10; Chen 2002: 2, MS).
func Ham(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// parabolaFit returns coefficients of the polynomial c1 + c2*x + c3*x*x fitted to the data vector, where abscissa is 0, 1, 2, ... , n-1
func parabolaFit(v Vector64) (c1, c2, c3 float64, err bool) {
	eps := 1e-6
	n := float64(v.Len())

	sumX := 0.0
	sumY := 0.0
	sumXY := 0.0
	sumX2 := 0.0
	sumX2Y := 0.0
	sumX3 := 0.0
	sumX4 := 0.0
	for i, y := range v {
		x := float64(i)
		sumX += x
		sumY += y
		sumXY += x * y
		sumX2 += x * x
		sumX2Y += x * x * y
		sumX3 += x * x * x
		sumX4 += x * x * x * x
	}

	// The normal equations
	// sumY = c1*n + c2*sumX + c3*sumX2
	// sumXY = c1*sumX + c2*sumX2  + c3*sumX3
	// sumX2Y = c1*sumX2 + c2*sumX3  + c3*sumX4

	a := n
	b := sumX
	c := sumX2
	d := sumX
	e := sumX2
	f := sumX3
	g := sumX2
	h := sumX3
	k := sumX4

	b0 := sumY
	b1 := sumXY
	b2 := sumX2Y

	// determinant of a 3x3 matrix can be computed by applying the rule of Sarrus as follows:
	det := a*(e*k-f*h) - b*(k*d-f*g) + c*(d*h-e*g)
	if det > -eps && det < eps { // determinant is zero, matrix not invertible
		c1 = 0.0
		c2 = 0.0
		c3 = 0.0
		err = true
		return
	}

	// inverse 3x3 matrix
	// http://en.wikipedia.org/wiki/Matrix_inverse
	i00 := (e*k - f*h) / det
	i10 := (f*g - d*k) / det
	i20 := (d*h - e*g) / det
	i01 := (c*h - b*k) / det
	i11 := (a*k - c*g) / det
	i21 := (b*g - a*h) / det
	i02 := (b*f - c*e) / det
	i12 := (c*d - a*f) / det
	i22 := (a*e - b*d) / det

	c1 = i00*b0 + i01*b1 + i02*b2
	c2 = i10*b0 + i11*b1 + i12*b2
	c3 = i20*b0 + i21*b1 + i22*b2
	err = false
	return
}

// Par sum of squared residues of two parabolas fitted to the rows of a permuted similarity or distance matrix as the loss criterion.
func Par(sim Matrix64, p IntVector) float64 {
	if !sim.IsSymmetric() {
		panic("similarity matrix not symmetric")
	}
	n := p.Len()
	if sim.Rows() != n {
		panic("bad permutation vector length")
	}

	loss := 0.0
	rows := p.Len()
	cols := p.Len()

	v := NewVector64(cols)
	for i := 0; i < rows; i++ {
		// unload the permuted row to a vector
		for j := 0; j < cols; j++ {
			v[j] = sim[p[i]][p[j]]
		}

		// find position of the maximum
		mx := -inf
		pos := 0
		for j := 0; j < cols; j++ {
			if v[j] > mx {
				mx = v[j]
				pos = j
			}
		}

		if pos < 3 {
			// fit parabola to upper part
			m := cols - pos
			//			m:= cols-pos-1 // do not include pos
			w := NewVector64(m)
			for j := 0; j < m; j++ {

				// fmt.Println(m, j, pos)
				w[j] = v[j+pos]
				//w[j] = v[j+pos+1] // do not include pos
			}
			c1, c2, c3, err := parabolaFit(w)
			if !err {
				for j := 0; j < m; j++ {
					k := float64(j)
					x := w[j]
					y := c1 + c2*k + c3*k*k
					z := x - y
					//fmt.Println(x, yyy, z)
					z *= z
					loss += z
				}
			}
		} else if pos > n-4 {
			// fit parabola to lower part
			m := pos + 1
			//	m:= pos // do not include pos
			w := NewVector64(m)
			for j := 0; j < m; j++ {
				w[j] = v[j]
				//			w[j] = v[j+1] // do not include pos
			}
			c1, c2, c3, err := parabolaFit(w)
			if !err {
				for j := 0; j < m; j++ {
					k := float64(j)
					x := w[j]
					y := c1 + c2*k + c3*k*k
					z := x - y
					//fmt.Println(x, yyy, z)
					z *= z
					loss += z
				}
			}
		} else {

			// fit parabola to both: lower part
			m := pos + 1
			//m := pos // do not include pos

			w := NewVector64(m)
			for j := 0; j < m; j++ {
				w[j] = v[j]
				//w[j] = v[j+1] // do not include pos
			}

			//w.Print()
			c1, c2, c3, err := parabolaFit(w)
			if !err {
				for j := 0; j < m; j++ {
					k := float64(j)
					x := w[j]
					y := c1 + c2*k + c3*k*k
					z := x - y
					z *= z
					loss += z
				}
			}

			// fit parabola to upper part
			//fmt.Println("---")
			m = cols - pos
			//m = cols - pos - 1 // do not include pos
			w = NewVector64(m)
			for j := 0; j < m; j++ {

				// fmt.Println(m, j, pos)
				w[j] = v[j+pos]
				//w[j] = v[j+pos+1] // do not include pos
			}
			// w.Print()
			c1, c2, c3, err = parabolaFit(w)
			if !err {
				for j := 0; j < m; j++ {
					k := float64(j)
					x := w[j]
					y := c1 + c2*k + c3*k*k
					z := x - y
					//fmt.Println(x, yyy, z)
					z *= z
					loss += z
				}
			}

		}

	}
	return loss
}

// Are returns anti-Robinson events violation loss criterion of a permuted distance matrix (Hahsler et al. 2008: 4, Eq. 7, 8; Streng 1991; Streng and Schönfelder 1978; Chen 2002: 2, AR(i)  ; Wu et al. 2010: 773, AR_{n}).
func Are(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Ware returns weighted anti-Robinson events violation loss criterion of a permuted distance matrix (Hahsler et al. 2008; Streng 1991; Tien et al. 2008; Streng and Schönfelder 1978; Chen 2002:21, AR(s); Wu et al. 2010: 773, AR_{s}).
func Ware(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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

// Dware returns doubly weighted anti-Robinson violation loss criterion of a permuted distance matrix (Hahsler et al. 2008; Streng 1991; Tien et al. 2008; Chen 2002: 21, AR(w); Wu et al. 2010: 773, AR_{w}).
func Dware(dis Matrix64, p IntVector) float64 {
	if !dis.IsSymmetric() {
		panic("distance matrix not symmetric")
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
