// Copyright 2012 - 2014 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Objective functions for similarity matrices.

import (
	"math"
)

// Psis computes energy ψ(p) of the permuted similarity matrix according to Podani (1994); see  Miklos (2005), Eq. 4. It is a loss function.
func Psis(sim Matrix64, p IntVector) float64 {
	if !sim.IsSymmetric() {
		panic("similarity matrix not symmetric")
	}
	n := p.Len()
	if sim.Rows() != n {
		panic("dimensions not equal")
	}

	loss := 0.0
	rows := p.Len()
	cols := p.Len()
	for i := 0; i < p.Len(); i++ {
		for j := 0; j < p.Len(); j++ {
			x := sim[p[i]][p[j]]
			a := math.Abs(float64(cols*(i+1))/float64(rows) - float64(j+1))
			b := math.Abs(float64(rows*(j+1))/float64(cols) - float64(i+1))
			loss += x*a + b
		}
	}
	return loss
}

// Bers returns loss of the permuted matrix according to Kostopoulos & Goulermas(2014, in preparation). It is a loss function.
func Bers(sim Matrix64, p IntVector) float64 {
	if !sim.IsSymmetric() {
		panic("simtance matrix not symmetric")
	}
	n := p.Len()
	if sim.Rows() != n {
		panic("bad permutation vector length")
	}
	return Ber(sim, p, p)
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
