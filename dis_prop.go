// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Distance matrices of compositional data (proportions, frequencies). Tested.

import (
	"math"
)

// Dat2Prop turns data martix to matrix of row proportions. Replaces zeros with a small number.
func Dat2Prop(dat Matrix64) {
	smallNum := 0.0
	for i, row := range dat {
		rowSum := 0.0
		for j, _ := range row { // columns of "dat"
			rowSum += dat[i][j]
		}
		sum := 0.0
		for j, _ := range row { // columns of "dat"
			dat[i][j] /= rowSum
			sum += dat[i][j]

			if dat[i][j] < smallNum {
				dat[i][j] = smallNum // replace zeros with a small number
			}

		}

		if math.Abs(sum-1) > 0.1 {
			panic("bad rowsum")
		}
	}
	return
}

// Logratio turns data martix to matrix of logarithms of row proportions. Replaces zeros with a small number.
func Logratio(dat Matrix64, col int) {
	if col < 0 {

		// get number of column with max nonzero entries
		rows, cols := dat.Dims()
		nonzero := 0
		maxNow := -1
		col = -1
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if dat[i][j] > 0 {
					nonzero++
				}
			}
			if nonzero > maxNow {
				col = j
				maxNow = nonzero
			}
		}
	}

	// turn data to proportions

	Dat2Prop(dat)

	x := dat.Clone()
	for i, row := range dat {
		for j, _ := range row { // columns of "dat"
			x[i][j] = -math.Log(dat[i][j] / dat[i][col])
		}
	}
	dat.CopyFrom(x)
	return
}

// Manly returns a Manly distance matrix.
func Manly(dat Matrix64) (dis Matrix64) {
	// Manly, B. F. (1994) Multivariate Statistical Methods. A primer. Second edition. Chapman & Hall, London.
	// Same as Prevosti (one-locus) genetic distance.
	// Prevosti A. (1974) La distancia genética entre poblaciones. Miscellanea Alcobé, 68, 109–118.
	// Prevosti A., Ocana J. and Alonso G. (1975) Distances between populations of Drosophila sub-
	// obscura, based on chromosome arrangements frequencies. Theoretical and Applied Genetics, 45, 231–241.
	// Equation from R: ade4:dist.prop.

	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum += math.Abs(p - q)
			}
			dis[i][j] = sum / 2
		}
	}
	return
}

// ManlyOverlap returns a matrix of the Overlap index of Manly.
func ManlyOverlap(dat Matrix64) (dis Matrix64) {
	// Manly, B. F. (1994) Multivariate Statistical Methods. A primer. Second edition. Chapman & Hall, London.
	// Equation from R: ade4:dist.prop.
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			sum3 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += p * q
				sum2 += p * p
				sum3 += q * q
			}
			dis[i][j] = 1 - sum1/(math.Sqrt(sum2)*math.Sqrt(sum3))
		}
	}
	return
}

// Nei returns a matrix of Nei (one locus) genetic distances.
func Nei(dat Matrix64) (dis Matrix64) {
	// Nei, M. (1972) Genetic distances between populations. The American Naturalist, 106, 283–292.
	// Equation from R: ade4:dist.prop.
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			sum3 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += p * q
				sum2 += p * p
				sum3 += q * q
			}
			/*
				if sum1*sum2*sum3 == 0 {
					panic("bad data")
				}
			*/
			dis[i][j] = -math.Log(sum1 / (math.Sqrt(sum2) * math.Sqrt(sum3)))
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
	}
	return
}

// Rogers72 returns a matrix of Rogers 1972 (one locus) genetic distances. Behaves as EuclideanD() on proportions.
func Rogers72(dat Matrix64) (dis Matrix64) {
	// Equation from R: ade4:dist.prop.
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum += (p - q) * (p - q)
			}
			dis[i][j] = math.Sqrt(sum / 2)
		}
	}
	return
}

// Edwards returns a matrix of Edwards 1971 (one locus) genetic distances. Edwards, A. W. F. (1971) Distance between populations on the basis of gene frequencies. Biometrics, 27, 873–881.
func Edwards(dat Matrix64) (dis Matrix64) {
	// Has problems, returns some NaNs. And is not always symmetric.
	// Def. from R:ade4: Edwards 1971 (one locus)}{\eqn{d_5=\sqrt{1-\sum_{i=1}^{K}{\sqrt{p_1 q_i}}}}{d5= sqrt (1 - (Sum(sqrt(p(i)q(i)))))}}
	// Equation from R: ade4:dist.prop.

	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum += math.Sqrt(p * q)
			}
			dis[i][j] = math.Sqrt(1 - sum)
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
		// zero out the diagonal (because the problems with NaNs)
		dis[i][i] = 0.0
	}
	return
}

/* Definitions
   1 = Manly d1 = sum|p(i) - q(i)|/2

   2 = Overlap index Manly d2 = 1 - Sum(p(i)q(i))/sqrt(Sum(p(i)^2))/sqrt(Sum(q(i)^2))

   3 = Rogers 1972 (one locus) d3 = sqrt(0.5*Sum(p(i)-q(i)^2))

   4 = Nei 1972 (one locus) d4 = -ln(Sum(p(i)q(i))/sqrt(Sum(p(i)^2))/sqrt(Sum(q(i)^2)))

   5 = Edwards 1971 (one locus) d5= sqrt (1 - (Sum(sqrt(p(i)q(i)))))

*/

// Reynolds returns a matrix of Reynolds 1971 (one locus) genetic distances.
func Reynolds(dat Matrix64) (dis Matrix64) {
	// Reynolds, J. B., B. S. Weir, and C. C. Cockerham. (1983) Estimation of the coancestry coefficient:
	// basis for a short-term genetic distance. Genetics, 105, 767–779.

	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += (p - q) * (p - q)
				sum2 += (p * q)
			}
			dis[i][j] = math.Sqrt(sum1 / (2 * (1 - sum2)))
		}
	}
	return
}

// Sforza returns a matrix of Cavalli-Sforza (one locus) genetic distances.
func Sforza(dat Matrix64) (dis Matrix64) {
	// L.L. Cavalli-Sforza, A.W.F. Edwards (1967). “Phylogenetic Analysis -Models and Estimation Procedures”.
	// The American Journal of Human Genetics 19 (3 Part I (May)).

	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += math.Sqrt(p * q)
			}
			dis[i][j] = 2 / math.Pi * math.Sqrt(2*(1-(sum1)))
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
		// zero out the main diagonal
		dis[i][i] = 0.0
	}
	return
}

// Nei83 returns a matrix of Nei-Tajima-Tateno-1983 (one locus) genetic distances.
func Nei83(dat Matrix64) (dis Matrix64) {
	// Nei, M., F. Tajima, & Y. Tateno (1983) Accuracy of estimated phylogenetic trees from molecular data.
	// II. Gene frequency data. J. Mol. Evol. 19:153-170.
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += math.Sqrt(p * q)
			}
			dis[i][j] = (1 - (sum1))
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
	}
	return
}

// Euclid returns a matrix of Euclidean (one locus) genetic distances.
func Euclid(dat Matrix64) (dis Matrix64) {
	// Nei, M. (1987). Molecular Evolutionary Genetics. (Chapter 9).
	// New York: Columbia University Press.
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += (p - q) * (p - q)
			}
			dis[i][j] = math.Sqrt(sum1)
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
	}
	return
}

// Nei73 returns a matrix of Nei’s 1973 minimum genetic distance (one locus).
func Nei73(dat Matrix64) (dis Matrix64) {
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum1 := 0.0
			sum2 := 0.0
			sum3 := 0.0
			for k, _ := range row { // columns of "dat"
				p := dat[i][k]
				q := dat[j][k]
				sum1 += p * p
				sum2 += q * q
				sum3 += p * q
			}
			dis[i][j] = (sum1*sum2)/2 - sum3
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
	}
	return
}

// Percentage computes a "Percentage dissimilarity" between the rows of a data matrix.
func Percentage(dat Matrix64) (dis Matrix64) {
	// Ref.: http://ordination.okstate.edu/distsim.htm
	// Gauch (1982).
	// Not tested yet. It is asymmetric, should it be??

	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
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
			dis[i][j] = 1 - 200*(sum1/sum2)
			if dis[i][j] == math.NaN() {
				dis[i][j] = 0.0
			}
		}
	}
	return
}
// Fixation index
// Prevosti = Manly
