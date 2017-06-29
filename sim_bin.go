// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Similarity matrices, binary.
// Thanks to Gerald Jurasinski <gerald.jurasinski@uni-rostock.de> 2012-03-29, for overview (R simba).
// Not tested yet.

import (
	"math"
)

func abcd(dat Matrix64, i, h int, norm bool) (a, b, c, d, n, n1, n2 float64) {
	a, b, c, d = 0, 0, 0, 0
	cols := dat.Cols()
	for j := 0; j < cols; j++ { // columns of "dat"
		if dat[i][j] > 0 && dat[h][j] > 0 {
			a++
		} else if dat[i][j] > 0 && dat[h][j] == 0 {
			b++
		} else if dat[i][j] == 0 && dat[h][j] > 0 {
			c++
		} else {
			d++
		}
	}

	if norm {
		a = a / (a + b + c)
		b = b / (a + b + c)
		c = c / (a + b + c)
	}

	n = a + b + c + d
	if a+b <= a+c {
		n1 = a + b
	} else {
		n1 = a + c
	}

	if a+b > a+c {
		n2 = a + b
	} else {
		n2 = a + c
	}

	return
}

//// Computable asymmetric indices.

// Soerensen (1948)
func Soerensen(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat { // rows of "dat"
		for h, _ := range dat { // rows of "dat"
			// a, b, c, d, n, n1, n2:= dat.abcd(i,h,norm)
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 2 * a / (2*a + b + c)
		}
	}
	return
}

// Jaccard (1912)
func Jaccard(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (a + b + c)
		}
	}
	return
}

// Ochiai (1957), Shi (1993)
func Ochiai(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / math.Sqrt((a+b)*(a+c))
		}
	}
	return
}

// Mountford (1962), Shi (1993)
func Mountford(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 2 * a / (a*(b+c) + 2*(b+c))
		}
	}
	return
}

// Whittaker (1960), Magurran (1988)
func Whittaker(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (2 * (a + b + c) / (2*a + b + c)) - 1
		}
	}
	return
}

// Lande (1996)
func Lande(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			_, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (b + c) / 2
		}
	}
	return
}

// Wilson & Shmida (1984)
func WilsonShmida(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (b + c) / (2*a + b + c)
		}
	}
	return
}

// Colwell & Coddington (1948), Gaston et al. (2001)
func CoCoGaston(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (b + c) / (a + b + c)
		}
	}
	return
}

// Magurran (1988)
func Magurran(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (2*a + b + c) * (1 - a/(a+b+c))
		}
	}
	return
}

// Harrison et al. (1992), Koleff et al. (2003)
func Harrison(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = min(b, c) / (max(b, c) + a)
		}
	}
	return
}

// Cody (1993)
func Cody(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 1 - (a * (2*a + b + c) / (2 * (a + b) * (a + c)))
		}
	}
	return
}

// Williams (1996), Koleff et al. (2003)
func Williams1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = min(b, c) / (a + b + c)
		}
	}
	return
}

// Williams (1996), Koleff et al. (2003)
func Williams2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 2 * (b*c + 1) / ((a+b+c)*(a+b+c) - (a + b + c))
		}
	}
	return
}

// Harte & Kinzig (1997), Koleff et al. (2003)
func Harte(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (2*a + b + c) * (1 - a/(a+b+c))
		}
	}
	return
}

// Simpson (1949), Koleff et al. (2003)
func Simpson1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = min(b, c) / (min(b, c) + a)
		}
	}
	return
}

// Lennon et al. (2001), Koleff et al. (2003)
func Lennon1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 2 * math.Abs(b-c) / (2*a + b + c)
		}
	}
	return
}

// Weiher & Boylen (1994)
func Weiher(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			_, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (b + c)
		}
	}
	return
}

// Ruggiero et al. (1998), Koleff et al. (2003)
func Ruggiero(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, _, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (a + c)
		}
	}
	return
}

// Lennon et al. (2001), Koleff et al. (2003)
func Lennon2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 1 - (math.Log((2*a+b+c)/(a+b+c)) / math.Log(2.0))
		}
	}
	return
}

// Routledge (1977), Magurran (1988)
func Routledge1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = ((a+b+c)*(a+b+c)/((a+b+c)*(a+b+c)-(2*b+c)) - 1)
		}
	}
	return
}

// Routledge, 1977; Koleff et al. 2003
func Routledge2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = math.Log(2*a+b+c) - 2*a*math.Log(2)/(2*a+b+c) - ((a+b)*math.Log(a+b)+(a+c)*math.Log(a+c))/(2*a+b+c)
		}
	}
	return
}

// Routledge (1977)
func Routledge3(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = math.Exp(math.Log(2*a+b+c)-2*a*math.Log(2)/(2*a+b+c)-((a+b)*math.Log(a+b)+(a+c)*math.Log(a+c))/(2*a+b+c)) - 1
		}
	}
	return
}

// Sokal & Sneath (1963)
func Sokal1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (a + 2*(b+c))
		}
	}
	return
}

// Association index of Dice (1945), Wolda (1981)
func Dice(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / min((a+b), (a+c))
		}
	}
	return
}

// Oosting (1956), Southwood (1978)
func Kulczinsky1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (b + c)
		}
	}
	return
}

// Kulczinsky2; Oosting (1956), Southwood (1978)
func Kulczinsky2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a * (2*a + b + c) / (2 * (a + b) * (a + c))
		}
	}
	return
}

// McConnagh; Hubalek (1982)
func McConnagh(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a*a - b*c) / ((a + b) * (a + c))
		}
	}
	return
}

// Simpson (1960), Shi (1993)
func Simpson2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, _, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (a + b)
		}
	}
	return
}

// Legendre & Legendre (1998)
func Legendre2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 3 * a / (3*a + b + c)
		}
	}
	return
}

// Fager (1957), Shi (1993)
func Fager(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, _, _, _, _, n1, n2 := abcd(dat, i, h, norm)
			s[i][h] = a/math.Sqrt(n1*n2) - 1/(2*math.Sqrt(n2))
		}
	}
	return
}

// van der Maarel (1969)
func Maarel(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 2*a - (b+c)/(2*a+b+c)
		}
	}
	return
}

// Lamont and Grant (1979)
func Lamont(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (2*a + b + c)
		}
	}
	return
}

// Johnson (1971)
func Johnson1(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, _, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (2 * b)
		}
	}
	return
}

// Sorgenfrei (1959)
func Sorgenfrei(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a * a / ((a + b) * (a + c))
		}
	}
	return
}

// Johnson (1967)
func Johnson2(dat Matrix64, norm bool) (s Matrix64) {
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, _, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a/(a+b) + a/(a+c)
		}
	}
	return
}

//// Computable symmetric indices (including unshared species).

// Mean Manhattan, Legendre & Legendre (1998)
func Manhattan(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (b + c) / (a + b + c + d)
		}
	}
	return
}

// Sokal & Michener 1958
func SimpleMatching(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a + d) / (a + b + c + d)
		}
	}
	return
}

// Clifford & Stevenson (1975)
func Margaleff(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a * (a + b + c + d) / ((a + b) * (a + c))
		}
	}
	return
}

// Phi of Pearson, Gower & Legendre (1986), Yule (1912)
func Yule1(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a*d - b*c) / math.Sqrt((a+b)*(a+c)*(d+b)*(d+c))
		}
	}
	return
}

// Rogers & Tanimoto (1960), Gower & Legendre (1986)
func Rogers(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a + d) / (a + 2*(b+c) + d)
		}
	}
	return
}

// Baroni-Urbani & Buser (1976), Wolda (1981)
func Baroni(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (math.Sqrt(a*d) + a) / (math.Sqrt(a*d) + a + b + c)
		}
	}
	return
}

// Holliday et al. (2002), Ellis et al. (1993)
func Dennis(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a*d - b*c) / math.Sqrt((a+b+c+d)*(a+b)*(a+c))
		}
	}
	return
}

// Holliday et al. (2002), Ellis et al. (1993)
func Fossum(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = ((a + b + c + d) * -a * a / 4) / ((a + b) * (a + c))
			//((a+b+c+d) * (-1 * ((a/2)^2))) / ((a+b)*(a+c))
		}
	}
	return
}

// Gower & Legendre (1986)
func Gower(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a - (b + c) + d) / (a + b + c + d)
		}
	}
	return
}

// Gower & Legendre (1986), Russell/Rao in Ellis et al. (1993)
func Legendre(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = a / (a + b + c + d)
		}
	}
	return
}

// Sokal & Sneath (1963)
func Sokal2(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a * d) / math.Sqrt((a+b)*(a+c)*(d+b)*(d+c))
		}
	}
	return
}

// Sokal & Sneath (1963)
func Sokal3(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (2*a + 2*d) / (2*a + 2*d + b + c)
		}
	}
	return
}

// Sokal & Sneath (1963)
func Sokal4(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a + d) / (b + c)
		}
	}
	return
}

// Stiles (1946)
func Stiles(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			x := a + b + c + d
			y := math.Abs(a*d-b*c) - (a+b+c+d)/2
			y = y * y
			z := (a + b) * (a + c) * (b + d) * (c + d)
			s[i][h] = math.Log(x * y / z)
		}
	}
	return
}

// Yule & Kendall (1973)
func Yule2(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a*d - b*c) / (a*d + b*c)
		}
	}
	return
}

// Michael (1920), Shi (1993)
func TEMPLATE(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = 4 * (a*d - b*c) / ((a+d)*(a+d) + (b+c)*(b+c))
		}
	}
	return
}

// Hamann (1961)
func Hamann(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, n, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a + d - b - c) / n
		}
	}
	return
}

// Forbes (1925), Shi (1993)
func Forbes(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, _, _, _, n, n1, n2 := abcd(dat, i, h, norm)
			s[i][h] = (a*n - 2*n2) / (n*n1 - 2*n2)
		}
	}
	return
}

// Yule & Kendall (1950)
func ChiSquare(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			x := a + b + c + d
			y := (a*d - b*c) * (a*d - b*c)
			z := (a + b) * (a + c) * (b + d) * (c + d)
			s[i][h] = x * y / z
		}
	}
	return
}

// Peirce (1884)
func Peirce(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			s[i][h] = (a*d - b*c) / ((a + c) * (b + d))
		}
	}
	return
}

// Eyraud (1936) in Shi (1993)
func Eyraud(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			y := a - (a+b)*(a+c)
			z := (a + b) * (a + c) * (b + d) * (c + d)
			s[i][h] = y / z
		}
	}
	return
}

// Mean Euclidean in Ellis et al. (1993)
func BinEuclidean(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			x := math.Sqrt(b + c)
			y := a + b + c + d
			s[i][h] = x / y
		}
	}
	return
}

// Divergence in Ellis et al. (1993)
func Ellis(dat Matrix64) (s Matrix64) {
	norm := false
	rows := dat.Rows()
	s = NewMatrix64(rows, rows)
	for i, _ := range dat {
		for h, _ := range dat {
			a, b, c, d, _, _, _ := abcd(dat, i, h, norm)
			x := math.Sqrt(b + c)
			y := math.Sqrt(a + b + c + d)
			s[i][h] = x / y
		}
	}
	return
}

/*
=== References ===

Albatineh, A. N., Niewiadomska-Bugaj, M. & Mihalko, D. (2006) On Similarity Indices and Correction for Chance Agreement. Journal of Classification V23: 301-313.
Baroni-Urbani, C. & Buser, M. W. (1976) Similarity of Binary Data. Systematic Zoology 25: 251-259.
Clifford, H. T. & Stephenson, W. (1975) An introduction to numerical classification. Academic Press, New York, San Francisco, London.
Cody, M. L. (1993) Bird diversity components within and between habitats in Australia. - In: Ricklefs, R. E. & Schluter, D. (eds.), Species Diversity in Ecological Communities: historical and geographical perspectives, pp. 147-158, University of Chicago Press, Chicago.
Colwell, R. K. & Coddington, J. A. (1994) Estimating terrestrial biodiversity through extrapolation.
Philosophical Transactions of the Royal Society of London Series B-Biological Sciences 345: 101-118.
Dice, L. R. (1945) Measures of the amount of ecological association between species. Ecology 26: 297-302.
Ellis, D., Furner-Hines, J., Willett, P. (1993) Measuring the degree of similarity between objects in text retrieval systems. Perspectives in Information Management 3(2): 128-149
Fager, E. W. (1957) Determination and analysis of recurrent groups. Ecology 38: 586-595.
Faith, D. P., Minchin, P. R. & Belbin, L. (1987) Compositional dissimilarity as a robust measure of ecological distance. Plant Ecology 69: 57-68.
Gaston, K. J., Rodrigues, A. S. L., van Rensburg, B. J., Koleff, P. & Chown, S. L. (2001) Complementary representation and zones of ecological transition. Ecology Letters 4: 4-9.
Gower, J. C. & Legendre, P. (1986) Metric and Euclidean properties of dissimilarity coefficients. Journal of Classification 3: 5-48.
Hajdu, L. J. (1981) Graphical comparison of resemblance measures in phytosociology. Plant Ecology V48: 47-59.
Harrison, S., Ross, S. J. & Lawton, J. H. (1992) Beta diversity on geographic gradients in Britain. Journal of Animal Ecology 61: 151-158.
Harte, J. & Kinzig, A. (1997) On the implications of species-area relationships for endemism, spatial turnover and food web patterns. Oikos 80.
Holliday, J. D., Hu, C.-Y. & Willett, P. (2002) Grouping of Coefficients for the Calculation of InterMolecular Similarity and Dissimilarity using 2D Fragment Bit-Strings. Combinatorial Chemistry & High Throughput Screening 5: 155-166.
Hubalek, Z. (1982) Coefficients of association and similarity, based on binary (presence-absence) data: An evaluation. Biological Reviews of the Cambridge Philosophical Society 57: 669-689.
Huhta, V. (1979) Evaluation of different similarity indices as measures of succession in arthropod communities of the forest floor after clear-cutting. Oecologia V41: 11-23.
Jaccard, P. (1901) Etude comparative de la distribution florale d’une portion des Alpes et du Jura. Bulletin de la Societé Vaudoise des Sciences Naturelles 37: 547-579.
Jaccard, P. (1912) The distribution of the flora of the alpine zone. New Phytologist 11: 37-50.
Johnson, J. G. (1971) A quantitative approach to faunal province analysis. American Journal of Science 270: 257-280.
Johnson, S. C. (1967) Hierarchical clustering schemes. Psychometrika 32: 241-254.
Jurasinski, G. & Beierkuhnlein, C. (2006) Spatial patterns of biodiversity - assessing vegetation using hexagonal grids. Proceedings of the Royal Irish Academy - Biology and Environment 106B: 401-411.
Jurasinski, G. & Beierkuhnlein, C. (submitted) Distance decay and non-stationarity in a semi-arid Mediterranean ecosystem. Journal of Vegetation Science.
Koleff, P., Gaston, K. J. & Lennon, J. J. (2003) Measuring beta diversity for presence-absence data. Journal of Animal Ecology 72: 367-382.
Lamont, B. B. & Grant, K. J. (1979) A comparison of twenty-one measures of site dissimilarity. - In: Orlóci, L., Rao, C. R. & Stiteler, W. M. (eds.), Multivariate Methods in Ecological Work, pp. 101-126, Int. Coop. Publ. House, Fairland, MD
Lande, R. (1996) Statistics and partitioning of species diversity and similarity along multiple communities. Oikos 76: 25-39.
Legendre, P. & Legendre, L. (1998) Numerical Ecology. Elsevier, Amsterdam.
Lennon, J. J., Koleff, P., Greenwood, J. J. D. & Gaston, K. J. (2001) The geographical structure of British bird distributions: diversity, spatial turnover and scale. J Anim Ecology 70: 966-979.
Magurran, A. E. (1988) Ecological Diversity and its Measurement. Chapman & Hall, London.
Mountford, M. D. (1962) An index of similarity and its application to classification problems. - In: Murphy, P. W. (ed.) Progress in Soil Zoology, pp. 43-50, Butterworths.
Ochiai, A. (1957) Zoogeographical studies on the soleoid fishes found in Japan and its neighbouring regions. Bulletin of the Japanese Society of Fisheries Science 22(9): pp. 526-530.
Oosting, H. J. (1956) The study of plant communities: an introduction to plant ecology. W. H. Freeman, San Francisco.
Rogers, D. J. & Tanimoto, T. T. (1960) A computer program for classifying plants. Science 132: 1115-1118.
Routledge, R. D. (1977) On Whittaker’s components of diversity. Ecology 58: 1120-1127.
Ruggiero, A., Lawton, J. H. & Blackburn, T. M. (1998) The geographic ranges of mammalian species in South America: spatial patterns in environmental resistance and anisotropy. Journal of Biogeography 25: 1093-1103.
Shi, G. R. (1993) Multivariate data analysis in palaeoecology and palaeobiogeography–a review. Palaeogeography, Palaeoclimatology, Palaeoecology 105: 199-234.
Simpson, E. H. (1949) The measurement of diversity. Nature 163: 688.
Simpson, G. G. (1960) Notes on the measurement of faunal resemblance. American Journal of Science 258-A: 300-311.
Sokal, R. R. & Michener, C. D. (1958) A statistical method for evaluating systematic relationships. University of Kansas Science Bulletin 38: 1409-1438.
Sokal, R. R. & Sneath, P. H. A. (1963) Principles of numerical taxonomy. W. H. Freeman, San Francisco.
Sørensen, T. (1948) A method of establishing groups of equal amplitude in plant sociology based on similarity of species content. Biologiske Skrifter 5: 1-34.
Sorgenfrei, T. (1959) Molluscan assemblages from the marine middle Miocene of South Jutland and their environments. Danmark Geologiske Undersøgelse. Serie 2 79: 403-408.
Southwood, T. S. (1978) Ecological Methods. Chapman and Hall, London.
Stiles (1961) The association factor in information retrieval. Journal of the Association for Computing Machinery. 8: 271-279
Weiher, E. & Boylen, C. W. (1994) Patterns and prediction of alpha and beta diversity of aquatic plants in Adirondack (New York) lakes. Canadian Journal of Botany - Revue Canadienne De Botanique 72: 1797-1804.
Whittaker, R. H. (1960) Vegetation of the Siskiyou Mountains, Orgeon and California. Ecological Monographs 30: 279-338.
Williams, P. H. (1996) Mapping variations in the strength and breadth of biogeographic transition zones using species turnover. Proceedings of the Royal Society of London Series B-Biological Sciences 263: 579-588.
Williams, P. H., Klerk, H. M. & Crowe, T. M. (1999) Interpreting biogeographical boundaries among Afrotropical birds: spatial patterns in richness gradients and species replacement. J Biogeography 26: 459-474.
Wilson, M. V. & Shmida, A. (1984) Measuring beta-diversity with presence-absence data. Journal of Ecology 72: 1055-1064.
Wolda, H. (1981) Similarity indices, sample size and diversity. Oecologia 50: 296-302.
Yule, G. U. & Kendall, M. G. (1973) An introduction to the theory of statistics. Griffin, London.
Yule, G. U. (1912) On the methods of measuring association between two attributes. Journal of the Royal Statistical Society 75(6): 579-642

*/
