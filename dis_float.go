// Copyright 2012 - 2015 The Seriation Authors. All rights reserved. See the LICENSE file.

package ser

// Distance matrices, float.

import (
	"math"
)

func ManhattanD(dat Matrix64) (dis Matrix64) {
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			dis[i][j] = 0
			for k, _ := range row { // columns of "dat"
				dis[i][j] += math.Abs(dat[i][k] - dat[j][k])
			}
		}
	}
	return
}

func EuclideanD(dat Matrix64) (dis Matrix64) {
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum := 0.0
			for k, _ := range row { // columns of "dat"
				sum += (dat[i][k] - dat[j][k]) * (dat[i][k] - dat[j][k])
			}
			dis[i][j] = math.Sqrt(sum)
		}
	}
	return
}

/*
func CanberraD(dat Matrix64) (dis Matrix64) {
	rows := dat.Rows()
	dis = NewMatrix64(rows, rows)
	for i, row := range dat {
		for j, _ := range dat {
			sum := 0.0
			for k, _ := range row { // columns of "dat"
*/
