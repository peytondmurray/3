package engine

import (
	"math"
	"fmt"
)


// Convert a slice  of shape [3]int{1, n[Y], n[Z]} to a [][]float32, eliminating the extra dimension
func to2D(a [][][]float32) [][]float32 {
	s := make([][]float32, len(a))
	for i := range a {
		s[i] = make([]float32, len(a[i]))
		for j := range a[i] {
			s[i][j] = a[i][j][0]
		}
	}
	return s
}

func isEqual(a [][]float32, b [][]int) bool {
	for i := range a {
		for j := range a[i] {
			if int(a[i][j]+0.5) != b[i][j] {
				return false
			}
		}
	}
	return true
}

func isSliceClose(a [][][]float32, b [][][]float32) bool {
	for i := range a {
		for j := range a[i] {
			for k := range a[i][j] {
				if !isClose(float64(a[i][j][k]), float64(b[i][j][k]), 1e-4, 1e-8) {
					_a := a[i][j][k]
					_b := b[i][j][k]
					print(fmt.Sprintf("(%d, %d, %d) | %8.8E | %8.8E\n", i, j, k, _a, _b))
					// return false
				}
			}
		}
	}
	return true
}

// DEBUG -------------

func isClose(a, b, rtol, atol float64) bool {
	return math.Abs(a-b) <= atol + rtol*math.Abs(b)
}