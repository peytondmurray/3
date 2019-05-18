package engine

import (
	// "fmt"
	"math"
)

func traceWall3D(mz [][][]float32) [][][2]int {
	path := make([][][2]int, len(mz))
	for i := range mz {
		path[i] = traceWall2D(mz[i])
	}
	return path
}

// Some extremely contrived domain wall configurations will not be identified correctly with this method. For example,
// if an extended bubble domain wall is traced which spans the entire simulation, but isn't the real domain wall,
// this function *might* trace that wall instead of the actual one you want.
func traceWall2D(mz [][]float32) [][2]int {

	buffer := [2]int{}
	path := make([][2]int, 1)

	// zcs := xRasterZCs(mz)
	zcs := allZCs(mz)
	zcs, buffer = pop(zcs)
	path[0] = buffer

	for len(zcs) > 0 {
		if pathSpansY(path, mz) {
			return path
		}
		zcs, buffer = popMin(zcs, path[len(path)-1], len(mz))

		// Append is fine here because path only ever expands.
		path = append(path, buffer)
	}
	return path
}

// Remove the first element from a, returning a tuple: (remaining list, first element)
func pop(a [][2]int) ([][2]int, [2]int) {
	return a[:len(a)-1], a[len(a)-1]
}

// Remove the element from zcs which is the nearest to b. Returns a tuple: (remaining zcs, closest element)
func popMin(zcs [][2]int, b [2]int, ny int) ([][2]int, [2]int) {
	dsts := allDistances(zcs, b, ny)
	minIndex := 0

	for i := range zcs {
		if dsts[i] < dsts[minIndex] {
			minIndex = i
		}
	}

	// Hideous bodge avoids unexpected behavior when using append(), which can modify the input parameters. In this
	// case, append modifies the array underlying the input zcs slice, which is needed later on in the next iteration.
	minZC := zcs[minIndex]
	remainingZCs := make([][2]int, len(zcs)-1)
	copy(remainingZCs[:minIndex], zcs[:minIndex])
	copy(remainingZCs[minIndex:], zcs[minIndex+1:])
	return remainingZCs, minZC
}

// Checks if path explicitly covers every y-index of the simulation region
func pathSpansY(path [][2]int, mz [][]float32) bool {
	for i := range mz {
		if !iyInPath(path, i) {
			return false
		}
	}
	return true
}

// Check whether a y-index iy is in the pathz
func iyInPath(path [][2]int, iy int) bool {
	for j := range path {
		if path[j][0] == iy {
			return true
		}
	}
	return false
}

// Return an array of {iy, ix} pairs of indices for each zero crossing in mz, scanning across rows
func xRasterZCs(mz [][]float32) [][2]int {
	zcs := make([][2]int, 0)
	for iy := range mz {
		for ix := len(mz[iy]) - 1; ix > 1; ix-- {
			if mz[iy][ix]*mz[iy][ix-1] < 0 {

				// Append is fine here, because the array is expanding, and we have no need to access the input array
				zcs = append(zcs, [2]int{iy, ix})
			}
		}
	}
	return zcs
}

// Return an array of {iy, ix} pairs of indices for each zero crossing in mz, scanning across columns
func yRasterZCs(mz [][]float32) [][2]int {
	zcs := make([][2]int, 0)
	for ix := range mz[0] {
		for iy := len(mz) - 1; iy > 1; iy-- {
			if mz[iy][ix]*mz[iy-1][ix] < 0 {

				// Append is fine here, because the array is expanding, and we have no need to access the input array
				zcs = append(zcs, [2]int{iy, ix})
			}
		}
	}
	return zcs
}

func allZCs(mz [][]float32) [][2]int {
	return uniqueZCs(xRasterZCs(mz), yRasterZCs(mz))
}

// Sort a by the cartesian distance between a[i] and b using insertion sort.
func insertionSort(a [][2]int, b [2]int, ny int) {

	dists := allDistances(a, b, ny)

	i := 1
	for i < len(a) {
		_a := a[i]
		_dist := dists[i]

		j := i - 1
		for j >= 0 && dists[j] > _dist {
			a[j+1] = a[j]
			dists[j+1] = dists[j]
			j--
		}
		a[j+1] = _a
		dists[j+1] = _dist
		i++
	}

	return
}

// Distances are wrapped in the y-direction.
func allDistances(a [][2]int, b [2]int, ny int) []float64 {

	ret := make([]float64, len(a))
	for i := range a {
		ret[i] = distance(a[i], b, ny)
	}
	return ret
}

func distance(a, b [2]int, ny int) float64 {
	dyAbs := math.Abs(float64(a[0] - b[0]))
	dy := math.Min(float64(ny)-dyAbs, dyAbs) // Wrap the distance around in the y-direction
	dx := float64(a[1] - b[1])
	return math.Sqrt(dx*dx + dy*dy)
}

func intMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// Given two lists of tuples, return a list of all the unique tuples.
func uniqueZCs(a, b [][2]int) [][2]int {
	zcs := make([][2]int, len(a)) // Allocate memory for the smallest possible list size
	for i := range a {
		zcs[i] = a[i]
	}
	for i := range b {
		if !isIn(zcs, b[i]) {
			zcs = append(zcs, b[i])
		}
	}
	return zcs
}

// Return true if tuple b is in list of tuples a; otherwise false.
func isIn(a [][2]int, b [2]int) bool {
	for i := range a {
		if b[0] == a[i][0] && b[1] == a[i][1] {
			return true
		}
	}
	return false
}
