package engine

import (
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
// if an extened bubble domain wall is traced which spans the entire simulation, but isn't the real domain wall,
// this function *might* trace that wall instead of the actual one you want.
func traceWall2D(mz [][]float32) [][2]int {

	zcs := allZCs(mz)
	buffer := [2]int{}
	zcs, buffer = pop(zcs)
	path := make([][2]int, 1)
	path[0] = buffer

	for len(zcs) > 0 {
		if !isPathLoop(path) {

			// This can be done by doing an insertion sort. This is pretty efficient in terms of number of swaps
			// after the first sort, but sorting every time is probably not as efficient as just popping the nearest
			// zc from the list of zcs.
			// insertionSort(zcs, path[len(path)-1], len(mz))
			// path = append(path, pop(zcs))

			// Pop the zc nearest to the current zc.
			zcs, buffer = popMin(zcs, path[len(path)-1], len(mz))
			path = append(path, buffer)
		} else if pathSpansY(path, mz) {
			// The domain wall forms a loop, and all the y-indices in the simulation are covered. There are leftover
			// zcs, so there must be a bubble somewhere. Ignore.
			return path
		} else {
			// The domain wall forms a loop, but does not span every y-index in the simulation. Therefore a bubble must
			// have been traced, not the actual domain wall. Throw out all the ZCs in the current path
			// (which form the bubble), and start tracing a new path using a different ZC from the remaining zcs.
			zcs, buffer = pop(zcs)
			path = make([][2]int, 1)
			path[0] = buffer
		}
	}
	return path
}

func pop(a [][2]int) ([][2]int, [2]int) {
	ret := a[len(a)-1]
	a = a[:len(a)-1]
	return a, ret
}

func popMin(zcs [][2]int, b [2]int, ny int) ([][2]int, [2]int) {
	minDist := distance(zcs[0], b, ny)
	minIndex := 0

	for i := range zcs {
		_dist := distance(zcs[i], b, ny)

		if _dist < minDist {
			minDist = _dist
			minIndex = i
		}

	}
	minZC := zcs[minIndex]
	zcs = append(zcs[:minIndex], zcs[minIndex+1:]...)

	return zcs, minZC
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

// Checks if path forms a loop, path[0] == path[-1]
func isPathLoop(path [][2]int) bool {
	return len(path) > 1 && path[len(path)-1][0] == path[0][0] && path[len(path)-1][1] == path[0][1]
}

// Check whether a y-index iy is in the path
func iyInPath(path [][2]int, iy int) bool {
	for j := range path {
		if path[j][0] == iy {
			return true
		}
	}
	return false
}

// Return an array of {iy, ix} pairs of indices for each zero crossing in mz
func allZCs(mz [][]float32) [][2]int {
	zcs := make([][2]int, 0)
	for iy := range mz {
		for ix := len(mz[iy]) - 1; ix > 1; ix-- {
			if mz[iy][ix]*mz[iy][ix-1] < 0 {
				zcs = append(zcs, [2]int{iy, ix})
			}
		}
	}
	return zcs
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
	_dy := float64(a[0] - b[0])
	dy := math.Min(float64(ny)-math.Abs(_dy), math.Abs(_dy)) // Wrap the distance around in the y-direction
	dx := float64(a[1] - b[1])
	return math.Sqrt(dx*dx + dy*dy)
}

func intMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}
