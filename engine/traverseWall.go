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

func traceWall2D(mz [][]float32) [][2]int {

	path := [][2]int{{0, zc(mz[0])}}
	for len(path) > 1 && isPathLoop(path) {
		path = append(path, nextWallStep(mz, path))
	}
	return path
}

func isPathLoop(path [][2]int) bool {
	return path[len(path)-1][0] != path[0][0] && path[len(path)-1][1] != path[0][1]
}

// Search across adjacent rows. Find the next zero crossing closest in space to current location which is not
// already in the path.
func nextWallStep(mz [][]float32, path [][2]int) [2]int {

	zcs := allZCs(mz, path[len(path)-1])
	distances := distances(zcs, path[len(path)-1])

	zcs, distances = insertionSort(zcs, distances)


	return tup
}

// Search the path back to front for a tuple in the path.
func inPath(tup [2]int, path [][2]int) bool {
	for i := len(path)-1; i > 0; i-- {
		if tup[0] == path[i][0] && tup[1] == path[i][1] {
			return true
		}
	}
	return false
}

// Find all ZCs in rows adjacent to tup, which contains an {iy, ix} index pair of latest point in path
func allZCs(mz [][]float32, tup [2]int) [][2]int {

	k := 0
	tups := make([][2]int, 1)

	for i := tup[0]-1; i <= tup[0]+1; i++ {
		for j := len(mz[i])-1; j > 1; j-- {
			if mz[i][j]*mz[tup[0]][j-1] < 0 {
				tups[k] = [2]int{i, j}
			}
		}
	}
	return tups
}

// For a set of tuples containing {iy, ix} pairs, return the cartesian distances to loc
func distances(pts [][2]int, loc [2]int) []float64 {

	dists := make([]float64, len(pts))
	for i := 0; i < len(pts); i++ {
		dy := pts[i][0] - loc[0]
		dx := pts[i][1] - loc[1]
		dists[i] = math.Sqrt(float64(dy*dy) + float64(dx*dx))
	}
	return dists
}

// Sort a and b by the values in b using insertion sort.
func insertionSort(a [][2]int, b []float64) ([][2]int, []float64) {

	i := 1
	for i < len(b) {
		_a := a[i]
		_b := b[i]

		j := i-1
		for j >= 0 && b[j] > _b {
			a[j+1] = a[j]
			b[j+1] = b[j]
			j--
		}
		a[j+1] = _a
		b[j+1] = _b
		i++
	}

	return a, b
}




// func nextWallStep(mz [][]float32, path [][2]int) [2]int {

// 	last := path[len(path)-1]
// 	llast := path[len(path)-2]

// 	iMinAbsMz := [2]int{0, 0}
// 	minAbsMz := 1.0

// 	for _iy := last[0] - 1; _iy <= last[0]+1; _iy++ {

// 		iy := _iy % len(mz)

// 		// Don't iterate out of bounds; limit the range using intMax and intMin
// 		for ix := intMax(last[1]-1, 0); ix <= intMin(last[1]+1, len(mz[0])-1); ix++ {

// 			// If the path can be connected back onto itself in a loop, do so.
// 			if path[0][0] == iy && path[0][1] == ix {
// 				return [2]int{iy, ix}
// 			}

// 			// Don't include the location of the last point, or the next-to-last point
// 			if (last[0] == iy && last[1] == ix) || (llast[0] == iy && llast[1] == ix) {
// 				continue
// 			}

// 			// If _minAbsMz is closer to 0 than we've seen yet, store the associated indices
// 			_minAbsMz := math.Abs(float64(mz[iy][ix]))
// 			if _minAbsMz < minAbsMz {
// 				minAbsMz = _minAbsMz
// 				iMinAbsMz = [2]int{iy, ix}
// 			}
// 		}
// 	}
// 	return iMinAbsMz
// }

func intMax(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func intMin(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func zc(mz []float32) int {
	for i := len(mz) - 1; i > 1; i-- {
		if mz[i]*mz[i-1] < 0 {
			return i - 1
		}
	}
	panic("No zero crossing: domain wall not found.")
}
