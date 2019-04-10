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

func nextWallStep(mz [][]float32, path [][2]int) [2]int {

	last := path[len(path)-1]
	llast := path[len(path)-2]

	iMinAbsMz := [2]int{0, 0}
	minAbsMz := 1.0

	for _iy := last[0] - 1; _iy <= last[0]+1; _iy++ {

		iy := _iy % len(mz)

		// Don't iterate out of bounds; limit the range using intMax and intMin
		for ix := intMax(last[1]-1, 0); ix <= intMin(last[1]+1, len(mz[0])-1); ix++ {

			// If the path can be connected back onto itself in a loop, do so.
			if path[0][0] == iy && path[0][1] == ix {
				return [2]int{iy, ix}
			}

			// Don't include the location of the last point, or the next-to-last point
			if (last[0] == iy && last[1] == ix) || (llast[0] == iy && llast[1] == ix) {
				continue
			}

			// If _minAbsMz is closer to 0 than we've seen yet, store the associated indices
			_minAbsMz := math.Abs(float64(mz[iy][ix]))
			if _minAbsMz < minAbsMz {
				minAbsMz = _minAbsMz
				iMinAbsMz = [2]int{iy, ix}
			}
		}
	}
	return iMinAbsMz
}

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
