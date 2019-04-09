package engine

// type WallTree2D struct {
// 	mz [][]float32
// 	loc [2]int
// 	candidates [][2]int
// 	child *WallTree2D
// }

// func (w *WallTree2D) getWall() [][2]int {

// 	// If there are no children, simply return the location
// 	if w.child == nil {
// 		return [][2]int{{w.loc[0], w.loc[1]}}
// 	}

// 	// Otherwise, get the location of all the children, and append the location
// 	return append(w.child.getWall(), w.loc)
// }



// func plant(mz [][]float32) *WallTree2D {
// 	var w WallTree2D
// 	w.mz = mz
// 	w.loc = [2]int{0, zc(mz[0])}
// 	w.candidates = make([][2]int, 5)

// 	k := 0
// 	for i := 0; i <= 1; i++ {
// 		for j := -1; j <= 1; j++ {
// 			if i == 0 && j == 0 {
// 				continue
// 			}
// 			w.candidates[k] = [2]int{w.loc[0]+i, w.loc[1]+j}
// 		}
// 	}

// 	w.children = traceWall2D(mz, llast, last)

// 	return &w
// }


// func traceWall2D(mz [][]float32, path [][2]int) [][2]int {

// 	// If we are at the root, trace the path below
// 	if path == nil {
// 		return traceWall2D(mz, [1][2]int{0, zc(mz[0])})
// 	}




// 	return
// }


func traceWall2D(mz [][]float32) [][2]int {

	path := [][2]int{{0, zc(mz[0])}}
	for len(path) > 1 && isPathLoop(path)  {
		path = append(path, nextWallStep(mz, path))
	}
	return path
}

func isPathLoop(path [][2]int) bool {
	return path[len(path)-1][0] != path[0][0] && path[len(path)-1][1] != path[0][1]
}

func nextWallStep(mz [][]float32, path [][2]int) [2]int {

	last := path[len(path)-1]

	options := make([][2]int, 8)

	k := 0
	for iy := last[0]-1; iy <= last[0]+1; iy++ {

		// Don't iterate out of bounds; limit the range using intMax and intMin
		for ix := intMax(last[1]-1, 0); ix <= intMin(last[1]+1, len(mz[0])-1); ix++ {

			// If the path can be connected back onto itself in a loop, do so.
			if path[0][0] == iy % len(mz) && path[0][1] == ix {
				return [2]int{iy % len(mz), ix}
			}

			// Don't include the location of the last point
			if last[0] == iy && last[1] == ix {
				continue
			}

			options[k] = [2]int{iy, ix}
			k++
		}
	}

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








// func (w *WallTree2D) traverse() [][6]float64 {

// 	indices := make([][2]int, 1)
// 	for !wallComplete2D(indices) {
// 		indices = append(indices, w)
// 	}

// }

// func traverseWall2D() [][6]float64 {
// 	m := M.Buffer().HostCopy().Vectors()

// 	indices := make([][2]int, 1)		// Each element is a pair of indices [iy, ix] giving location of domain wall
// 	for !wallComplete2D(indices) {
// 		indices = append(indices, wallStep2D(m[2][0], indices))
// 	}

// 	return evaluateWall2D(indices)
// }

func zc(mz []float32) int {
	for i := len(mz)-1; i > 1; i-- {
		if mz[i]*mz[i-1] < 0 {
			return i-1
		}
	}
	panic("No zero crossing: domain wall not found.")
}

func wallComplete2D(indices [][2]int) bool {
	return len(indices) > 1 && indices[len(indices)-1][0] == indices[0][0] && indices[len(indices)-1][1] == indices[0][1]
}

func wallStep2D(mz [][]float32, indices [][2]int) [2]int {

	candidates := candidateSteps2D(mz, indices[len(indices)-2], indices[len(indices)-1])


	return
}

func candidateSteps2D(mz [][]float32, llast, last [2]int) [][2]int {

	candidates := make([][2]int, 7) // At most there are 7 possible points the next step can take

	k := 0
	for i := last[0]-1; i <= last[0]+1; i++ {
		for j := last[1]-1; j <= last[1]+1; j++ {

			// Exclude the points which are already in the domain wall
			if (i == last[0] && j == last[1]) || (i == llast[0] && j == llast[1]) {
				continue
			}

			// Wrap the y index if it is near the top or bottom
			if i < 0 {
				candidates[k][0] = len(mz) - i
			} else {
				candidates[k][0] = i % len(mz)
			}

			// Exclude points outside the simulation region in x
			if j < 0 || j > len(mz[0]) - 1 {
				continue
			} else {
				candidates[k][1] = j
			}
			k++
		}
	}

	return candidates[:k]

}