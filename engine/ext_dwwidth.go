package engine

import (
	"github.com/mumax/3/cuda"
	"math"
	// "github.com/mumax/3/data"
)

var (
	DWWidth     = NewScalarValue("ext_dwwidth", "m", "Width of the domain wall, averaged along y.", setDWWidth)
	DWHalfWidth = 10
)

func init() {
	DeclVar("ext_dwwhalfwidth", &DWHalfWidth, "Number of cells on either side of the domain wall to fit to find the width")
}

func setDWWidth() float64 {
	n := MeshSize()
	c := Mesh().CellSize()

	widths := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	cuda.SetDWWidth(widths, M.Buffer(), DWHalfWidth)
	avgWidth := c[X] * float64(cuda.Sum(widths)) / float64(n[Y]*n[Z])
	widths.Free()
	// cpuAvgWidth := avg(tanhFit(M.Comp(Z).HostCopy().Scalars(), DWHalfWidth))

	return avgWidth
}

// func toMeters(x [][]float32) [][]float32 {
// 	c := Mesh().CellSize()
// 	ret := make([][]float32, len(x))
// 	for i := range x {
// 		ret[i] = make([]float32, len(x[i]))
// 		for j := range x[i] {
// 			ret[i][j] = float32(c[X]) * x[i][j]
// 		}
// 	}
// 	return ret
// }

// func to2D(x [][][]float32) [][]float32 {
// 	c := Mesh().CellSize()
// 	ret := make([][]float32, len(x))
// 	for i := range x {
// 		ret[i] = make([]float32, len(x[i]))
// 		for j := range x[i] {
// 			ret[i][j] = x[i][j][0] * float32(c[X])
// 		}
// 	}
// 	return ret
// }

func avg(a [][]float64) float64 {
	ret := 0.0
	for i := range a {
		for j := range a[i] {
			ret += a[i][j]
		}
	}
	return ret / float64(len(a)*len(a[0]))
}

func tanhFit(mz [][][]float32, hw int) [][]float64 {

	c := Mesh().CellSize()
	n := MeshSize()
	ret := make([][]float64, n[Z])
	for iz := range mz {
		ret[iz] = make([]float64, n[Y])
		for iy := range mz[iz] {
			izc := zc(mz[iz][iy])
			ret[iz][iy] = (1.0 / fs(mz[iz][iy], izc, hw)) * c[X]
		}
	}

	return ret

}

func fs(mz []float32, izc int, hw int) float64 {

	N := float64(2*hw + 1) //# of cells running across the domain wall

	x := make([]float64, 2*hw+1)
	y := make([]float64, 2*hw+1)
	for i := range x {
		x[i] = float64(i)
		y[i] = math.Atanh(float64(mz[izc-hw+i]))
	}

	t1 := N * sum(mul(x, y))
	t2 := sum(y) * (N - 1) * N / 2
	t3 := N * (N - 1) * N * (2*N - 1) / 6.0
	t4 := ((N - 1) * N / 2) * ((N - 1) * N / 2)

	// t2 := outersum(x, y)
	// t3 := N * sum(mul(x, x))
	// t4 := outersum(x, x)

	return (t1 - t2) / (t3 - t4)
}

// func indexMulSumAtanh(x []float32, izc, hw int) float64 {
// 	ret := float64(0)
// 	k := 0
// 	for i := izc - hw; i < izc+hw+1; i++ {
// 		// ret += math.Atanh(float64(x[i]))
// 		ret += float64(k) * math.Atanh(float64(x[i]))
// 		// ret += float64(i) * math.Atanh(float64(x[i]))
// 		k++
// 	}
// 	return ret
// }

// func sumAtanh(x []float32, izc, hw int) float64 {
// 	ret := 0.0
// 	for i := izc - hw; i < izc+hw+1; i++ {
// 		ret += math.Atanh(float64(x[i]))
// 	}
// 	return ret
// }

func zc(arr []float32) int {
	for i := len(arr) - 1; i > 1; i-- {
		if arr[i]*arr[i-1] < 0 {
			return i - 1
		}
	}
	panic("hey")
}
