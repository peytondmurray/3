package engine

import (
	"github.com/mumax/3/cuda"
	"math"
	"fmt"
)

var (
	DWWidth     = NewScalarValue("ext_dwwidth", "m", "Width of the domain wall, averaged along y.", setDWWidthCPU)
	// DWWidth     = NewScalarValue("ext_dwwidth", "m", "Width of the domain wall, averaged along y.", setDWWidth)
	DWHalfWidth = 10 //5 is usually a safe number. If you increase this, you get better domain wall fits but fitting will fail if Mz >= 1 (which can happen due to numerical error), which gets more likely the futher away from the domain wall you are.
)

func init() {
	DeclVar("dwhalfwidth", &DWHalfWidth, "Number of cells on either side of the domain wall to fit to find the width")
}

func setDWWidth() float64 {
	n := MeshSize()
	c := Mesh().CellSize()

	widths := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	cuda.SetDWWidth(widths, M.Buffer(), DWHalfWidth)
	avgWidth := c[X] * float64(cuda.Sum(widths)) / float64(n[Y]*n[Z])
	widths.Free()
	return avgWidth
}

func setDWWidthCPU() float64 {
	n := MeshSize()

	indices := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	defer indices.Free()
	cuda.SetDomainWallIndices(indices, M.Buffer())
	ret := avg(domainWallWidth(M.Buffer().Comp(Z).HostCopy().Scalars(), float32to2dint(indices.HostCopy().Scalars()), DWHalfWidth))


	print(fmt.Sprintf("CPU: %6.6e | GPU: %6.6e\n", ret, setDWWidth()))
	return ret
}

func avg(a [][]float64) float64 {
	sum := 0.0
	for i := range a {
		for j := range a[i] {
			sum += a[i][j]
		}
	}
	return sum/float64(len(a)*len(a[0]))
}

func domainWallWidth(mz [][][]float32, intPos [][]int, halfwidth int) [][]float64 {
	c := Mesh().CellSize()
	N := MeshSize()

	x := make([]float64, 2*halfwidth+1)
	for i := range x {
		x[i] = float64(i) * c[X]
	}

	width := make([][]float64, N[Z])
	for i := range mz {
		width[i] = make([]float64, N[Y])
		for j := range mz[i] {
			leftBound := intPos[i][j] - halfwidth
			rightBound := intPos[i][j] + halfwidth + 1
			width[i][j] = 1.0/fitSlope1D(x, atanh32to64(mz[i][j][leftBound:rightBound]))
		}
	}
	return width
}

func atanh32to64(m []float32) []float64 {
	ret := make([]float64, len(m))
	for i := range(m) {
		ret[i] = math.Atanh(float64(m[i]))
	}
	return ret
}

func sumSlice(s []float64) float64 {
	ret := 0.0
	for i := range s {
		ret += s[i]
	}
	return ret
}

func mulSlice(a, b []float64) []float64 {
	ret := make([]float64, len(a))
	for i := range a {
		ret[i] = a[i] * b[i]
	}
	return ret
}

func fitSlope1D(x, y []float64) float64 {
	N := float64(len(x))

	t4 := sumSlice(x)
	t1 := N * sumSlice(mulSlice(x, y)) / t4
	t2 := sumSlice(y)
	t3 := N * sumSlice(mulSlice(x, x)) / t4

	return (t1 - t2) / (t3 - t4)
}

func float32to2dint(a [][][]float32) [][]int {
	s := make([][]int, len(a))
	for i := range a {
		s[i] = make([]int, len(a[i]))
		for j := range a[i] {
			s[i][j] = int(a[i][j][0]+0.5)
		}
	}
	return s
}