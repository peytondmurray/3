package engine

import (
	"github.com/mumax/3/cuda"
)

var (
	DWWidth     = NewScalarValue("ext_dwwidth", "m", "Width of the domain wall, averaged along y.", setDWWidth)
	DWHalfWidth = 10
)

func init() {
	DeclVar("dwwhalfwidth", &DWHalfWidth, "Number of cells on either side of the domain wall to fit to find the width")
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
