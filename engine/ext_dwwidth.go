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
	defer widths.Free()
	cuda.SetDWWidth(widths, M.Buffer(), DWHalfWidth)
	_dww := widths.HostCopy().Host()
	_M := M.Buffer().HostCopy().Host()[2][1015:1036]
	print(len(_dww[0]))
	print(len(_M))
	ret := math.Abs(c[X] * float64(cuda.Sum(widths)/float32(n[Y]*n[Z])))
	return ret
}
