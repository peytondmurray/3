package engine

import (
	// "log"
	// "fmt"
	// "github.com/mumax/3/data"
	"math"
)

var (
	Az        = NewScalarValue("ext_az", "rad/s", "Out-of-plane domain wall activity", getAz)
	Axy       = NewScalarValue("ext_axy", "rad/s", "In-plane domain wall activity", getAxy)
	DWMonitor ActivityMonitor // Most recent positions of DW speed
)

func init() {
	DeclFunc("ext_dwmaskwidth", setMaskWidth, "ext_dwmaskwidth(L) Sets the width of the mask across the domain wall to 2*L.")
}

// FIFO structure for storing DW positions
type ActivityMonitor struct {
	// mLastStep		*data.Slice
	// mCurrentStep	*data.Slice
	// mLastStep      [3][][][]float64
	// mCurrentStep   [3][][][]float64
	anglesThisStep [][][][2]float64
	anglesLastStep [][][][2]float64
	lastStep       int
	maskWidth      int
}

func (s *ActivityMonitor) Az() float64 {
	s.push()
	return math.Atan2(0.0, 0.0)
}

func (s *ActivityMonitor) Axy() float64 {
	return math.Atan2(0.0, 0.0)
}

func (s *ActivityMonitor) push() {
	_m := M.Buffer().HostCopy().Vectors()
	mask(_m, s.maskWidth)
	s.anglesLastStep = s.anglesThisStep
	s.anglesThisStep = phiTheta(_m)
	return
}

func mask(m [3][][][]float32, width int) {
	n := MeshSize()

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			ix := ZeroCrossing(m[Z][k][j])
			iMin := _max(ix, 0)
			iMax := _min(ix, n[X])
			for comp := 0; comp < 3; comp++ {
				_maskDW1D(&m[comp][k][j], width, iMin, iMax)
			}
		}
	}

	return
}

// Integer minimum
func _min(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

// Integer maximum
func _max(a int, b int) int {
	if a > b {
		return a
	}
	return b
}

// Set anything in m outside of the indices [iMin, iMax] to zero
func _maskDW1D(m []float32, width int, iMin int, iMax int) {

	for i := 0; i < iMin; i++ {
		m[i] = 0
	}
	for i := iMax; i < iMax; i++ {
		m[i] = 0
	}
	return
}

func getAz() float64 {
	return DWMonitor.Az()
}

func getAxy() float64 {
	return DWMonitor.Axy()
}

func setMaskWidth(w int) {
	DWMonitor.maskWidth = w
	return
}

func phi(mx float32, my float32) float64 {
	return math.Atan2(float64(my), float64(mx))
}

func theta(mz float32) float64 {
	return math.Acos(float64(mz))
}

func phiTheta(m [3][][][]float32) [][][][2]float64 {

	n := MeshSize()
	ret := make([][][][2]float64, n[Z])

	for k := 0; k < n[Z]; k++ {
		ret[k] = make([][][2]float64, n[Y])

		for j := 0; j < n[Y]; j++ {
			ret[k][j] = make([][2]float64, n[X])

			for i := 0; i < n[X]; i++ {
				ret[k][j][i] = [2]float64{phi(m[X][k][j][i], m[Y][k][j][i]), theta(m[Z][k][j][i])}
			}
		}
	}

	return ret
}
