package engine

import (
	// "log"
	// "fmt"
	// "github.com/mumax/3/data"
	"math"
)

var (
	Az	  = NewScalarValue("ext_az", "rad/s", "Out-of-plane domain wall activity", getAz)
	Axy   = NewScalarValue("ext_axy", "rad/s", "In-plane domain wall activity", getAxy)
	DWMonitor  ActivityMonitor // Most recent positions of DW speed
)

func init() {
	DeclFunc("ext_dwmaskwidth", setMaskWidth, "ext_dwmaskwidth(L) Sets the width of the mask across the domain wall to 2*L.")
}

// FIFO structure for storing DW positions
type ActivityMonitor struct {
	// mLastStep		*data.Slice
	// mCurrentStep	*data.Slice
	mLastStep		[3][][][]float64
	mCurrentStep	[3][][][]float64
	lastStep		int
	maskWidth		int
}

func (s *ActivityMonitor) Az() float64 {
	s.push()
	return math.Atan2(0.0, 0.0)
}

func (s *ActivityMonitor) Axy() float64 {
	return math.Atan2(0.0, 0.0)
}

func (s *ActivityMonitor) push() {
	s.anglesLastStep = s.anglesThisStep
	s.anglesThisStep = s.mask(M.Buffer().HostCopy().Vectors())
	return
}

func (s *ActivityMonitor) mask(m [3][][][]float32) [3][][][]float64 {
	n := MeshSize()

	mx :=

	ret := [][][]float64{}
	for comp:=0; comp<3; comp++ {
		ret[comp] = make([][][]float64, n[Z])
		for k:=0; k<n[Z]; k++ {
			ret[comp][k] = make([][]float64, n[Y])
			for j:=0; j<n[Y]; j++ {
				ret[comp][k][j] = make([]float64, n[X])
				for i:=0; i<n[X]; i++ {
					ret[k][j][i] = math.Atan2(float64(m[Y][k][j][i]), float64(m[X][k][j][i]))
				}
			}
		}
	}
	return ret
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

func phi(m [3][][][]float32) [][][]float64 {
	n := MeshSize()
	ret := make([][][]float64, n[Z])
	for k:=0; k<n[Z]; k++ {
		ret[k] = make([][]float64, n[Y])
		for j:=0; j<n[Y]; j++ {
			ret[k][j] = make([]float64, n[X])
			for i:=0; i<n[X]; i++ {
				ret[k][j][i] = math.Atan2(float64(m[Y][k][j][i]), float64(m[X][k][j][i]))
			}
		}
	}
	return ret
}

func theta(m [3][][][]float32) [][][]float64 {
	n := MeshSize()
	ret := make([][][]float64, n[Z])
	for k:=0; k<n[Z]; k++ {
		ret[k] = make([][]float64, n[Y])
		for j:=0; j<n[Y]; j++ {
			ret[k][j] = make([]float64, n[X])
			for i:=0; i<n[X]; i++ {
				ret[k][j][i] = math.Acos(float64(m[Z][k][j][i]))
			}
		}
	}
	return ret
}