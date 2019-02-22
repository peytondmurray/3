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
	DWMonitor activityStack // Most recent positions of DW speed
)

func init() {
	DeclFunc("ext_dwactivityinit", DWActivityInit, "ext_dwactivityinit(a, w) sets the order of the DW activity calculation to a and the mask width to w.")
}

// DWActivityInit(w) sets the mask width to apply to the domain wall; only values of the magnetization within w cells
// of the domain wall are included in the domain wall activity
func DWActivityInit(w int) {
	DWMonitor.maskWidth = w
	return
}

type activityStack struct {
	phi       [][][]float64
	theta     [][][]float64
	t         float64
	maskWidth int
	_AzCache  float64
	_AxyCache float64
}

func (s *activityStack) Az() float64 {
	if s.t != Time {
		s.cacheAzxy()
	}
	return s._AzCache
}

func (s *activityStack) Axy() float64 {
	if s.t != Time {
		s.cacheAzxy()
	}
	return s._AxyCache
}

func getAz() float64 {
	return DWMonitor.Az()
}

func getAxy() float64 {
	return DWMonitor.Axy()
}

func (s *activityStack) cacheAzxy() {
	_m := M.Buffer().HostCopy().Vectors()
	maskDW(_m, s.maskWidth)
	_phi, _theta := phiTheta(_m)
	s._AzCache = calcAz(_theta, s.theta, Time, s.t)
	s._AxyCache = calcAxy(_phi, s.phi, Time, s.t, _m)
	return
}

func calcAz(thetaNew [][][]float64, thetaOld [][][]float64, tNew float64, tOld float64) float64 {

	n := MeshSize()
	dt := tNew - tOld
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < n[X]; i++ {
				ret += (thetaNew[k][j][i] - thetaOld[k][j][i]) / dt
			}
		}
	}
	return ret
}

func calcAxy(phiNew [][][]float64, phiOld [][][]float64, tNew float64, tOld float64, mNew [3][][][]float32) float64 {

	n := MeshSize()
	dt := tNew - tOld
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < n[X]; i++ {
				_magnitude := inPlaneMagnitude(mNew[X][k][j][i], mNew[Y][k][j][i])
				ret += (phiNew[k][j][i] - phiOld[k][j][i]) * _magnitude / dt
			}
		}
	}
	return ret
}

func inPlaneMagnitude(mx float32, my float32) float64 {
	return math.Pow(float64(mx*mx+my*my), 0.5)
}

func maskDW(m [3][][][]float32, width int) {
	n := MeshSize()

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			ix := ZeroCrossing(m[Z][k][j])
			iMin := _max(ix-width, 0)
			iMax := _min(ix+width, n[X])
			for comp := 0; comp < 3; comp++ {
				_maskDW1D(m[comp][k][j], iMin, iMax)
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
func _maskDW1D(m []float32, iMin int, iMax int) {

	for i := 0; i < iMin; i++ {
		m[i] = 0
	}
	for i := iMax; i < len(m); i++ {
		m[i] = 0
	}
	return
}

func phi(mx float32, my float32) float64 {
	return math.Atan2(float64(my), float64(mx))
}

func theta(mz float32) float64 {
	return math.Acos(float64(mz))
}

func phiTheta(m [3][][][]float32) ([][][]float64, [][][]float64) {

	n := MeshSize()
	_phi := make([][][]float64, n[Z])
	_theta := make([][][]float64, n[Z])

	for k := 0; k < n[Z]; k++ {
		_phi[k] = make([][]float64, n[Y])
		_theta[k] = make([][]float64, n[Y])

		for j := 0; j < n[Y]; j++ {
			_phi[k][j] = make([]float64, n[X])
			_theta[k][j] = make([]float64, n[X])

			for i := 0; i < n[X]; i++ {
				_phi[k][j][i] = phi(m[X][k][j][i], m[Y][k][j][i])
				_theta[k][j][i] = theta(m[Z][k][j][i])
			}
		}
	}
	return _phi, _theta
}
