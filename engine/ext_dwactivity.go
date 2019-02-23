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
	DeclFunc("ext_getphi", getPhi, "Get the current phi angle as a slice.")
	DeclFunc("ext_gettheta", getTheta, "Get the current theta angle as a slice.")
}

// DWActivityInit(w) sets the mask width to apply to the domain wall; only values of the magnetization within w cells
// of the domain wall are included in the domain wall activity
func DWActivityInit(w int) {
	DWMonitor.maskWidth = w
	DWMonitor.initialized = false
	return
}

type activityStack struct {
	phi         [][][]float64
	theta       [][][]float64
	t           float64
	maskWidth   int
	_AzCache    float64
	_AxyCache   float64
	initialized bool
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
	_phi, _theta := phiTheta(_m)
	if s.initialized {
		s._AzCache = calcAz(_theta, s.theta, Time-s.t, s.maskWidth, _m)
		s._AxyCache = calcAxy(_phi, s.phi, Time-s.t, s.maskWidth, _m)
	} else {
		s._AzCache = 0
		s._AxyCache = 0
		s.initialized = true
	}
	s.phi = _phi
	s.theta = _theta
	return
}

// Calculate the change in angle for two angles. This only makes sense to use if the input angles range from 0 to 2pi,
// so only use this for phi coordinate, NOT theta.
func deltaAngle(a, b float64) float64 {
	dA := a - b
	if dA < -math.Pi {
		return 2*math.Pi + dA
	} else if dA > math.Pi {
		return 2*math.Pi - dA
	}
	return dA
}

func getBandIndices(mz []float32, width int, xmin int, xmax int) (int, int) {
	ix := ZeroCrossing(mz)
	iMin := _max(ix-width, xmin)
	iMax := _min(ix+width, xmax)
	return iMin, iMax
}

func calcAz(thetaNew [][][]float64, thetaOld [][][]float64, dt float64, width int, m [3][][][]float32) float64 {

	n := MeshSize()
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			iMin, iMax := getBandIndices(m[Z][k][j], width, 0, n[X])
			// print(fmt.Sprintf("%i\t%i\n", iMin, iMax))
			for i := iMin; i <= iMax; i++ {
				ret += (thetaNew[k][j][i] - thetaOld[k][j][i]) / dt
			}
		}
	}
	return ret
}

func calcAxy(phiNew [][][]float64, phiOld [][][]float64, dt float64, width int, m [3][][][]float32) float64 {

	n := MeshSize()
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			iMin, iMax := getBandIndices(m[Z][k][j], width, 0, n[X])
			for i := iMin; i < iMax; i++ {
				_magnitude := inPlaneMagnitude(m[X][k][j][i], m[Y][k][j][i])
				ret += deltaAngle(phiNew[k][j][i], phiOld[k][j][i]) * _magnitude / dt
			}
		}
	}
	return ret
}

// Get in-plane magnitude
func inPlaneMagnitude(mx float32, my float32) float64 {
	return math.Sqrt(float64(mx*mx + my*my))
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

func phi(mx float32, my float32) float64 {
	return math.Atan2(float64(my), float64(mx))
}

func theta(mz float32) float64 {
	return math.Acos(float64(mz))
}

func getPhi() [][][]float64 {
	return DWMonitor.phi
}

func getTheta() [][][]float64 {
	return DWMonitor.theta
}
