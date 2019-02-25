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
	DeclFunc("ext_getWindowShift", GetWindowShift, "Get the shift in the window since last step.")
	DeclFunc("ext_activitydwpos", GetActivityDWPos, "Get the domain wall position from the activityStack")
}

// DWActivityInit(w) sets the mask width to apply to the domain wall; only values of the magnetization within w cells
// of the domain wall are included in the domain wall activity
func DWActivityInit(w int) {
	DWMonitor.maskWidth = w
	DWMonitor.initialized = false
	return
}

type activityStack struct {
	windowpos   [2]int
	dwpos       [2][][]int
	rxy         [2][][][]float64
	phi         [2][][][]float64
	theta       [2][][][]float64
	t           [2]float64
	maskWidth   int
	_AzCache    float64
	_AxyCache   float64
	initialized bool
}

func (s *activityStack) Az() float64 {
	if s.t[1] != Time {
		s.push()
	}
	return s._AzCache
}

func (s *activityStack) Axy() float64 {
	if s.t[1] != Time {
		s.push()
	}
	return s._AxyCache
}

func getAz() float64 {
	return DWMonitor.Az()
}

func getAxy() float64 {
	return DWMonitor.Axy()
}

func (s *activityStack) lastTime() float64 {
	return s.t[1]
}

func (s *activityStack) push() {

	// Move the values of the angles, dw positions, time, and window shift to the old slots
	s.rxy[0] = s.rxy[1]
	s.phi[0] = s.phi[1]
	s.theta[0] = s.theta[1]
	s.t[0] = s.t[1]
	s.dwpos[0] = s.dwpos[1]
	s.windowpos[0] = s.windowpos[1]

	// Update the DW and window positions; put new values of the angles into the new slots.
	if s.initialized {
		s.windowpos[1] = GetIntWindowPos()
		s.dwpos[1] = GetIntDWPos()

		// Get r, phi, theta corresponding to the previous domain wall configuration. These in principle will have
		// evolved in time since the last time this function was called, so while these should be the same spins as
		// were stored in the last time step, they won't be pointing in the same direction.
		s.rxy[1], s.phi[1], s.theta[1] = rxyPhiThetaBand(s.dwpos[0], s.maskWidth, s.WindowShift())
		s.t[1] = Time

		s._AzCache = calcAz(s.theta[1], s.theta[0], s.t[1]-s.t[0])
		s._AxyCache = calcAxy(s.rxy[1], s.rxy[0], s.phi[1], s.phi[0], s.t[1]-s.t[0])

		// Get r, phi, theta from the new domain wall configuration (therefore no window correction necessary).
		s.rxy[0], s.phi[0], s.theta[0] = rxyPhiThetaBand(s.dwpos[1], s.maskWidth, 0)

	} else {
		s.windowpos[1] = GetIntWindowPos()
		s.dwpos[1] = GetIntDWPos()
		s.rxy[1], s.phi[1], s.theta[1] = rxyPhiThetaBand(s.dwpos[1], s.maskWidth, 0)
		s.t[1] = Time
		s.initialized = true
	}
	return
}

// Calculate the change in angle for two angles, taking into account the fact that 2*pi = 0. If the difference in angles
// (a-b) is large, the vector is assumed to have wrapped around this boundary.
func deltaAngle(a, b float64) float64 {
	dA := a - b
	if dA < -math.Pi {
		return 2*math.Pi + dA
	} else if dA > math.Pi {
		return 2*math.Pi - dA
	}
	return dA
}

func calcAz(thetaNew, thetaOld [][][]float64, dt float64) float64 {

	n := MeshSize()
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < len(thetaNew[k][j]); i++ {
				ret += deltaAngle(thetaNew[k][j][i], thetaOld[k][j][i]) / dt
			}
		}
	}
	return ret
}

func calcAxy(rxyNew, rxyOld, phiNew, phiOld [][][]float64, dt float64) float64 {

	n := MeshSize()
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < len(phiNew[k][j]); i++ {
				ret += deltaAngle(phiNew[k][j][i], phiOld[k][j][i]) * (rxyNew[k][j][i] + rxyOld[k][j][i]) * 0.5 / dt
			}
		}
	}
	return ret
}

func rxyPhiThetaBand(_dwpos [][]int, width int, windowShift int) ([][][]float64, [][][]float64, [][][]float64) {

	n := MeshSize()
	m := M.Buffer().HostCopy().Vectors()

	_rxy := make([][][]float64, n[Z])
	_phi := make([][][]float64, n[Z])
	_theta := make([][][]float64, n[Z])

	for k := 0; k < n[Z]; k++ {
		_rxy[k] = make([][]float64, n[Y])
		_phi[k] = make([][]float64, n[Y])
		_theta[k] = make([][]float64, n[Y])

		for j := 0; j < n[Y]; j++ {
			_rxy[k][j] = make([]float64, 2*width+1)
			_phi[k][j] = make([]float64, 2*width+1)
			_theta[k][j] = make([]float64, 2*width+1)

			iMin := _dwpos[k][j] - width - windowShift
			iMax := _dwpos[k][j] + 1 + width - windowShift

			im := 0
			for i := iMin; i < iMax; i++ {
				_rxy[k][j][im] = rxy(float64(m[X][k][j][i]), float64(m[Y][k][j][i]))
				_phi[k][j][im] = phi(float64(m[X][k][j][i]), float64(m[Y][k][j][i]))
				_theta[k][j][im] = theta(float64(m[Z][k][j][i]))
				im++
			}
		}
	}
	return _rxy, _phi, _theta
}

func rxy(mx float64, my float64) float64 {
	return math.Sqrt(float64(mx*mx + my*my))
}

func phi(mx float64, my float64) float64 {
	return math.Atan2(float64(my), float64(mx))
}

func theta(mz float64) float64 {
	return math.Acos(float64(mz))
}

func getPhi() [][][]float64 {
	return DWMonitor.phi[1]
}

func getTheta() [][][]float64 {
	return DWMonitor.theta[1]
}

func IntRound(x float64) int {
	return int(x + 0.5)
}

// Get the indices of the domain wall within the simulation window
func GetIntDWPos() [][]int {

	n := MeshSize()
	mz := M.Buffer().Comp(Z).HostCopy().Scalars()
	pos := make([][]int, n[Z])

	for i := 0; i < len(mz); i++ {
		pos[i] = make([]int, n[Y])
		for j := 0; j < n[Y]; j++ {
			pos[i][j] = _getIntDWPos1D(mz[i][j])
		}
	}
	return pos
}

func _getIntDWPos1D(mz []float32) int {
	ix := ZeroCrossing(mz)
	if math.Abs(float64(mz[ix])) < math.Abs(float64(mz[ix+1])) {
		return ix
	}
	return ix + 1
}

func (s *activityStack) DWShift() [][]int {
	n := MeshSize()
	ret := make([][]int, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([]int, n[Y])
		for j := 0; j < n[Y]; j++ {
			ret[i][j] = s.dwpos[1][i][j] - s.dwpos[0][i][j]
		}
	}
	return ret
}

func (s *activityStack) WindowShift() int {
	return s.windowpos[1] - s.windowpos[0]
}

// Number of cells the window has shifted
func GetIntWindowPos() int {
	c := Mesh().CellSize()
	windowPos := GetShiftPos()
	return IntRound(windowPos / c[Y])
}

// Generate a zeroed slice for initializing the DW shift
func _zeroDWShift() [][]int {
	n := MeshSize()
	ret := make([][]int, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([]int, n[Y])
	}
	return ret
}

func GetWindowShift() int {
	return DWMonitor.WindowShift()
}

func GetActivityDWPos() [][][]float64 {

	n := MeshSize()
	pos := DWMonitor.dwpos[1]
	width := DWMonitor.maskWidth

	ret := make([][][]float64, n[Z])
	for i:=0; i<n[Z]; i++ {
		ret[i] = make([][]float64, n[Y])
		for j:=0; j<n[Y]; j++ {
			ret[i][j] = make([]float64, n[X])
			for k:=0; k<n[X]; k++ {
				if k > pos[i][j]-width && k < pos[i][j]+width+1 {
					ret[i][j][k] = 1.0
				} else {
					ret[i][j][k] = 0.0
				}
			}
		}
	}
	return ret
}