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
	DeclFunc("ext_dphidt", getDPhiDt, "Find dphi/dt")
	DeclFunc("ext_dthetadt", getDThetaDt, "Find dtheta/dt")
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
	s.dwpos[1] = GetIntDWPos()
	s.windowpos[1] = GetIntWindowPos()
	_m := magnetizationInBand(s.dwpos[1], s.maskWidth, s.WindowShift(), s.DWShift())
	s.rxy[1], s.phi[1], s.theta[1] = rxyPhiTheta(_m)
	s.t[1] = Time

	if s.initialized {
		s._AzCache = calcAz(s.theta[1], s.theta[0], s.t[1]-s.t[0])
		s._AxyCache = calcAxy(s.rxy[1], s.phi[1], s.phi[0], s.t[1]-s.t[0])
	} else {
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
			for i := 0; i <= len(thetaNew); i++ {
				ret += deltaAngle(thetaNew[k][j][i], thetaOld[k][j][i]) / dt
			}
		}
	}
	return ret
}

func calcAxy(rxyNew, phiNew, phiOld [][][]float64, dt float64) float64 {

	n := MeshSize()
	ret := float64(0)

	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < len(phiNew); i++ {
				ret += deltaAngle(phiNew[k][j][i], phiOld[k][j][i]) * rxyNew[k][j][i] / dt
			}
		}
	}
	return ret
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

// Return the magnetization within the region near the DW in (x, y, z, component) order.
func magnetizationInBand(dwpos [][]int, width int, windowShift int, dwShift [][]int) [][][][3]float64 {
	n := MeshSize()

	m := make([][][][3]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		m[i] = make([][][3]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			m[i][j] = make([][3]float64, 2*width+1)

			iMin := dwpos[i][j] - width - windowShift - dwShift[i][j]
			iMax := dwpos[i][j] + 1 + width - windowShift - dwShift[i][j]

			for k := iMin; k < iMax; k++ {
				_m := M.GetCell(k, j, i)
				m[i][j][k] = [3]float64{_m[X], _m[Y], _m[Z]}
			}
		}
	}
	return m
}

// For a given magnetization slice in (x, y, z, component) order, return the in-plane magnitude, phi angle, and theta
// angle.
func rxyPhiTheta(m [][][][3]float64) ([][][]float64, [][][]float64, [][][]float64) {

	_rxy := make([][][]float64, len(m))
	_phi := make([][][]float64, len(m))
	_theta := make([][][]float64, len(m))

	for k := 0; k < len(m); k++ {
		_rxy[k] = make([][]float64, len(m[k]))
		_phi[k] = make([][]float64, len(m[k]))
		_theta[k] = make([][]float64, len(m[k]))

		for j := 0; j < len(m[k]); j++ {
			_rxy[k][j] = make([]float64, len(m[k][j]))
			_phi[k][j] = make([]float64, len(m[k][j]))
			_theta[k][j] = make([]float64, len(m[k][j]))

			for i := 0; i < len(m[k][j]); i++ {
				_rxy[k][j][i] = rxy(m[k][j][i][X], m[k][j][i][Y])
				_phi[k][j][i] = phi(m[k][j][i][X], m[k][j][i][Y])
				_theta[k][j][i] = theta(m[k][j][i][Z])
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

func (s *activityStack) dphidt() [][][]float64 {
	n := MeshSize()

	diff := make([][][]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		diff[i] = make([][]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			diff[i][j] = make([]float64, n[X])
			for k := 0; k < n[X]; k++ {
				diff[i][j][k] = deltaAngle(s.phi[1][i][j][k], s.phi[0][i][j][k]) / (s.t[1] - s.t[0])
			}
		}
	}
	return diff
}

func (s *activityStack) dthetadt() [][][]float64 {
	n := MeshSize()

	diff := make([][][]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		diff[i] = make([][]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			diff[i][j] = make([]float64, n[X])
			for k := 0; k < n[X]; k++ {
				diff[i][j][k] = deltaAngle(s.theta[1][i][j][k], s.theta[0][i][j][k]) / (s.t[1] - s.t[0])
			}
		}
	}
	return diff
}

func getDPhiDt() [][][]float64 {
	if !DWMonitor.initialized {
		DWMonitor.push()
		return zeroedGlobalSlice()
	} else if DWMonitor.lastTime() != Time {
		DWMonitor.push()
	}
	return DWMonitor.dphidt()
}

func getDThetaDt() [][][]float64 {
	if !DWMonitor.initialized {
		DWMonitor.push()
		return zeroedGlobalSlice()
	} else if DWMonitor.lastTime() != Time {
		DWMonitor.push()
	}
	return DWMonitor.dthetadt()
}

func zeroedGlobalSlice() [][][]float64 {
	n := MeshSize()

	zero := make([][][]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		zero[i] = make([][]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			zero[i][j] = make([]float64, n[X])
		}
	}
	return zero
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
