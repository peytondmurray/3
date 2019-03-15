package engine

import (
	// "log"
	// "fmt"
	// "github.com/mumax/3/data"
	"math"
)

var (
	Az         = NewScalarValue("ext_az", "rad/s", "Out-of-plane domain wall activity", getAz)
	Axy        = NewScalarValue("ext_axy", "rad/s", "In-plane domain wall activity", getAxy)
	exactDWVelAvg = NewScalarValue("ext_exactdwvelavg", "m/s", "Domain wall velocity, found via average value of mz", getExactDWVelAvg)
	exactDWPosAvg = NewScalarValue("ext_exactdwposavg", "m", "Domain wall position, found via average value of mz", getExactDWPosAvg)
	exactDWVelZC = NewScalarValue("ext_exactdwvelzc", "m/s", "Domain wall velocity, found via zero crossing of mz", getExactDWVelZC)
	exactDWPosZC = NewScalarValue("ext_exactdwposzc", "m", "Domain wall position, found via zero crossing of mz", getExactDWPosZC)
	DWMonitor  activityStack // Most recent positions of DW speed
)

func init() {
	DeclFunc("ext_dwactivityinit", DWActivityInit, "ext_dwactivityinit(w, l, r, kind) sets the mask width to w, the sign of M which is being inserted at the left l and right r sides of the simulation, and the kind of DW position calculation ('zc' or 'avg').")
	DeclFunc("ext_getphi", getPhi, "Get the current phi angle as a slice.")
	DeclFunc("ext_gettheta", getTheta, "Get the current theta angle as a slice.")
	DeclFunc("ext_getphidot", getPhiDot, "Get the current phi angle as a slice.")
	DeclFunc("ext_getthetadot", getThetaDot, "Get the current theta angle as a slice.")
}

// DWActivityInit(w) sets the mask width to apply to the domain wall; only values of the magnetization within w cells
// of the domain wall are included in the domain wall activity
func DWActivityInit(w int, l int, r int, kind string) {
	DWMonitor.maskWidth = w
	DWMonitor.signL = l
	DWMonitor.signR = r
	DWMonitor.initialized = false
	DWMonitor.postype = kind
	return
}

type dwAvgPosStack struct {
	t float64
	dwexactvel float64
	dwexactpos float64
}

func(s *dwAvgPosStack) init(t float64, mz [][][]float32) {
	s.dwexactpos = exactPosAvg(mz)
	s.dwexactvel = 0
	return
}

func (s *dwAvgPosStack) update(t float64, mz [][][]float32) {
	dwexactpos := exactPosAvg(mz)
	dwexactvel := exactVel(dwexactpos, s.dwexactpos, t, s.t)
	s.dwexactpos = dwexactpos
	s.dwexactvel = dwexactvel
	return
}

type dwZCPosStack struct {
	t float64
	dwexactvel float64
	dwexactpos float64
}

func(s *dwZCPosStack) init(t float64, mz [][][]float32, dwPos [][]int) {
	s.dwexactpos = exactPosZC(mz, dwPos)
	s.dwexactvel = 0
	return
}

func (s *dwZCPosStack) update(t float64, mz [][][]float32, dwPos [][]int) {
	dwexactpos := exactPosZC(mz, dwPos)
	dwexactvel := exactVel(dwexactpos, s.dwexactpos, t, s.t)
	s.dwexactpos = dwexactpos
	s.dwexactvel = dwexactvel
	return
}

type activityStack struct {
	signL       int
	signR       int
	postype     string
	windowpos   int
	dwexactMonitorZC dwZCPosStack
	dwexactMonitorAvg dwAvgPosStack
	dwpos       [][]int
	rxy         [][][]float64
	phi         [][][]float64
	theta       [][][]float64
	phidot      [][][]float64
	thetadot    [][][]float64
	t           float64
	maskWidth   int
	Az          float64
	Axy         float64
	initialized bool
}

func (s *activityStack) getExactDWPosZC() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.dwexactMonitorZC.dwexactpos
}

func (s *activityStack) getExactDWVelZC() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.dwexactMonitorZC.dwexactvel
}

func (s *activityStack) getExactDWPosAvg() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.dwexactMonitorAvg.dwexactpos
}

func (s *activityStack) getExactDWVelAvg() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.dwexactMonitorAvg.dwexactvel
}

func (s *activityStack) getAz() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.Az
}

func (s *activityStack) getAxy() float64 {
	if s.t != Time || !s.initialized {
		s.push()
	}
	return s.Axy
}

func getExactDWPosZC() float64 {
	return DWMonitor.getExactDWPosZC()
}

func getExactDWVelZC() float64 {
	return DWMonitor.getExactDWVelZC()
}
func getExactDWPosAvg() float64 {
	return DWMonitor.getExactDWPosAvg()
}

func getExactDWVelAvg() float64 {
	return DWMonitor.getExactDWVelAvg()
}

func getAz() float64 {
	return DWMonitor.getAz()
}

func getAxy() float64 {
	return DWMonitor.getAxy()
}

func (s *activityStack) lastTime() float64 {
	return s.t
}

func (s *activityStack) push() {

	// Update the DW and window positions; put new values of the angles into the new slots.
	if s.initialized {

		// Get the new magnetization configuration
		_m := M.Buffer().HostCopy().Vectors()

		// Get new window and domain wall positions. Get the newest rxy, phi, theta values.
		_windowpos := GetIntWindowPos()
		_dwpos := GetIntDWPos(_m[Z])
		_rxy, _phi, _theta := rxyPhiTheta(_m)
		_t := Time

		if s.postype == "zc" {
			s.dwexactMonitorZC.update(_t, _m[Z], _dwpos)
		} else if s.postype == "avg" {
			s.dwexactMonitorAvg.update(_t, _m[Z])
		} else if s.postype == "both" {
			s.dwexactMonitorZC.update(_t, _m[Z], _dwpos)
			s.dwexactMonitorAvg.update(_t, _m[Z])
		} else {
			panic("Invalid DW position calculation type; must be 'zc' or 'avg' or 'both'.")
		}
		_rxyAvg := averageRxy(_rxy, s.rxy)
		_phidot := angularVel(_phi, s.phi, _windowpos, s.windowpos, _t, s.t)
		_thetadot := angularVel(_theta, s.theta, _windowpos, s.windowpos, _t, s.t)

		// Calculate Axy and Az by summing the angular velocity of the cells which were near the domain wall at the
		// last step (i.e., within maskWidth cells of the old dwpos).
		s.Axy = calcAxy(_phidot, _rxyAvg, s.dwpos, s.maskWidth)
		s.Az = calcAz(_thetadot, s.dwpos, s.maskWidth)

		s.windowpos = _windowpos
		s.dwpos = _dwpos
		s.rxy = _rxy
		s.phi = _phi
		s.theta = _theta
		s.t = _t
		s.phidot = _phidot
		s.thetadot = _thetadot

	} else {

		_m := M.Buffer().HostCopy().Vectors()

		s.t = Time
		s.windowpos = GetIntWindowPos()
		s.dwpos = GetIntDWPos(_m[Z])

		if s.postype == "zc" {
			s.dwexactMonitorZC.init(s.t, _m[Z], s.dwpos)
		} else if s.postype == "avg" {
			s.dwexactMonitorAvg.init(s.t, _m[Z])
		} else if s.postype == "both" {
			s.dwexactMonitorZC.init(s.t, _m[Z], s.dwpos)
			s.dwexactMonitorAvg.init(s.t, _m[Z])
		} else {
			panic("Invalid DW position calculation type; must be 'zc' or 'avg' or 'both'.")
		}

		s.rxy, s.phi, s.theta = rxyPhiTheta(_m)

		s.phidot = ZeroWorld()
		s.thetadot = ZeroWorld()
		s.initialized = true
	}
	return
}

func exactPosAvg(mz [][][]float32) float64 {

	ws := Mesh().WorldSize()

	// Get average magnetization
	avg := avgMz(mz)

	// Percentage of the magnetization which is flipped up gives the position of the domain wall, for example if
	// 50% are flipped up, the DW is 50% from the left side of the simulation window
	pct := 1.0 - (1.0-avg)/2.0

	// Convert to actual position in window, then add on window shift
	return pct*ws[X] + GetShiftPos()
}

func avgMz(mz [][][]float32) float64 {
	n := MeshSize()
	sum := float32(0.0)
	for i := range mz {
		for j := range mz[i] {
			for k := range mz[i][j] {
				sum += mz[i][j][k]
			}
		}
	}
	return float64(sum) / float64(n[X]*n[Y]*n[Z])
}

// Get the exact position of the domain wall. Takes into account the shift in the simulation window, as well as the
// position of the domain wall within the window.
func exactPosZC(mz [][][]float32, intPos [][]int) float64 {
	c := Mesh().CellSize()
	return exactPosInWindow(mz, intPos)*c[X] + GetShiftPos()
}

// Find the exact average DW position in the simulation space in units of Mesh().CellSize()[X]
func exactPosInWindow(mz [][][]float32, intPos [][]int) float64 {
	pos := float64(0)
	for iz := range mz {
		for iy := range mz[iz] {
			pos += float64(_interpolateZeroCrossing(mz[iz][iy], intPos[iz][iy]))
		}
	}
	return pos / float64(len(intPos)*len(intPos[0]))
}

// Get the exact domain wall velocity. Takes into account the shift in the simulation window.
func exactVel(posNew, posOld, tNew, tOld float64) float64 {
	return (posNew - posOld) / (tNew - tOld)
}

// Interpolate the index of the 1D slice wehre the zero crossing of the Mz component occurs, between x-index i and i+1.
// Returns a float32, so you can't use this for indexing an array at the zero crossing point.
func _interpolateZeroCrossing(mz []float32, i int) float32 {
	return float32(i) - (mz[i] / (mz[i+1] - mz[i]))
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

func calcAxy(phiDot, rxyAvg [][][]float64, dwpos [][]int, maskWidth int) float64 {

	n := MeshSize()
	ret := float64(0)
	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := dwpos[i][j] - maskWidth; k < dwpos[i][j]+maskWidth+1; k++ {
				ret += phiDot[i][j][k] * rxyAvg[i][j][k]
			}
		}
	}
	return ret
}

func calcAz(thetaDot [][][]float64, dwpos [][]int, maskWidth int) float64 {

	n := MeshSize()
	ret := float64(0)
	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := dwpos[i][j] - maskWidth; k < dwpos[i][j]+maskWidth+1; k++ {
				ret += thetaDot[i][j][k]
			}
		}
	}
	return ret
}

func rxyPhiTheta(m [3][][][]float32) ([][][]float64, [][][]float64, [][][]float64) {

	n := MeshSize()

	_rxy := make([][][]float64, n[Z])
	_phi := make([][][]float64, n[Z])
	_theta := make([][][]float64, n[Z])

	for k := 0; k < n[Z]; k++ {
		_rxy[k] = make([][]float64, n[Y])
		_phi[k] = make([][]float64, n[Y])
		_theta[k] = make([][]float64, n[Y])

		for j := 0; j < n[Y]; j++ {
			_rxy[k][j] = make([]float64, n[X])
			_phi[k][j] = make([]float64, n[X])
			_theta[k][j] = make([]float64, n[X])

			for i := 0; i < n[X]; i++ {
				_rxy[k][j][i] = rxy(float64(m[X][k][j][i]), float64(m[Y][k][j][i]))
				_phi[k][j][i] = phi(float64(m[X][k][j][i]), float64(m[Y][k][j][i]))
				_theta[k][j][i] = theta(float64(m[Z][k][j][i]))
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
	return DWMonitor.phi
}

func getTheta() [][][]float64 {
	return DWMonitor.theta
}

func IntRound(x float64) int {
	return int(x + 0.5)
}

// Get the indices of the domain wall within the simulation window.
func GetIntDWPos(mz [][][]float32) [][]int {

	n := MeshSize()
	pos := make([][]int, n[Z])

	for i := 0; i < len(mz); i++ {
		pos[i] = make([]int, n[Y])
		for j := 0; j < n[Y]; j++ {
			pos[i][j] = _getIntDWPos1D(mz[i][j])
		}
	}
	return pos
}

// Check from right to left, since the domain wall moves in that direction; this means that bubbles will only be
// forming to the left, and we will not be susceptible to them with this method.
func _getIntDWPos1D(mz []float32) int {
	ix := zeroCrossing(mz)
	if math.Abs(float64(mz[ix-1])) < math.Abs(float64(mz[ix])) {
		return ix - 1
	}
	return ix
}

// Find the index of the 1D slice where the zero crossing of the Mz component occurs. Check from right to left.
func zeroCrossing(mz []float32) int {
	for ix := len(mz) - 1; ix > 1; ix-- {
		if sign32(mz[ix-1]) == DWMonitor.signL && sign32(mz[ix]) == DWMonitor.signR {
			return ix - 1
		}
	}
	panic("Can't find domain wall position")
}

func sign32(a float32) int {
	if a < 0 {
		return -1
	}
	return 1
}

// Number of cells the window has shifted
func GetIntWindowPos() int {
	c := Mesh().CellSize()
	windowPos := GetShiftPos()
	return IntRound(windowPos / c[Y])
}

func angularVel(aNew, aOld [][][]float64, windowposNew, windowposOld int, tNew, tOld float64) [][][]float64 {

	shift := windowposNew - windowposOld
	dt := tNew - tOld
	n := MeshSize()

	ret := make([][][]float64, n[Z])

	if shift >= 0 {
		for i := 0; i < n[Z]; i++ {
			ret[i] = make([][]float64, n[Y])
			for j := 0; j < n[Y]; j++ {
				ret[i][j] = make([]float64, n[X])
				for k := 0; k < n[X]-shift; k++ {
					ret[i][j][k] = deltaAngle(aNew[i][j][k], aOld[i][j][k+shift]) / dt
				}
			}
		}
	} else {
		for i := 0; i < n[Z]; i++ {
			ret[i] = make([][]float64, n[Y])
			for j := 0; j < n[Y]; j++ {
				ret[i][j] = make([]float64, n[X])
				for k := -shift; k < n[X]; k++ {
					ret[i][j][k] = deltaAngle(aNew[i][j][k], aOld[i][j][k+shift]) / dt
				}
			}
		}
	}
	return ret
}

// Make a slice the same size as the simulation, initialized with zeros.
func ZeroWorld() [][][]float64 {

	n := MeshSize()

	ret := make([][][]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([][]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			ret[i][j] = make([]float64, n[X])
			for k := 0; k < n[X]; k++ {
				ret[i][j][k] = 0
			}
		}
	}
	return ret
}

func getPhiDot() [][][]float64 {
	return DWMonitor.phidot
}

func getThetaDot() [][][]float64 {
	return DWMonitor.thetadot
}

// Find the average of old and new rxy slices.
func averageRxy(rxyNew, rxyOld [][][]float64) [][][]float64 {
	n := MeshSize()

	ret := make([][][]float64, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([][]float64, n[Y])
		for j := 0; j < n[Y]; j++ {
			ret[i][j] = make([]float64, n[X])
			for k := 0; k < n[X]; k++ {
				ret[i][j][k] = 0.5 * (rxyNew[i][j][k] + rxyOld[i][j][k])
			}
		}
	}
	return ret
}
