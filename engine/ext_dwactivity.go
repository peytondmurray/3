package engine

import (
	"github.com/mumax/3/data"
	// "github.com/mumax/3/mag"
	// "log"
	// "fmt"
	"github.com/mumax/3/cuda"
	"math"
)

var (
	Az            = NewScalarValue("ext_az", "rad/s", "Out-of-plane domain wall activity", getAz)
	Axy           = NewScalarValue("ext_axy", "rad/s", "In-plane domain wall activity", getAxy)
	ExactDWVelAvg = NewScalarValue("ext_exactdwvelavg", "m/s", "Speed of domain wall", getExactVelAvg)
	ExactDWPosAvg = NewScalarValue("ext_exactdwposavg", "m", "Position of domain wall from start", getExactPosAvg)
	ExactDWPosZC  = NewScalarValue("ext_exactdwposzc", "m", "Position of the domain wall from start", getExactPosZC)
	DWMonitor     activityStack // Most recent positions of DW speed
)

func init() {
	DeclFunc("ext_dwactivityinit", DWActivityInit, "ext_dwactivityinit(w, l, r) sets the mask width to w, the sign of M which is being inserted at the left l and right r sides of the simulation.")
	DeclFunc("ext_getphi", getPhi, "Get the current phi angle as a slice.")
	DeclFunc("ext_gettheta", getTheta, "Get the current theta angle as a slice.")
	DeclFunc("ext_getphidot", getPhiDot, "Get the current phi angle as a slice.")
	DeclFunc("ext_getthetadot", getThetaDot, "Get the current theta angle as a slice.")
	DeclFunc("ext_debug_setdw", debugSetDWMonitor, "Set DW parameters")
}

func avgDiff() float64 {
	averageCPU := float64(avgMz(M.Comp(Z).HostCopy().Scalars()))
	averageGPU := M.Comp(Z).Average()
	return averageCPU - averageGPU
}

// DWActivityInit(w) sets the mask width to apply to the domain wall; only values of the magnetization within w cells
// of the domain wall are included in the domain wall activity
func DWActivityInit(w int, l int, r int) {
	DWMonitor.maskWidth = w
	DWMonitor.signL = l
	DWMonitor.signR = r
	DWMonitor.initialized = false
	return
}

type activityStack struct {

	// DWVel
	posAvg float64
	velAvg float64

	// DWActivity
	signL       int
	signR       int
	windowpos   int
	dwpos       [][]int
	rxy         [][][]float32
	phi         [][][]float32
	theta       [][][]float32
	phidot      [][][]float32
	thetadot    [][][]float32
	t           float64
	lastStep    int
	maskWidth   int
	Az          float32
	Axy         float32
	initialized bool

	// Pointers to rxy, phi, theta, and the angular velocities
	p_rpt *data.Slice
	p_pd  *data.Slice
	p_td  *data.Slice
}

func (s *activityStack) update() {
	if !s.initialized {
		s.init()
	} else if s.t != Time || s.lastStep != NSteps {
		s.push()
	}
	return
}

func debugSetDWMonitor(vel float64) {

	_rpt := ext_rxyphitheta.HostCopy().Vectors()

	DWMonitor.t = Time
	DWMonitor.lastStep = NSteps
	DWMonitor.windowpos = GetIntWindowPos()

	_intPosZC := getIntDWPos(_rpt[2]) // Move up above GetNearestIntDWPos
	DWMonitor.posAvg = exactPosAvg()
	DWMonitor.velAvg = vel
	DWMonitor.dwpos = GetNearestIntDWPos(_rpt[2], _intPosZC)

	DWMonitor.rxy = _rpt[0]
	DWMonitor.phi = _rpt[1]
	DWMonitor.theta = _rpt[2]

	DWMonitor.phidot = ZeroWorldScalar32()
	DWMonitor.thetadot = ZeroWorldScalar32()
	DWMonitor.initialized = true

	return

}

func getExactPosZC() float64 {
	_rpt := ext_rxyphitheta.HostCopy().Vectors()
	_mz := M.Comp(Z).HostCopy().Scalars()
	_intPosZC := getIntDWPos(_rpt[2])

	return float64(exactPosZC(_mz, _intPosZC))
}

func getExactPosAvg() float64 {
	DWMonitor.update()
	return DWMonitor.posAvg
}

func getExactVelAvg() float64 {
	DWMonitor.update()
	return DWMonitor.velAvg
}

func getAz() float64 {
	DWMonitor.update()
	return float64(DWMonitor.Az)
}

func getAxy() float64 {
	DWMonitor.update()
	return float64(DWMonitor.Axy)
}

func (s *activityStack) lastTime() float64 {
	return s.t
}

func (s *activityStack) init() {

	_rpt := ext_rxyphitheta.HostCopy().Vectors()

	s.t = Time
	s.lastStep = NSteps
	s.windowpos = GetIntWindowPos()

	_intPosZC := getIntDWPos(_rpt[2]) // Move up above GetNearestIntDWPos
	s.posAvg = exactPosAvg()
	s.velAvg = 0.0
	s.dwpos = GetNearestIntDWPos(_rpt[2], _intPosZC)

	s.rxy = _rpt[0]
	s.phi = _rpt[1]
	s.theta = _rpt[2]

	s.phidot = ZeroWorldScalar32()
	s.thetadot = ZeroWorldScalar32()
	s.initialized = true

	// Set aside buffers for holding r, phi, and theta
	s.p_rpt = cuda.Buffer(3, MeshSize()) // Possibly remove?
	s.p_pd = cuda.Buffer(1, MeshSize())
	s.p_td = cuda.Buffer(1, MeshSize())

	return
}

func (s *activityStack) push() {

	// Update the DW and window positions; put new values of the angles into the new slots.
	// Get the new magnetization configuration
	_rpt := ext_rxyphitheta.HostCopy().Vectors()

	_t := Time

	// DWVel_________________________
	_intPosZC := getIntDWPos(_rpt[2])
	_posAvg := exactPosAvg()
	s.velAvg = (_posAvg - s.posAvg) / (_t - s.t)
	s.posAvg = _posAvg
	// DWVel_________________________

	// Get new window and domain wall positions. Get the newest rxy, phi, theta values.
	_windowpos := GetIntWindowPos()
	_dwpos := GetNearestIntDWPos(_rpt[2], _intPosZC)

	_rxyAvg := averageRxy(_rpt[0], s.rxy)
	_phidot := angularVel(_rpt[1], s.phi, _windowpos, s.windowpos, _t, s.t)
	_thetadot := angularVel(_rpt[2], s.theta, _windowpos, s.windowpos, _t, s.t)

	// Calculate Axy and Az by summing the angular velocity of the cells which were near the domain wall at the
	// last step (i.e., within maskWidth cells of the old dwpos).
	s.Axy = calcAxy(_phidot, _rxyAvg, s.dwpos, s.maskWidth)
	s.Az = calcAz(_thetadot, s.dwpos, s.maskWidth)

	s.windowpos = _windowpos
	s.dwpos = _dwpos
	s.rxy = _rpt[0]
	s.phi = _rpt[1]
	s.theta = _rpt[2]
	s.t = _t
	s.lastStep = NSteps
	s.phidot = _phidot
	s.thetadot = _thetadot

	return
}

// Calculate the change in angle for two angles, taking into account the fact that 2*pi = 0. If the difference in angles
// (a-b) is large, the vector is assumed to have wrapped around this boundary.
func deltaAngle(a, b float32) float32 {
	dA := a - b
	if dA < -math.Pi {
		return 2*math.Pi + dA
	} else if dA > math.Pi {
		return 2*math.Pi - dA
	}
	return dA
}

func calcAxy(phiDot, rxyAvg [][][]float32, dwpos [][]int, maskWidth int) float32 {

	n := MeshSize()
	ret := float32(0)
	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := dwpos[i][j] - maskWidth; k < dwpos[i][j]+maskWidth+1; k++ {
				ret += phiDot[i][j][k] * rxyAvg[i][j][k]
			}
		}
	}
	return ret
}

func calcAz(thetaDot [][][]float32, dwpos [][]int, maskWidth int) float32 {

	n := MeshSize()
	ret := float32(0)
	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := dwpos[i][j] - maskWidth; k < dwpos[i][j]+maskWidth+1; k++ {
				ret += thetaDot[i][j][k]
			}
		}
	}
	return ret
}

func rxyPhiTheta(m [3][][][]float32) ([][][]float32, [][][]float32, [][][]float32) {
	_rxyphitheta := ext_rxyphitheta.HostCopy().Vectors()
	return _rxyphitheta[X], _rxyphitheta[Y], _rxyphitheta[Z]
}

func rxy(mx float32, my float32) float32 {
	return float32(math.Sqrt(float64(mx*mx + my*my)))
}

func phi(mx float32, my float32) float32 {
	return float32(math.Atan2(float64(my), float64(mx)))
}

func theta(mz float32) float32 {
	return float32(math.Acos(float64(mz)))
}

func getPhi() [][][]float32 {
	return DWMonitor.phi
}

func getTheta() [][][]float32 {
	return DWMonitor.theta
}

func IntRound(x float32) int {
	return int(x + 0.5)
}

func GetNearestIntDWPos(theta [][][]float32, dw [][]int) [][]int {

	nearest := make([][]int, len(dw))
	for i := range theta {
		nearest[i] = make([]int, len(dw[i]))
		for j := range theta[i] {
			if sign64(math.Cos(float64(theta[i][j][dw[i][j]-1]))) == DWMonitor.signL && sign64(math.Cos(float64(theta[i][j][dw[i][j]]))) == DWMonitor.signR {
				nearest[i][j] = dw[i][j] - 1
			} else {
				nearest[i][j] = dw[i][j]
			}
		}
	}
	return nearest
}

// Find the index of the 1D slice just before the zero crossing of the Mz component. Check from right to left;
// since bubbles will only be forming to the left, and we will not be susceptible to them with this method.
func zeroCrossing(theta []float32) int {
	for ix := len(theta) - 1; ix > 1; ix-- {
		if sign64(math.Cos(float64(theta[ix-1]))) == DWMonitor.signL && sign64(math.Cos(float64(theta[ix]))) == DWMonitor.signR {
			return ix - 1
		}
	}
	panic("Can't find domain wall position")
}

func sign64(a float64) int {
	if a < 0 {
		return -1
	}
	return 1
}

// Number of cells the window has shifted
func GetIntWindowPos() int {
	c := Mesh().CellSize()
	windowPos := GetShiftPos()
	return IntRound(float32(windowPos) / float32(c[Y]))
}

func angularVel(aNew, aOld [][][]float32, windowposNew, windowposOld int, tNew, tOld float64) [][][]float32 {

	shift := windowposNew - windowposOld
	dt := float32(tNew - tOld)
	n := MeshSize()

	ret := make([][][]float32, n[Z])

	if shift >= 0 {
		for i := 0; i < n[Z]; i++ {
			ret[i] = make([][]float32, n[Y])
			for j := 0; j < n[Y]; j++ {
				ret[i][j] = make([]float32, n[X])
				for k := 0; k < n[X]-shift; k++ {
					ret[i][j][k] = deltaAngle(aNew[i][j][k], aOld[i][j][k+shift]) / dt
				}
			}
		}
	} else {
		for i := 0; i < n[Z]; i++ {
			ret[i] = make([][]float32, n[Y])
			for j := 0; j < n[Y]; j++ {
				ret[i][j] = make([]float32, n[X])
				for k := -shift; k < n[X]; k++ {
					ret[i][j][k] = deltaAngle(aNew[i][j][k], aOld[i][j][k+shift]) / dt
				}
			}
		}
	}
	return ret
}

// Make a slice the same size as the simulation, initialized with zeros.
func ZeroWorldScalar32() [][][]float32 {

	n := MeshSize()

	ret := make([][][]float32, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([][]float32, n[Y])
		for j := 0; j < n[Y]; j++ {
			ret[i][j] = make([]float32, n[X])
			for k := 0; k < n[X]; k++ {
				ret[i][j][k] = 0
			}
		}
	}
	return ret
}

// Make a slice the same size as the simulation, initialized with zeros.
func ZeroWorldScalar64() [][][]float64 {

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

// Make a slice the same size as the simulation, initialized with zeros.
func ZeroWorldVector64() [3][][][]float64 {

	ret := [3][][][]float64{}
	for i := 0; i < 3; i++ {
		ret[i] = ZeroWorldScalar64()
	}
	return ret
}

func getPhiDot() [][][]float32 {
	return DWMonitor.phidot
}

func getThetaDot() [][][]float32 {
	return DWMonitor.thetadot
}

// Find the average of old and new rxy slices.
func averageRxy(rxyNew, rxyOld [][][]float32) [][][]float32 {
	n := MeshSize()

	ret := make([][][]float32, n[Z])
	for i := 0; i < n[Z]; i++ {
		ret[i] = make([][]float32, n[Y])
		for j := 0; j < n[Y]; j++ {
			ret[i][j] = make([]float32, n[X])
			for k := 0; k < n[X]; k++ {
				ret[i][j][k] = 0.5 * (rxyNew[i][j][k] + rxyOld[i][j][k])
			}
		}
	}
	return ret
}

// Find the exact dw position by searching for zero crossings of mz along the x and y directions.
// Take the average of the x-coordinate of these zero crossings to get the DW position.
func exactPosTrace() float64 {

	c := Mesh().CellSize()
	wall := traceWall3D(M.Comp(Z).HostCopy().Scalars())
	sum := 0
	for i := range wall {
		for j := range wall[i] {
			sum += wall[i][j][1]
		}
	}
	return GetShiftPos() + (c[X] * float64(sum) / float64(len(wall)*len(wall[0])))
}

func exactPosAvg() float64 {

	ws := Mesh().WorldSize()

	// Get average magnetization; M.Comp(Z).Average() is ~ 2x faster than using my avgMz function. They don't return
	// exactly the same values, however...?
	avg := M.Comp(Z).Average()

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

// Get the exact position of the domain wall from the zero crossing of Mz.
func exactPosZC(mz [][][]float32, intPosZC [][]int) float32 {
	c := Mesh().CellSize()
	return exactPosInWindow(mz, intPosZC)*float32(c[X]) + float32(GetShiftPos())
}

// Get the indices on the left side of the domain wall within the simulation window.
func getIntDWPos(theta [][][]float32) [][]int {

	n := MeshSize()
	pos := make([][]int, n[Z])

	for i := range theta {
		pos[i] = make([]int, n[Y])
		for j := range theta[i] {
			pos[i][j] = zeroCrossing(theta[i][j])
		}
	}
	return pos
}

// Find the exact average DW position in the simulation space in units of Mesh().CellSize()[X]
func exactPosInWindow(mz [][][]float32, intPos [][]int) float32 {
	pos := float32(0.0)
	for iz := range mz {
		for iy := range mz[iz] {
			pos += interpolateZeroCrossing(mz[iz][iy], intPos[iz][iy])
		}
	}
	return pos / float32(len(intPos)*len(intPos[0]))
}

// Interpolate the index of the 1D slice wehre the zero crossing of the Mz component occurs, between x-index i and i+1.
// Returns a float32, so you can't use this for indexing an array at the zero crossing point.
func interpolateZeroCrossing(mz []float32, i int) float32 {
	return float32(i) - (mz[i] / (mz[i+1] - mz[i]))
}

func avg2D(a [][]float64) float64 {
	ret := 0.0
	for i := range a {
		for j := range a[i] {
			ret += a[i][j]
		}
	}
	return ret / float64(len(a)*len(a[0]))
}
