package engine

import (
	"github.com/mumax/3/data"
	// "github.com/mumax/3/mag"
	// "log"
	"fmt"
	"github.com/mumax/3/cuda"
	"math"
	// "time"
)

var (
	Az            = NewScalarValue("ext_az", "rad/s", "Out-of-plane domain wall activity", getAz)
	Axy           = NewScalarValue("ext_axy", "rad/s", "In-plane domain wall activity", getAxy)
	ExactDWVelAvg = NewScalarValue("ext_exactdwvelavg", "m/s", "Speed of domain wall", getExactVelAvg)
	ExactDWPosAvg = NewScalarValue("ext_exactdwposavg", "m", "Position of domain wall from start", getExactPosAvg)
	ExactDWPosZC  = NewScalarValue("ext_exactdwposzc", "m", "Position of the domain wall from start", getExactPosZC)
	DWMonitor     activityStack // Most recent positions of DW speed
	MaskWidth     = 10
	EdgeSpacing   = 100 // Number of cells to avoid along the left and right edges when calculating domain wall width
	AvoidEdges    = false
)

func init() {
	DeclVar("activitymaskwidth", &MaskWidth, "Number of cells on either side of the domain wall to include in activity calculations")
	DeclVar("edgespacing", &EdgeSpacing, "Number of cells to avoid on the left and right edges when calculating DW velocity")
	DeclVar("avoidedges", &AvoidEdges, "When true, DW position (used to calculate DW activity and velocity) will not include cells on the left and right edges of the simulation; the number of cells to avoid can be set by changing EdgeSpacing")

	DWMonitor.initialized = false
}

type activityStack struct {

	// DWVel
	posAvg float64
	velAvg float64

	// DWActivity
	windowpos   int
	t           float64
	lastStep    int
	maskWidth   int
	Az          float32
	Axy         float32
	initialized bool

	// Pointers to rxy, ϕ, θ, and the angular velocities ϕDot and θDot
	rϕθ   *data.Slice
	ϕDot  *data.Slice
	θDot  *data.Slice
	dwpos *data.Slice
}

func (s *activityStack) update() {
	if !s.initialized {
		s.init()
	} else if s.t != Time || s.lastStep != NSteps {
		s.push()
	}
	return
}

func getExactPosZC() float64 {
	// n := MeshSize()
	// c := Mesh().CellSize()
	// _dwposZCs := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	// _dwposExact := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	// defer _dwposZCs.Free()
	// defer _dwposExact.Free()
	// cuda.SetDomainWallIndices(_dwposZCs, M.Buffer())
	// cuda.SetExactDWPosZCInterpolated(_dwposExact, M.Buffer(), _dwposZCs)
	// ret := c[X]*float64(cuda.Sum(_dwposExact))/float64(n[Y]*n[Z]) + GetShiftPos()

	ret := getExactPosZCCPU()

	return ret
}

func getExactPosZCCPU() float64 {
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

func checkEdgeCutoff() {
	if float64(EdgeSpacing) > 0.5*float64(MeshSize()[X]) {
		print(fmt.Sprintf("WARNING: The amount of simulation being ignored for domain wall velocity/activity calculations is large compared to the size of the system! EdgeSpacing: %d, MeshSize[X]: %d", EdgeSpacing, MeshSize()[X]))
		if float64(EdgeSpacing) > float64(MeshSize()[X]) {
			panic("EdgeSpacing > MeshSize[X]")
		}
	}
	return
}

func (s *activityStack) init() {

	// Set aside buffers for holding r, phi, and theta, and phidot and thetadot
	n := MeshSize()
	s.t = Time
	s.lastStep = NSteps
	s.windowpos = GetIntWindowPos()
	if AvoidEdges {
		checkEdgeCutoff()
		s.posAvg = exactPosAvgAvoidEdges(EdgeSpacing)
	} else {
		s.posAvg = exactPosAvg()
	}
	s.velAvg = 0.0

	s.rϕθ = cuda.Buffer(3, MeshSize())
	ext_rxyphitheta.EvalTo(s.rϕθ)
	s.ϕDot = cuda.Buffer(1, MeshSize())
	cuda.Zero(s.ϕDot)
	s.θDot = cuda.Buffer(1, MeshSize())
	cuda.Zero(s.θDot)
	s.dwpos = cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	cuda.SetDomainWallIndices(s.dwpos, M.Buffer())

	s.initialized = true

	return
}

func (s *activityStack) push() {

	n := MeshSize()
	_t := Time

	_posAvg := 0.0
	// DWVel_________________________
	if AvoidEdges {
		_posAvg = exactPosAvgAvoidEdges(EdgeSpacing)
	} else {
		_posAvg = exactPosAvg()
	}
	s.velAvg = (_posAvg - s.posAvg) / (_t - s.t)
	s.posAvg = _posAvg
	// DWVel_________________________

	// Get new window and domain wall positions. Get the newest rxy, phi, theta values.
	_windowpos := GetIntWindowPos()

	// Allocate GPU memory to hold intermediate quantities
	_dwpos := cuda.Buffer(1, [3]int{1, n[Y], n[Z]})
	_rxyAvgGPU := cuda.Buffer(1, n)
	_ϕDot := cuda.Buffer(1, n)
	_θDot := cuda.Buffer(1, n)
	_rϕθ := cuda.Buffer(3, n)
	defer cuda.Recycle(_rxyAvgGPU) // Free it when we return from the function, we don't use it after calculating Axy

	ext_rxyphitheta.EvalTo(_rϕθ)
	cuda.SetDomainWallIndices(_dwpos, M.Buffer())

	// Average the rxy from last step and this step
	cuda.AvgSlices(_rxyAvgGPU, _rϕθ.Comp(0), s.rϕθ.Comp(0))
	cuda.SubDivAngle(_ϕDot, _rϕθ.Comp(1), s.rϕθ.Comp(1), _windowpos-s.windowpos, _t-s.t) // ang. vel. phi
	cuda.SubDivAngle(_θDot, _rϕθ.Comp(2), s.rϕθ.Comp(2), _windowpos-s.windowpos, _t-s.t) // ang. vel. theta

	s.Axy = cuda.Axy(_ϕDot, _rxyAvgGPU, s.dwpos, MaskWidth)
	s.Az = cuda.Az(_θDot, s.dwpos, MaskWidth)

	s.dwpos.Free()
	s.rϕθ.Free()
	s.ϕDot.Free()
	s.θDot.Free()

	// To be done after the other stuff!
	s.windowpos = _windowpos
	s.dwpos = _dwpos
	s.rϕθ = _rϕθ
	s.t = _t
	s.lastStep = NSteps
	s.ϕDot = _ϕDot
	s.θDot = _θDot

	return
}

func IntRound(x float32) int {
	return int(x + 0.5)
}

// Find the index of the 1D slice just before the zero crossing of the Mz component. Check from right to left;
// since bubbles will only be forming to the left, and we will not be susceptible to them with this method.
func zeroCrossing(theta []float32) int {
	for ix := len(theta) - 1; ix > 1; ix-- {
		if math.Cos(float64(theta[ix-1]))*math.Cos(float64(theta[ix])) < 0 {
			return ix - 1
		}
	}
	panic("Can't find domain wall position")
}

// Number of cells the window has shifted
func GetIntWindowPos() int {
	c := Mesh().CellSize()
	windowPos := GetShiftPos()
	return IntRound(float32(windowPos) / float32(c[Y]))
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

func exactPosAvgAvoidEdges(edgeSpacing int) float64 {

	n := MeshSize()
	c := Mesh().CellSize()

	_tmp := cuda.SubsetXRange(M.Buffer().Comp(Z), edgeSpacing, n[X]-edgeSpacing)
	defer _tmp.Free()

	avg := float64(cuda.Sum(_tmp)) / float64(n[Z]*n[Y]*(n[X]-2*edgeSpacing))
	pct := 1.0 - (1.0-avg)/2.0
	posInMiddle := pct * float64(n[X]-2*edgeSpacing) * c[X]
	lowerPosBound := float64(edgeSpacing) * c[X]

	return posInMiddle + lowerPosBound + GetShiftPos()
}

// Given the average value of Mz, compute the absolute domain wall position.
func exactPosAvg() float64 {

	ws := Mesh().WorldSize()

	// Get average magnetization; M.Comp(Z).Average() is ~ 2x faster than using my avgMz function
	avg := M.Comp(Z).Average()

	// Percentage of the magnetization which is flipped up gives the position of the domain wall, for example if
	// 50% are flipped up, the DW is 50% from the left side of the simulation window
	pct := 1.0 - (1.0-avg)/2.0

	// Convert to actual position in window, then add on window shift
	return pct*ws[X] + GetShiftPos()
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
