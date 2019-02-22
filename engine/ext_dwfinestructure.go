package engine

// "log"
// "fmt"
// "github.com/mumax/3/data"

var (
	DWFineSpeed = NewScalarValue("ext_dwfinespeed", "m/s", "Speed of domain wall", getDWFineSpeed) // Speed of DW
	DWFinePos   = NewScalarValue("ext_dwfinepos", "m", "Position of domain wall from start", getDWxFinePos)
	DWPosStack  posStack // Most recent positions of DW speed
	SignL       int      // Store the sign of the z-component of the magnetization inserted at the Left side when sim window shfits
	SignR       int      // Store the sign of the z-component of the magnetization inserted at the Right side when sim window shfits
)

func init() {
	DeclFunc("ext_dwfineposinit", DWFinePosInit, "DWFinePosInit(q, l, r)  sets the order of the DW velocity calculation to order q, and the sign of the magnetization which is being inserted at the left and right sides of the simulation. (order, signL, signR)")
}

// FIFO structure for storing DW positions
type posStack struct {
	t		[2]float64
	pos     [2]float64
	initialized bool
}

func DWFinePosInit(signL int, signR int) {
	SignL = signL
	SignR = signR
	DWPosStack.initialized = false
	return
}

// Preserving the length of the stack, append an item, removing the first item. If the stack isn't full, just append.
func (s *posStack) push(t float64, pos float64) {

	s.t[0] = s.t[1]
	s.t[1] = t

	s.pos[0] = s.pos[1]
	s.pos[1] = pos
	return
}

func (s *posStack) speed() float64 {
	return (s.pos[1]-s.pos[0])/(s.t[1]-s.t[0])
}

func (s *posStack) lastTime() float64 {
	return s.t[1]
}

// Get the minimum of two integers
func min(i int, j int) int {
	if i < j {
		return i
	}
	return j
}

// Return the speed of the DW. If the DW position hasn't been sampled enough, sample it and return 0 for the speed.
func getDWFineSpeed() float64 {
	if !DWPosStack.initialized {
		DWPosStack.initialized = true
		DWPosStack.push(Time, getDWxFinePos())
		return 0
	} else if DWPosStack.lastTime() != Time {
		DWPosStack.push(Time, getDWxFinePos())
	}
	return DWPosStack.speed()
}

// The fine position of the DW is the position of the DW in the window plus the shift of the window relative to the lab
func getDWxFinePos() float64 {
	return _window2DDWxPos() + GetShiftPos()
}

// _window2DDWxPos finds the position of the domain wall within the simulation window
func _window2DDWxPos() float64 {
	mz := M.Buffer().Comp(Z).HostCopy().Scalars()[0]
	c := Mesh().CellSize()
	pos := _avg(DWFinexPos2D(mz))
	return c[0] * float64(pos)
}

// Find the DW position at each row of a 2D simulation space
func DWFinexPos2D(mz [][]float32) []float32 {
	pos := make([]float32, len(mz))
	for iy := range mz {
		pos[iy] = DWFinexPos1D(mz[iy])
	}
	return pos
}

// Find the DW position along a 1D slice of the simulation region, interpolating across the zero crossing.
func DWFinexPos1D(mz []float32) float32 {
	return _interpolateZeroCrossing(mz, ZeroCrossing(mz))
	panic("Can't find domain wall position.")
}

// Find the index of the 1D slice where the zero crossing of the Mz component occurs.
func ZeroCrossing(mz []float32) int {
	for ix := 0; ix < len(mz)-1; ix++ {
		if _sign32(mz[ix]) == SignL && _sign32(mz[ix+1]) == SignR {
			return ix
		}
	}
	panic("Can't find domain wall position")
}

func _interpolateZeroCrossing(mz []float32, i int) float32 {
	return float32(i) - (mz[i] / (mz[i+1] - mz[i]))
}

func _sign32(f float32) int {
	if f >= 0 {
		return 1
	}
	return -1
}

func _avg(s []float32) float32 {
	sum := float32(0)
	for v := range s {
		sum += s[v]
	}
	return sum / float32(len(s))
}
