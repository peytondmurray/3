package engine

import (
	// "log"
	// "fmt"
	"github.com/mumax/3/data"
)

var (
	DWFineSpeed = NewScalarValue("ext_dwfinespeed", "m/s", "Speed of domain wall", getDWFineSpeed) // Speed of DW
	DWFinePos   = NewScalarValue("ext_dwfinepos", "m", "Position of domain wall from start", getDWxFinePos)
	lastDWPos   float64  // Position of DW last time we checked DW speed
	lastTime    float64  // Time at the last time we checked DW speed
	lastDWSpeed float64  // Speed at the last time we checked DW speed
	lastStep    int      // Step at the last time we checked DW speed
	DWPosStack  posStack // Most recent positions of DW speed
)

// FIFO structure for storing DW positions
type posStack struct {
	data []float64
	size int
}

// Preserving the length of the stack, append an item, removing the first item.
func (s *posStack) push(v float64) {
	s.data = append(s.data[1:s.size-1], v)
	return
}

func (s *posStack) speed() float64 {

	weights := fornbergWeights(s.data[s.size]-1, s.data, 1)
	v := float64(0)
	for i := 0; i < s.size; i++ {
		v += weights[i] * s.data[i]
	}
	return v
}

// Gives the forward finite difference coefficients in a slice for a given differentiation order m
// and number of points n (which determines the order of accuracy).
// Sorry for the bad code, the notation in the original papers is just as bad.
// Fornberg, Bengt (1988), "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
// Mathematics of Computation, 51 (184): 699–706
func fornbergWeights(u float64, x []float64, k int) []float64 {

	n := len(x)
	C := make([][]float64, k+1)
	for i := 0; i < k+1; i++ {
		C[i] = make([]float64, n)
	}

	c1 := float64(1)
	c2 := float64(1)
	c3 := float64(0)
	c4 := x[1] - u
	c5 := float64(0)
	C[0][0] = 1.0

	for i := 0; i < n; i++ {
		mn := min(i, k)
		c2 = float64(1)
		c5 = c4
		c4 = x[i] - u

		for j := 0; j < i; j++ {
			c3 = x[i] - x[j]
			c2 *= c3

			if j == i-1 {
				for s := mn; s > 0; s-- {
					C[s][i] = c1 * (float64(s)*C[s-1][i-1] - c5*C[s][i-1]) / c2
				}
				C[0][i] = -c1 * c5 * C[0][i-1] / c2
			}
			for s := mn; s > 0; s-- {
				C[s][j] = (c4*C[s][j] - float64(s)*C[s-1][j]) / c3
			}
			C[0][j] = c4 * C[0][j] / c3
		}
		c1 = c2
	}

	return C[k]
}

func min(i int, j int) int {
	if i < j {
		return i
	}
	return j
}

func getDWFineSpeed() float64 {
	// fmt.Print(getDWxFinePos())
	if NSteps == 0 {
		// lastDWPos = getDWxFinePos()
		DWPosStack.push(getDWxFinePos())
	}
	if lastTime != Time {
		// currentDWPos := getDWxFinePos()
		// lastDWSpeed = (currentDWPos - lastDWPos)/(Time - lastTime)
		// lastTime = Time
		// lastStep = NSteps
		// lastDWPos = currentDWPos

		DWPosStack.push(getDWxFinePos())

	}
	// return lastDWSpeed
	return DWPosStack.speed()
}

func getDWxFinePos() float64 {
	// print(GetShiftPos())
	return _window2DDWxPos() + GetShiftPos()
}

// _window2DDWxPos finds the position of the domain wall within the simulation window
func _window2DDWxPos() float64 {
	mz := M.Buffer().Comp(Z).HostCopy().Scalars()[0]
	c := Mesh().CellSize()

	// If ShiftMagL and ShiftMagR are not specified, use the current magnetic configuration to determine what these
	// values should be; they are used in _1DDWxPos to find the zero crossing, where the domain wall is centered.
	zero := data.Vector{0, 0, 0}
	if ShiftMagL == zero || ShiftMagR == zero {
		sign := _sign32(mz[len(mz)/2][0])
		ShiftMagL[Z] = float64(sign)
		ShiftMagR[Z] = -float64(sign)
	}

	pos := _avg(_2DDWxPos(mz))
	return c[0] * float64(pos)
}

func _2DDWxPos(mz [][]float32) []float32 {
	pos := make([]float32, len(mz))
	for iy := range mz {
		pos[iy] = _1DDWxPos(mz[iy])
	}
	return pos
}

func _1DDWxPos(mz []float32) float32 {
	// Find the DW position by finding the index of the element with z-component nearest to 0
	// min := abs(mz[0])
	// pos := 0
	// for ix := range mz {
	// 	if abs(min) > abs(mz[ix]) {
	// 		pos = ix
	// 		min = abs(mz[ix])
	// 	}
	// }

	signR := _sign32(float32(ShiftMagR[Z]))
	signL := _sign32(float32(ShiftMagL[Z]))

	// print(signL, signR)

	for ix := 0; ix < len(mz)-1; ix++ {
		if _sign32(mz[ix]) == signL && _sign32(mz[ix+1]) == signR {
			return _interpolateZeroCrossing(mz, ix)
		}
	}
	panic("Can't find domain wall position.")
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

// func _backwardFiniteDifference1stOrder(arr []float64) {

// }

// func _backwardFiniteDifference8thOrder()
