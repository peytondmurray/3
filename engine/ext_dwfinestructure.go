package engine

import (
	//"fmt"
	// "github.com/mumax/3/data"
)

var (
	DWFineSpeed = NewScalarValue("ext_dwfinespeed", "m/s", "Speed of domain wall", getDWFineSpeed) // Speed of DW
	DWFinePos = NewScalarValue("ext_dwfinepos", "m", "Position of domain wall from start", getDWxFinePos)
	lastDWPos	float64			// Position of DW last time we checked DW speed
	lastTime	float64			// Time at the last time we checked DW speed
	lastDWSpeed float64			// Speed at the last time we checked DW speed
	lastStep	int				// Step at the last time we checked DW speed
)

func getDWFineSpeed() float64 {
	if NSteps == 0 {
		lastDWPos = getDWxFinePos()
		lastTime = Time
		lastStep = NSteps
		lastDWSpeed = 0
	}
	if lastStep != NSteps {
		currentDWPos := getDWxFinePos()
		lastDWSpeed = (currentDWPos - lastDWPos)/(Time - lastTime)
		lastTime = Time
		lastStep = NSteps
		lastDWPos = currentDWPos
	}
	return lastDWSpeed
}

func getDWxFinePos() float64 {
	return TotalShift + _window2DDWxPos()
}

// _window2DDWxPos finds the position of the domain wall within the simulation window
func _window2DDWxPos() float64 {
	mz := M.Buffer().Comp(Z).HostCopy().Scalars()[0]
	c := Mesh().CellSize()
	pos := _avg(_2DDWxPos(mz))
	return c[0]*float64(pos)
}

func _2DDWxPos(mz [][]float32) []int {
	pos := make([]int, len(mz))
	for iy := range mz {
		pos[iy] = _1DDWxPos(mz[iy])
	}
	return pos
}

func _1DDWxPos(mz []float32) int {
	min := abs(mz[0])
	pos := 0
	for ix := range mz {
		if abs(min) > abs(mz[ix]) {
			pos = ix
			min = abs(mz[ix])
		}
	}

	return pos
}

func _avg(s []int) float32 {
	sum := 0
	for v := range s {
		sum += s[v]
	}
	return float32(sum)/float32(len(s))
}