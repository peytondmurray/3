package engine

import (
	"github.com/mumax/3/data"
)

var (
	ext_dm = NewVectorField("ext_dm", "", "Change in the components of the magnetization since last time", setdm)
	SWMonitor swmonitor
)

type swmonitor struct {
	t float64
	m [3][][][]float32
	initialized bool
	lastStep int
}

func idx(ix, iy, iz, nx, ny, nz int) int {
	return ((iz)*(ny) + (iy))*nx + ix
}

func (s *swmonitor) update() {
	if !s.initialized {
		s.init()
	} else if s.t != Time || s.lastStep != NSteps {
		s.push()
	}
	return
}

func (s *swmonitor) init() {
	s.t = Time
	s.lastStep = NSteps
	s.m = M.Buffer().HostCopy().Vectors()
	return
}

func (s *swmonitor) push() {
	s.t = Time
	s.lastStep = NSteps
	s.m = M.Buffer().HostCopy().Vectors()
	return
}

func setDM(dst *data.Slice) {
	SWMonitor.update()
	n := MeshSize()
	_m := M.Buffer().HostCopy().Vectors()

	for iz := range(_m[0]) {
		for iy := range(_m[1]) {
			for ix := range(_m[2]) {
				_v := Vector(float64(_m[0][iz][iy][ix] - SWMonitor.m[0][iz][iy][ix]), float64(_m[1][iz][iy][ix] - SWMonitor.m[1][iz][iy][ix]), float64(_m[2][iz][iy][ix] - SWMonitor.m[2][iz][iy][ix]))
				dst.SetVector(ix, iy, iz, _v)
			}
		}
	}
	return
}