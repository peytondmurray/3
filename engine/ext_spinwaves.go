package engine

import (
	"github.com/mumax/3/data"
)

var (
	ext_dm    = NewVectorField("ext_dm", "", "Change in the components of the magnetization since last time", setDM)
	SWMonitor swmonitor
)

type swmonitor struct {
	t           float64
	m           [3][][][]float32
	dm          [3][][][]float64
	initialized bool
	lastStep    int
}

func idx(ix, iy, iz, nx, ny, nz int) int {
	return ((iz)*(ny)+(iy))*nx + ix
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
	s.dm = ZeroWorldVector64()
	return
}

func (s *swmonitor) push() {
	_m := M.Buffer().HostCopy().Vectors()
	_t := Time
	_dt := _t - s.t

	for comp := 0; comp < 3; comp++ {
		for iz := range s.dm {
			for iy := range s.dm[iz] {
				for ix := range s.dm[iz][iy] {
					s.dm[comp][iz][iy][ix] = float64(_m[comp][iz][iy][ix]-s.m[comp][iz][iy][ix]) / _dt
				}
			}
		}
	}

	s.m = _m
	s.t = _t
	s.lastStep = NSteps

	return
}

func setDM(dst *data.Slice) {
	SWMonitor.update()
	// dst.SetVectorArray(SWMonitor.dm)
	// dst = sliceVectorsFromList(SWMonitor.dm)
	n := MeshSize()
	dst = NewSlice(3, n[X], n[Y], n[Z]).HostCopy()
	dst.SetVectorArray(SWMonitor.dm)
	return
}
