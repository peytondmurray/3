package engine

import (
	"github.com/mumax/3/cuda"
	"github.com/mumax/3/data"
)

var (
	ext_dm    = NewVectorField("ext_dm", "", "Change in the components of the magnetization since last time", setDM)
	SWMonitor swmonitor
)

type swmonitor struct {
	t           float64
	// m           [3][][][]float32
	// dm          [3][][][]float64
	m  			*data.Slice
	initialized bool
	lastStep    int
}

// func idx(ix, iy, iz, nx, ny, nz int) int {
// 	return ((iz)*(ny)+(iy))*nx + ix
// }

// func (s *swmonitor) update() {
// 	if !s.initialized {
// 		s.init()
// 	} else if s.t != Time || s.lastStep != NSteps {
// 		s.push()
// 	}
// 	return
// }

// func (s *swmonitor) init() {
// 	s.t = Time
// 	s.lastStep = NSteps
// 	s.m = M.Buffer().HostCopy().Vectors()
// 	s.dm = ZeroWorldVector64()
// 	s.initialized = true
// 	return
// }

// func (s *swmonitor) push() {
// 	_m := M.Buffer().HostCopy().Vectors()
// 	_t := Time
// 	_dt := _t - s.t

// 	for comp := 0; comp < 3; comp++ {
// 		for iz := range s.dm[comp] {
// 			for iy := range s.dm[comp][iz] {
// 				for ix := range s.dm[comp][iz][iy] {
// 					s.dm[comp][iz][iy][ix] = float64(_m[comp][iz][iy][ix]-s.m[comp][iz][iy][ix]) / _dt
// 				}
// 			}
// 		}
// 	}

// 	s.m = _m
// 	s.t = _t
// 	s.lastStep = NSteps

// 	return
// }

func setDM(dst *data.Slice) {

	n := MeshSize()
	if !SWMonitor.initialized {
		SWMonitor.t = Time
		SWMonitor.m = cuda.Buffer(3, n)							// Generate a new buffer for storing old values of m
		data.Copy(SWMonitor.m, M.Buffer())						// Copy m to the buffer
		data.Copy(dst, NewSlice(3, n[X], n[Y], n[Z]))			// Copy an empty slice to dst
		SWMonitor.initialized = true
	} else if SWMonitor.t != Time || SWMonitor.lastStep != NSteps {
		_t := Time
		_m := M.Buffer()
		cuda.SubDiv(dst, _m, SWMonitor.m, _t - SWMonitor.t)
		data.Copy(SWMonitor.m, _m)						// Copy m to the buffer
		SWMonitor.t = _t
	} else {
		dst = NewSlice(3, n[X], n[Y], n[Z])
	}

	return
}
