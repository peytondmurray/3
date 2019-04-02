package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

// Set s to the azimuthal angle phi
// See angles.cu
func SetPhi(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	if cfg == nil {
		panic("hey")
	}
	k_setphi_async(s.DevPtr(X), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
}

// Set s to the polar angle theta
// See angles.cu
func SetTheta(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	if cfg == nil {
		panic("hey")
	}
	k_settheta_async(s.DevPtr(X), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
}
