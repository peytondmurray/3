package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

// Set s to the azimuthal angle phi
// See phi.cu
func SetPhi(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	k_setphi_async(s.DevPtr(X), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
	return
}

// Set s to the polar angle theta
// See theta.cu
func SetTheta(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	k_settheta_async(s.DevPtr(X), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
	return
}

func SetRxy(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	k_setrxy_async(s.DevPtr(X), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
	return
}

func SetRxyPhiTheta(s *data.Slice, m *data.Slice) {
	N := s.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	k_setrxyphitheta_async(s.DevPtr(X), s.DevPtr(Y), s.DevPtr(Z), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z), N[X], N[Y], N[Z], cfg)
	return
}
