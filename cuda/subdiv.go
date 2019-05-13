package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func SubDiv(s, a, b *data.Slice, dt float64) {
	N := s.Size()
	util.Argument(a.Size() == N && b.Size() == N)
	cfg := make3DConf(N)
	k_subDiv_async(s.DevPtr(X), s.DevPtr(Y), s.DevPtr(Z),
		a.DevPtr(X), a.DevPtr(Y), a.DevPtr(Z),
		b.DevPtr(X), b.DevPtr(Y), b.DevPtr(Z),
		float32(dt), N[X], N[Y], N[Z], cfg)
}

func SubDivAngle(s, a, b *data.Slice, shift int, dt float64) {
	N := s.Size()
	util.Argument(a.Size() == N && b.Size() == N)
	cfg := make3DConf(N)
	k_subDivAngle_async(s.DevPtr(X), a.DevPtr(X), b.DevPtr(X),
						shift, float32(dt), N[X], N[Y], N[Z], cfg)
}