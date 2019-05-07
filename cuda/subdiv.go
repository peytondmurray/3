package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func SubDiv(s, a, b *data.Slice, dt float64) {
	N := s.Size()
	util.Argument(a.Size() == N && b.Size() == N)
	cfg := make3DConf(N)
	if cfg == nil {
		panic("Invalid CUDA kernel config")
	}
	k_subdiv_async(s.DevPtr(X), s.DevPtr(Y), s.DevPtr(Z),
		a.DevPtr(X), a.DevPtr(Y), a.DevPtr(Z),
		b.DevPtr(X), b.DevPtr(Y), b.DevPtr(Z),
		float32(dt), N[X], N[Y], N[Z], cfg)
}
