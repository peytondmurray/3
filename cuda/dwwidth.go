package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func SetDWWidth(s *data.Slice, m *data.Slice, halfWidth int) {
	N := m.Size()
	util.Argument(s.Size()[X] == 1)
	util.Argument(s.Size()[Y] == N[Y])
	util.Argument(s.Size()[Z] == N[Z])
	cfg := make3DConf(N)
	k_avgDWWidth_async(s.DevPtr(X), m.DevPtr(Z), halfWidth, N[X], N[Y], N[Z], cfg)
	return
}
