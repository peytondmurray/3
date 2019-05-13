package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func SetDWWidth(s *data.Slice, m *data.Slice, halfWidth int) {
	n := m.Size()
	util.Argument(s.Size()[X] == 1)
	util.Argument(s.Size()[Y] == n[Y])
	util.Argument(s.Size()[Z] == n[Z])
	cfg := make3DConf([3]int{1, n[Y], n[Z]})
	k_DWWidth_async(s.DevPtr(X), m.DevPtr(Z), halfWidth, n[X], n[Y], n[Z], cfg)
	return
}
