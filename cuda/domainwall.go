package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func SetDomainWallIndices(s, m *data.Slice) {
	n := m.Size()
	util.Argument(s.NComp() == 1)
	util.Argument(s.Size()[X] == 1)
	util.Argument(s.Size()[Y] == n[Y])
	util.Argument(s.Size()[Z] == n[Z])
	cfg := make3DConf([3]int{1, n[Y], n[Z]})
	k_dwPos_async(s.DevPtr(X), m.DevPtr(Z), n[X], n[Y], n[Z], cfg)
	return
}