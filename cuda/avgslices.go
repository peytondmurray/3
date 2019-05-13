package cuda

import (
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

func AvgSlices(dst, src1, src2 *data.Slice) {
	N := dst.Size()
	util.Assert(src1.Size()[X] == N[X])
	util.Assert(src1.Size()[Y] == N[Y])
	util.Assert(src1.Size()[Z] == N[Z])
	util.Assert(src2.Size()[X] == N[X])
	util.Assert(src2.Size()[Y] == N[Y])
	util.Assert(src2.Size()[Z] == N[Z])
	cfg := make3DConf(N)
	k_avgSlices_async(dst.DevPtr(X), src1.DevPtr(X), src2.DevPtr(X), N[X], N[Y], N[Z], cfg)
	return
}