package cuda

import (
	"github.com/mumax/3/data"
)

func SubsetXRange(src *data.Slice, xmin, xmax int) *data.Slice {
	n := src.Size()
	_tmp := Buffer(1, [3]int{xmax - xmin, n[Y], n[Z]})
	cfg := make3DConf([3]int{xmax - xmin, n[Y], n[Z]})
	k_subsetXRange_async(_tmp.DevPtr(X), src.DevPtr(X), xmin, xmax, n[X], n[Y], n[Z], cfg)
	return _tmp
}
