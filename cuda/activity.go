package cuda

import (
	"github.com/mumax/3/data"
	// "github.com/mumax/3/util"
)

func Axy(phiDot, rxyAvg, oldDWPos *data.Slice, mw int) float32 {
	n := phiDot.Size()

	// phiDot[i]*rxyAvg[i]
	_tmp := Buffer(1, n)
	defer _tmp.Free()
	Mul(_tmp, phiDot, rxyAvg)

	// Extract the values of phiDot[i]*rxyAvg[i] which lie near where the old domain wall position was
	_tmp2 := Buffer(1, [3]int{2*mw + 1, n[Y], n[Z]})
	defer _tmp2.Free()
	cfg := make3DConf([3]int{2*mw + 1, n[Y], n[Z]})
	k_nearDW_async(_tmp2.DevPtr(X), _tmp.DevPtr(X), oldDWPos.DevPtr(X), mw, n[X], n[Y], n[Z], cfg)

	// Sum the resulting slice to get Axy
	return Sum(_tmp2)
}

func Az(thetaDot, oldDWPos *data.Slice, mw int) float32 {
	n := thetaDot.Size()

	// Extract the values of thetaDot[i] which lie near where the old domain wall position was
	_tmp := Buffer(1, [3]int{2*mw + 1, n[Y], n[Z]})
	defer _tmp.Free()
	cfg := make3DConf([3]int{2*mw + 1, n[Y], n[Z]})
	k_nearDW_async(_tmp.DevPtr(X), thetaDot.DevPtr(X), oldDWPos.DevPtr(X), mw, n[X], n[Y], n[Z], cfg)

	// Sum the resulting slice to get Axy
	return Sum(_tmp)
}

func GetNearDW(a, dwpos *data.Slice, mw int) *data.Slice {
	n := a.Size()
	_tmp := Buffer(1, [3]int{2*mw + 1, n[Y], n[Z]})
	cfg := make3DConf([3]int{2*mw + 1, n[Y], n[Z]})
	k_nearDW_async(_tmp.DevPtr(X), a.DevPtr(X), dwpos.DevPtr(X), mw, n[X], n[Y], n[Z], cfg)
	return _tmp
}
