package cuda

import (
	"unsafe"
	"github.com/mumax/3/data"
	"github.com/mumax/3/util"
)

// Add effective field of Dzyaloshinskii-Moriya interaction to Beff (Tesla).
// According to Bagdanov and Röβler, PRL 87, 3, 2001. eq.8 (out-of-plane symmetry breaking).
// See dmi.cu
func AddDMITensor(Beff *data.Slice, m *data.Slice, Aex_red,
				  Dxx_red SymmLUT, Dxy_red SymmLUT, Dxz_red SymmLUT, Dyy_red SymmLUT, Dyz_red SymmLUT, Dzz_red SymmLUT,
				  Msat MSlice, regions *Bytes, mesh *data.Mesh, OpenBC bool) {
	c := mesh.CellSize()
	N := Beff.Size()
	util.Argument(m.Size() == N)
	cfg := make3DConf(N)
	var openBC byte
	if OpenBC {
		openBC = 1
	}

	k_adddmitensor_async(Beff.DevPtr(X), Beff.DevPtr(Y), Beff.DevPtr(Z), m.DevPtr(X), m.DevPtr(Y), m.DevPtr(Z),
		Msat.DevPtr(0), Msat.Mul(0), unsafe.Pointer(Aex_red),
		unsafe.Pointer(Dxx_red), unsafe.Pointer(Dxy_red), unsafe.Pointer(Dxz_red),
		unsafe.Pointer(Dyy_red), unsafe.Pointer(Dyz_red), unsafe.Pointer(Dzz_red),
		regions.Ptr, float32(c[X]), float32(c[Y]), float32(c[Z]), N[X], N[Y], N[Z], mesh.PBC_code(), openBC, cfg)
}
