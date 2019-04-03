package engine

import (
	"github.com/mumax/3/cuda"
	"github.com/mumax/3/data"
)

var (
	ext_phi   = NewScalarField("ext_phi", "rad", "Azimuthal angle phi", SetPhiAngle)
	ext_theta = NewScalarField("ext_theta", "rad", "Polar angle theta", SetThetaAngle)
)

func SetPhiAngle(dst *data.Slice) {
	cuda.SetPhi(dst, M.Buffer())
	return
}

func SetThetaAngle(dst *data.Slice) {
	cuda.SetTheta(dst, M.Buffer())
	return
}
