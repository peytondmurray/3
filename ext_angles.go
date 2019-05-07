package engine

import (
	"github.com/mumax/3/cuda"
	"github.com/mumax/3/data"
)

var (
	ext_phi         = NewScalarField("ext_phi", "rad", "Azimuthal angle phi", SetPhiAngle)
	ext_theta       = NewScalarField("ext_theta", "rad", "Polar angle theta", SetThetaAngle)
	ext_rxy         = NewScalarField("ext_rxy", "m", "Magnitude of m in the xy plane", SetRxyMagnitude)
	ext_rxyphitheta = NewVectorField("ext_rxyphitheta", "m[rxy] rad[phi] rad[theta]", "Magnitude of m in xy plane, azimuthal angle, polar angle", SetRxyPhiTheta)
)

func SetPhiAngle(dst *data.Slice) {
	cuda.SetPhi(dst, M.Buffer())
	return
}

func SetThetaAngle(dst *data.Slice) {
	cuda.SetTheta(dst, M.Buffer())
	return
}

func SetRxyMagnitude(dst *data.Slice) {
	cuda.SetRxy(dst, M.Buffer())
	return
}

func SetRxyPhiTheta(dst *data.Slice) {
	cuda.SetRxyPhiTheta(dst, M.Buffer())
}

func GetPhi() *data.Slice {
	s := ValueOf(ext_phi)
	defer cuda.Recycle(s)
	return s
}

func GetTheta() *data.Slice {
	s := ValueOf(ext_theta)
	defer cuda.Recycle(s)
	return s
}

func GetRxy() *data.Slice {
	s := ValueOf(ext_rxy)
	defer cuda.Recycle(s)
	return s
}
