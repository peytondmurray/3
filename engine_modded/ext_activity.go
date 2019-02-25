package engine

import (
	"math"
	//"unsafe"
	//"fmt"
	//"github.com/mumax/3/cuda"
	//"github.com/mumax/3/data"
)

func init() {
	//DeclFunc("ext_mdwband", mDWBand, "magnetization of DW band cells")
	DeclFunc("ext_bandanglesxy", BandAnglesXY, "In-plane angle of DW band cells")
	DeclFunc("ext_bandanglesz", BandAnglesZ, "Out-of-plane angles of DW band cells")
	DeclFunc("ext_bandlenxy", BandLenXY, "DW band in-plane vector lengths")
	DeclFunc("ext_angularvelocity", AngularVelocity, "Angular velocity of cells")
	DeclFunc("ext_sumthetadot", SumThetaDot, "Sum angular velocity of cells")
	DeclFunc("ext_sumphidot", SumPhiDot, "Sum angular velocity of cells")
	DeclFunc("ext_sumthetadotabs", SumThetaDotAbs, "Sum angular velocity of cells")
	DeclFunc("ext_sumphidotabs", SumPhiDotAbs, "Sum angular velocity of cells")
	DeclFunc("ext_allband", anglesDwBand, "Return length of xy and angles xy and z")
	DeclFunc("ext_allband2", anglesDwBandNoMoveCorrection, "Return length of xy and angles xy and z")
	DeclFunc("ext_ftoi", FloatToInt, "Self-explanatory")
	DeclFunc("ext_round", Round, "Self-explanatory?")
	DeclFunc("ext_dwlocdiff", dwLocDiff, "Self-explanatory")
	DeclFunc("ext_dwposmoved", dwPosMoved, "Self-explanatory")
	
}

//func anglesDwBand(indMin [][]int, dx int) ([][]float64, [][]float64, [][]float64) {
	//if Mesh().Size()[Z] > 1 {
		//panic("ERROR: count zcell > 1. Implemented only for 2D.")
	//}
	//M := &M
	//ny := Mesh().Size()[Y]
	//absXYDwBand := make([][]float64, ny)
	//phiDwBand := make([][]float64, ny)
	//thetaDwBand := make([][]float64, ny)
	
	//for i:=0; i<ny; i++ {
		//absXYDwBand[i] = make([]float64, 2*dx)
		//phiDwBand[i] = make([]float64, 2*dx)
		//thetaDwBand[i] = make([]float64, 2*dx)
		//k := 0
		//for j:=indMin[0][i]-dx; j<indMin[0][i]+dx; j++ {
			//m := M.GetCell(j,i,0)
			//absXYDwBand[i][k] = math.Sqrt(math.Pow(m[X],2)+math.Pow(m[Y],2))
			//phiDwBand[i][k] = math.Atan2(m[Y],m[X])
			//thetaDwBand[i][k] = math.Acos(m[Z])
			//k++
		//}
	//}
	
	//return absXYDwBand, phiDwBand, thetaDwBand
//}
func FloatToInt(x float64) int {
	return int(x)
}

func Round(val float64, roundOn float64, places int ) (newVal float64) {
	var round float64
	pow := math.Pow(10, float64(places))
	digit := pow * val
	_, div := math.Modf(digit)
	_div := math.Copysign(div, val)
	_roundOn := math.Copysign(roundOn, val)
	if _div >= _roundOn {
		round = math.Ceil(digit)
	} else {
		round = math.Floor(digit)
	}
	newVal = round / pow
	return
}

func dwLocDiff(d0,d1 [][]int) []int {
	Ny := Mesh().Size()[Y]
	dwLocDiff := make([]int, Ny)
	for j:=0; j<Ny; j++{
		dwLocDiff[j] = d1[0][j]-d0[0][j]
	}
	return dwLocDiff
}
func anglesDwBand(indMin [][]int, moveCorrection []int, dx, windowCorrection int) [][][]float64 {
	if Mesh().Size()[Z] > 1 {
		panic("ERROR: count zcell > 1. Implemented only for 2D.")
	}
	M := &M
	ny := Mesh().Size()[Y]
	len_phi_theta := make([][][]float64, 3)
	for i:=0; i<3; i++ {
		len_phi_theta[i] = make([][]float64, ny)
	}
	
	for i:=0; i<ny; i++ {
		for ind:=0; ind<3; ind++{
			len_phi_theta[ind][i] = make([]float64, 2*dx)
		}
		k := 0
		for j:=indMin[0][i]-dx-windowCorrection-moveCorrection[i]; 
		j<indMin[0][i]+dx-windowCorrection-moveCorrection[i]; j++ {
			m := M.GetCell(j,i,0)
			len_phi_theta[0][i][k] = math.Sqrt(math.Pow(m[X],2)+math.Pow(m[Y],2))
			len_phi_theta[1][i][k] = math.Atan2(m[Y],m[X])
			len_phi_theta[2][i][k] = math.Acos(m[Z])
			k++
		}
	}
	
	return len_phi_theta
}

func anglesDwBandNoMoveCorrection(indMin [][]int, dx, windowCorrection int) [][][]float64 {
	if Mesh().Size()[Z] > 1 {
		panic("ERROR: count zcell > 1. Implemented only for 2D.")
	}
	M := &M
	ny := Mesh().Size()[Y]
	len_phi_theta := make([][][]float64, 3)
	for i:=0; i<3; i++ {
		len_phi_theta[i] = make([][]float64, ny)
	}
	
	for i:=0; i<ny; i++ {
		for ind:=0; ind<3; ind++{
			len_phi_theta[ind][i] = make([]float64, 2*dx)
		}
		k := 0
		for j:=indMin[0][i]-dx-windowCorrection; 
		j<indMin[0][i]+dx-windowCorrection; j++ {
			m := M.GetCell(j,i,0)
			len_phi_theta[0][i][k] = math.Sqrt(math.Pow(m[X],2)+math.Pow(m[Y],2))
			len_phi_theta[1][i][k] = math.Atan2(m[Y],m[X])
			len_phi_theta[2][i][k] = math.Acos(m[Z])
			k++
		}
	}
	
	return len_phi_theta
}

//func mDWBand(indMin [][]int, dx int) [][]data.Vector  {
	/* 	Input parameters: indMin = index array of DW position,
	 * 	dx = size of the band in cells (+/- dx around DW).
	 * 
	 *  Creates dx sized band around DW and returns magnetization
	 * 	slice.
	*/
	/*
	nx,ny := Mesh().Size()[X],Mesh().Size()[Y]
	DwBand := make([][]data.Vector, ny)
	
	for i:=0; i<ny; i++ {
		DwBand[i] = make([]data.Vector, 2*dx)
		
		for j:=0; j<2*dx; j++ {
			
			for j:=indMin[0][i]-dx; j<indMin[0][i]+dx; j++ {
				m := M.GetCell(j,i,0)
				DwBand[j][i] = m
			}
		}
	}
	return DwBand
}
*/
func BandAnglesXY(indMin [][]int, dx int) [][]float64 {
	
	if Mesh().Size()[Z] > 1 {
		panic("ERROR: count zcell > 1. Implemented only for 2D.")
	}
	
	M := &M
	ny := Mesh().Size()[Y]
	phiXYDwBand := make([][]float64, ny)

	for i:=0; i<ny; i++ {
		phiXYDwBand[i] = make([]float64, 2*dx)
		k := 0
		for j:=indMin[0][i]-dx; j<indMin[0][i]+dx; j++ {
				m := M.GetCell(j,i,0)
				phiXYDwBand[i][k] = math.Atan2(m[Y],m[X])
				k++
		}
	}
	return phiXYDwBand
}

func BandAnglesZ(indMin [][]int, dx int) [][]float64 {
	
	if Mesh().Size()[Z] > 1 {
		panic("ERROR: count zcell > 1. Implemented only for 2D.")
	}
	
	M := &M
	ny := Mesh().Size()[Y]
	phiZDwBand := make([][]float64, ny)

	for i:=0; i<ny; i++ {
		phiZDwBand[i] = make([]float64, 2*dx)
		k := 0
		for j:=indMin[0][i]-dx; j<indMin[0][i]+dx; j++ {
				m := M.GetCell(j,i,0)
				phiZDwBand[i][k] = math.Acos(m[Z])
				k++
		}
	}
	return phiZDwBand
}

func BandLenXY(indMin [][]int, dx int) [][]float64 {
	
	if Mesh().Size()[Z] > 1 {
		panic("ERROR: count zcell > 1. Implemented only for 2D.")
	}
	
	M := &M
	ny := Mesh().Size()[Y]
	LenXYDwBand := make([][]float64, ny)

	for i:=0; i<ny; i++ {
		LenXYDwBand[i] = make([]float64, 2*dx)
		k := 0
		for j:=indMin[0][i]-dx; j<indMin[0][i]+dx; j++ {
				m := M.GetCell(j,i,0)
				LenXYDwBand[i][k] = math.Sqrt(math.Pow(m[X],2)+math.Pow(m[Y],2))
				k++
		}
	}
	return LenXYDwBand
}

func AngularVelocity(angle0, angle1 [][]float64, timeDif float64, dx int) [][]float64 {
	ny := Mesh().Size()[Y]
	phidot := make([][]float64, ny)
	for i:=0; i<ny; i++ {
		phidot[i] = make([]float64, 2*dx)
		for j:=0; j<2*dx; j++ {
			dangle := angle1[i][j] - angle0[i][j]
			if dangle < -math.Pi {
				dangle += math.Pi*2
			}
			if dangle > math.Pi {
				dangle = math.Pi*2 - dangle
			}
			phidot[i][j] = (dangle)/timeDif
		}
	}
	return phidot
}

func SumThetaDotAbs(thetadot [][]float64, dx int) float64 {
	var sum float64 = 0.0
	ny := Mesh().Size()[Y]
	for i:=0; i<ny; i++ {
		for j:=0; j<2*dx; j++ {
				sum += math.Abs(thetadot[i][j])
		}
	}
	return sum
}

func SumPhiDotAbs(phidot, len0, len1 [][]float64, dx int) float64 {
	var sum float64 = 0.0
	ny := Mesh().Size()[Y]
	for i:=0; i<ny; i++ {
		for j:=0; j<2*dx; j++ {
				sum += math.Abs((len1[i][j]+len0[i][j])/2*phidot[i][j])
		}
	}
	return sum
}

func SumThetaDot(thetadot [][]float64, dx int) []float64 {
	//var sum float64 = 0.0
	sum := make([]float64, 2)
	sum[0] = 0.0
	sum[1] = 0.0
	ny := Mesh().Size()[Y]
	for i:=0; i<ny; i++ {
		for j:=0; j<2*dx; j++ {
				sum[0] += thetadot[i][j]
				sum[1] += math.Abs(thetadot[i][j])
		}
	}
	return sum
}

func SumPhiDot(phidot, len0, len1 [][]float64, dx int) []float64 {
	//var sum float64 = 0.0
	sum := make([]float64, 2)
	sum[0] = 0.0
	sum[1] = 0.0
	ny := Mesh().Size()[Y]
	for i:=0; i<ny; i++ {
		for j:=0; j<2*dx; j++ {
				sum[0] += (len1[i][j]+len0[i][j])/2*phidot[i][j]
				sum[1] += math.Abs((len1[i][j]+len0[i][j])/2*phidot[i][j])
		}
	}
	return sum
}

func dwPosMoved(indMin0, indMin1 [][]int) float64 {
	ny := Mesh().Size()[Y]
	d := 0.0
	for i:=0; i<ny; i++ {
		if indMin0[0][i] != indMin1[0][i] {
			d += 1.0
		}
	}
	return d
}
