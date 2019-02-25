package engine

import (
	"math"
	//"unsafe"
	//"fmt"
	//"github.com/mumax/3/cuda"
	//"github.com/mumax/3/data"
)

var (
	//allMoved	bool
	//initialPosition
)

func init() {
	DeclFunc("ext_exactdwpos", exactDwPos, "Define DW position by accuray of a cell")
	DeclFunc("ext_dwposproximity", findFromProximity, "Define DW position by scanning around the previous location of DW")
	DeclFunc("ext_dwposmeters", interpolatedLocation, "DW position in meters, interpolated around the center min(m_z)")
	//DeclVar("ext_allmoved", &allMoved, "All segments of DW have moved")
}



func exactDwPos() [][]int {
	//SetBusy(true)
	//defer SetBusy(false)
	//n := Mesh().Size()
	dummy := make([][]int, 0)
	return DwPosInd(dummy,0)
} 

func findFromProximity(indexArray [][]int, findArea int) [][]int {
	return DwPosInd(indexArray, findArea)

}

func DwPosInd(indMin [][]int, dx int) [][]int {
	M := &M
	nx,ny,nz := Mesh().Size()[X],Mesh().Size()[Y],Mesh().Size()[Z]
	dwIndex := make([][]int, nz)
	for i:=0; i<nz; i++{
		dwIndex[i] = make([]int, ny)
		for j:=0; j<ny; j++{
			minM := math.Inf(1)
			minMIndex := -1
			//dwX := float64(0.0)
			if len(indMin) == 0 {
				for k:=0; k<nx; k++{
					m := M.GetCell(k,j,i)
					if math.Abs(m[Z]) < minM {
						minM = math.Abs(m[Z])
						minMIndex = k
					}	
				}
			} else {
				for k:=indMin[i][j]-dx; k<indMin[i][j]+dx; k++{
					m := M.GetCell(k,j,i)
					if math.Abs(m[Z]) < minM {
						minM = math.Abs(m[Z])
						minMIndex = k
					}
				}
			}
			dwIndex[i][j] = minMIndex
		}
	}
	return dwIndex
}

func interpolatedLocation(indMin [][]int) [][]float64 {
		
	M := &M
	c := Mesh().CellSize()
	nx,ny,nz := Mesh().Size()[X],Mesh().Size()[Y],Mesh().Size()[Z]
	dwpos := make([][]float64, nz)
	var x0,x1,y0,y1 float64

	for i:=0; i<nz; i++{
		dwpos[i] = make([]float64, ny)
		for j:=0; j<ny; j++{
			minMIndex := indMin[i][j]
			mclose := M.GetCell(minMIndex,j,i)[Z]
			mprev := M.GetCell(minMIndex-1,j,i)[Z]
			mnext := M.GetCell(minMIndex+1,j,i)[Z]
			
			if mclose < 0.0 {
				x0 = float64(minMIndex-1-nx/2)*c[X]
				y0 = mprev
				x1 = float64(minMIndex-nx/2)*c[X]
				y1 = mclose
				
			} else {
				x0 = float64(minMIndex-nx/2)*c[X]
				y0 = mclose
				x1 = float64(minMIndex+1-nx/2)*c[X]
				y1 = mnext
			}

			dwX := x1 + y1*(x0-x1)/(y1-y0)
			dwpos[i][j] = dwX + GetShiftPos()

		}
	}
	return dwpos
}
