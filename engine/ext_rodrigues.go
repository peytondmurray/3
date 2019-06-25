package engine

import (
	"fmt"
	"os"
)

var (
	saveRodriguesIndex int
)

func init() {
	DeclFunc("ext_saverodrigues", saveRodrigues, "Save the set of rotation axes and angles needed to rotate the magnetization from (0, 0, 1) to the current orientation")
	DeclFunc("ext_savecellcenterlocs", saveCellCenterLocs, "Save the locations of the cell centers")
}

func saveRodrigues() {
	n := MeshSize()
	m := M.Buffer().HostCopy().Vectors()
	rϕθ := ext_rxyphitheta.HostCopy().Vectors()

	// Find current working dir
	f, err := os.Create(fmt.Sprintf(OD()+"rodrigues%06d.csv", saveRodriguesIndex)) //https://gobyexample.com/writing-files
	check(err)

	// Write timestamp
	f.WriteString(fmt.Sprintf("#time = %E\n", Time))
	f.WriteString(fmt.Sprintf("#ix,iy,iz,nx,ny,nz,angle\n")) //nx, ny, nz are the axis around which to rotate.

	// Write values
	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := 0; k < n[X]; k++ {
				f.WriteString(fmt.Sprintf("%d,%d,%d,%E,%E,%E,%E\n", k, j, i, m[Y][i][j][k]/rϕθ[0][i][j][k], -m[X][i][j][k]/rϕθ[0][i][j][k], 0.0, rϕθ[2][i][j][k]))
			}
		}
	}

	f.Sync()
	saveRodriguesIndex++
	return
}

func saveCellCenterLocs() {
	n := MeshSize()
	c := Mesh().CellSize()

	// Find current working dir
	f, err := os.Create(fmt.Sprintf(OD() + "centerlocs.csv")) //https://gobyexample.com/writing-files
	check(err)

	// Write timestamp
	f.WriteString(fmt.Sprintf("#ix,iy,iz,x,y,z\n"))

	for i := 0; i < n[Z]; i++ {
		for j := 0; j < n[Y]; j++ {
			for k := 0; k < n[X]; k++ {
				f.WriteString(fmt.Sprintf("%d,%d,%d,%E,%E,%E\n", k, j, i, (float64(k)+0.5)*c[X], (float64(j)+0.5)*c[Y], (float64(i)+0.5)*c[Z]))
			}
		}
	}
	f.Sync()
}
