package engine

import (
	"fmt"
	"os"
)

var (
	save2DIndex int
	saveXYZIndex int // Stores the number of the next file to be written. Increments each time a file is written.
	saveDWConfigIndex int
)

func init() {
	DeclFunc("ext_save2D", save2D, "Save scalar slice of [][]float64 values")
	DeclFunc("ext_saveVTK", saveVTK, "Save scalar slice of [][][]float64 values")
	DeclFunc("ext_saveDWConfig", saveDWConfig, "Save domain wall configuration as csv.")
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

// Writes a slice of scalars as an xyz file.
func saveVTK(s [][][]float64, name string) {

	n := MeshSize()
	ncells := n[X] * n[Y] * n[Z]
	c := Mesh().CellSize()

	// Find current working dir
	f, err := os.Create(fmt.Sprintf(OD()+name+"%06d.vtk", saveXYZIndex)) //https://gobyexample.com/writing-files
	check(err)

	// Write timestamp
	// f.WriteString(fmt.Sprintf("#time = %e\n", Time))
	// f.WriteString(fmt.Sprintf("xi,yi,zi,value\n"))
	f.WriteString("# vtk DataFile Version 3.1\n")
	f.WriteString("Written by ext_saveVTK\n")
	f.WriteString("ASCII\n")
	f.WriteString("DATASET STRUCTURED_GRID\n")
	f.WriteString(fmt.Sprintf("DIMENSIONS %d %d %d\n", n[X], n[Y], n[Z]))
	f.WriteString(fmt.Sprintf("POINTS %d float\n", ncells))

	// Write values
	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < n[X]; i++ {
				// f.WriteString(fmt.Sprintf("%d,%d,%d,%e\n", i, j, k, s[k][j][i]))
				f.WriteString(fmt.Sprintf("%E %E %E\n", float64(i)*c[X], float64(j)*c[Y], float64(k)*c[Z]))
			}
		}
	}

	f.WriteString(fmt.Sprintf("\nPOINT_DATA %d\n", ncells))
	f.WriteString(fmt.Sprintf("SCALARS %s float\n", name))
	f.WriteString("LOOKUP_TABLE default\n")
	for k := 0; k < n[Z]; k++ {
		for j := 0; j < n[Y]; j++ {
			for i := 0; i < n[X]; i++ {
				f.WriteString(fmt.Sprintf("%E\n", s[k][j][i]))
			}
		}
	}

	f.Sync()
	saveXYZIndex++
}

func saveDWConfig(name string) {

	n := MeshSize()
	c := Mesh().CellSize()
	dwpos := DWMonitor.dwpos

	m := M.Buffer().HostCopy().Vectors()
	dw := make([][6]float64, n[Z]*n[Y])		// (x, y, z, mx, my, mz) along the domain wall
	k := 0
	for i:=0; i<n[Z]; i++ {
		for j:=0; j<n[Y]; j++ {
			x := c[X]*float64(dwpos[i][j])
			y := c[Y]*float64(j)
			z := c[Z]*float64(i)
			mx := float64(m[X][i][j][dwpos[i][j]])
			my := float64(m[Y][i][j][dwpos[i][j]])
			mz := float64(m[Z][i][j][dwpos[i][j]])
			dw[k] = [6]float64{x, y, z, mx, my, mz}
			k++
		}
	}

	// Find current working dir
	f, err := os.Create(fmt.Sprintf(OD()+name+"%06d.csv", saveDWConfigIndex)) //https://gobyexample.com/writing-files
	check(err)

	// Write timestamp
	f.WriteString(fmt.Sprintf("#time = %e\n", Time))
	f.WriteString(fmt.Sprintf("x,y,z,mx,my,mz\n"))

	// Write values
	for i := 0; i < len(dw); i++ {
		f.WriteString(fmt.Sprintf("%E,%E,%E,%E,%E,%E\n", dw[i][0], dw[i][1], dw[i][2], dw[i][3], dw[i][4], dw[i][5]))
	}

	f.Sync()
	saveDWConfigIndex++
}

func save2D(s [][]float64, name string) {

	// Find current working dir
	f, err := os.Create(fmt.Sprintf(OD()+name+"%06d.csv", saveDWConfigIndex))
	check(err)

	// Write timestamp
	f.WriteString(fmt.Sprintf("#time = %e\n", Time))
	f.WriteString(fmt.Sprintf("i,j,value\n"))

	// Write values
	for i := range s {
		for j := range s[i] {
			f.WriteString(fmt.Sprintf("%d,%d,%E\n", i, j, s[i][j]))
		}
	}

	f.Sync()
	save2DIndex++
}