package engine

import (
	"fmt"
	"os"
)

var (
	saveXYZIndex int // Stores the number of the next file to be written. Increments each time a file is written.
	saveDWIndex int
)

func init() {
	DeclFunc("ext_saveVTK", saveVTK, "Save scalar slice of [][][]float64 values")
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