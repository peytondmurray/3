package engine

import (
	//"math"
	//"unsafe"
	"fmt"
	"os"
	//"github.com/mumax/3/cuda"
	//"github.com/mumax/3/data"
)

var (
	//initialPosition
	saveIndex int
)

func init() {
	DeclFunc("ext_saveslice", saveSlice, "Save slice")
	DeclFunc("ext_save_dw_m", saveDwM, "Saves dw configuration and magnetizations mx & my")
}

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func saveSlice(arr [][]float64) {
	
	filname := fmt.Sprintf(OD()+"dwpos%06d.txt",saveIndex)
	ny,nz := Mesh().Size()[Y],Mesh().Size()[Z]
	// Find current working dir
	f,err := os.Create(filname) //https://gobyexample.com/writing-files
	check(err)
	// Write timestamp
	line := fmt.Sprintf("#time = %e\n",Time)
	f.WriteString(line)

	for i:=0; i<nz; i++ {
		for j:=0; j<ny; j++{
			line := fmt.Sprintf("%e\n",arr[i][j])
			f.WriteString(line)
		}
	}
	f.Sync()
	saveIndex++
}


func saveDwM(posm [][]float64, posind [][]int) {
	
	filname := fmt.Sprintf(OD()+"dwposm%06d.txt",saveIndex)
	ny,nz := Mesh().Size()[Y],Mesh().Size()[Z]
	// Find current working dir
	f,err := os.Create(filname) //https://gobyexample.com/writing-files
	check(err)
	// Write timestamp
	line := fmt.Sprintf("#time = %e\n",Time)
	f.WriteString(line)
	// Write column headers
	line = fmt.Sprintf("#dwpos (m), mx, my\n")
	f.WriteString(line)

	M := &M
	for i:=0; i<nz; i++ {
		for j:=0; j<ny; j++{
			minIndex := posind[i][j]
			mcell := M.GetCell(minIndex,j,i)
			line := fmt.Sprintf("%[1]e,%[2]e,%[3]e\n",posm[i][j],mcell[X],mcell[Y])
			f.WriteString(line)
		}
	}
	f.Sync()
	saveIndex++
}
