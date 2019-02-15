package engine

import (
	"fmt"
	"math"
	"os"
	"testing"
	"time"
)

func TestFornberg(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	weights := fornbergWeights(x[len(x)-1], x, 1)
	for i := 0; i < len(weights); i++ {
		print("   ", weights[i])
		println("\n")
	}
}

func TestDifferentiation(t *testing.T) {

	file, err := os.Create("/tmp/benchmark.txt")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	nIter := 100000

	for k := 2; k < 50; k++ {
		var p posStack
		p.maxsize = k

		t0 := time.Now()
		// speed := float64(0)
		for i := 0; i < nIter; i++ {
			p.push(float64(i), math.Sin(float64(i)))
			if i > 1 {
				speed := p.speed()
				speed++
			}
		}

		file.WriteString(fmt.Sprintf("%d,%f\n", k, float64(nIter)/time.Since(t0).Seconds()))
	}
}
