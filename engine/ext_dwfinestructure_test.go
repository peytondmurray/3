package engine

import (
	"testing"
)

func TestFornberg(t *testing.T) {
	x := []float64{0, 1, 2, 3}
	weights := fornbergWeights(x[len(x)-1], x, 1)
	for i := 0; i < len(weights); i++ {
		print("   ", weights[i])
		println("\n")
	}
}
