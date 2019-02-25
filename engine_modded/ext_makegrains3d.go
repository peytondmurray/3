package engine

import (
	"math"
	"math/rand"
	//"fmt"
)

func init() {
	DeclFunc("ext_makegrains3d", Voronoi3d, "Voronoi tesselation (grain size, num regions)")
	DeclFunc("ext_setPBCvoronoi3d", SetPBCVoronoi3d, "Periodic Voronoi tesselation direction (X,Y,Z), 0 or 1")
}

func Voronoi3d(grainsize float64, numRegions, seed int) {
	SetBusy(true)
	defer SetBusy(false)
	
	t := newTesselation3d(grainsize, numRegions, int64(seed))
	regions.hist = append(regions.hist, t.RegionOf3d)
	regions.render(t.RegionOf3d)
}

func SetPBCVoronoi3d(bx,by,bz int) {
	
	arg("SetPBCVoronoi3d x", bx != 0 || bx !=1)
	arg("SetPBCVoronoi3d y", by != 0 || by !=1)
	arg("SetPBCVoronoi3d z", bz != 0 || bz !=1)
	//fmt.Println(bx,by,bz)
	vPBCx = bx
	vPBCy = by
	vPBCz = bz
}

type tesselation3d struct {
	grainsize float64
	tilesize  float64
	maxRegion int
	cache3d     map[int3][]center3d
	seed      int64
	rnd       *rand.Rand
}

var (
	vPBCx = 0
	vPBCy = 0
	vPBCz = 0
)

// integer tile coordinate
type int3 struct{ x, y, z int }

// Voronoi center info
type center3d struct {
	x, y, z   float64 // center position (m)
	region byte    // region for all cells near center
}

// nRegion exclusive
func newTesselation3d(grainsize float64, nRegion int, seed int64) *tesselation3d {
	return &tesselation3d{grainsize,
		float64(float32(grainsize * TILE3)), // expect 4 grains/block, 36 per 3x3 blocks = safe, relatively round number
		nRegion,
		make(map[int3][]center3d),
		seed,
		rand.New(rand.NewSource(0))}
}

const (
	TILE3   = 2           // tile size in grains
	LAMBDA3 = TILE3 * TILE3 * TILE3 // expected grains per tile
)

// Returns the region of the grain where cell at x,y,z belongs to
func (t *tesselation3d) RegionOf3d(x, y, z float64) int {
	tile := t.tileOf3d(x, y, z) // tile containing x,y

	// look for nearest center in tile + neighbors
	nearest := center3d{x, y, z, 0} // dummy initial value, but safe should the infinite impossibility strike.
	mindist := math.Inf(1)
	ws := Mesh().WorldSize() // globalmesh_.WorldSize()
	// Tile coordinates limiting the system used for filtering.
	xU := t.tileOf3d(ws[0]/2,0,0).x
	xL := t.tileOf3d(-ws[0]/2,0,0).x
	yU := t.tileOf3d(0,ws[1]/2,0).y
	yL := t.tileOf3d(0,-ws[1]/2,0).y
	zU := t.tileOf3d(0,0,ws[2]/2).z
	zL := t.tileOf3d(0,0,-ws[2]/2).z
	for tx := tile.x - 1; tx <= tile.x+1; tx++ {
		txPBC := tx
		if vPBCx == 1 {
			if tx > xU {
				txPBC = xL
			}
			if tx < xL {
				txPBC = xU
			}
		}
		for ty := tile.y - 1; ty <= tile.y+1; ty++ {
			tyPBC := ty
			if vPBCy == 1 {
				if ty > yU {
					tyPBC = yL
				}
				if ty < yL {
					tyPBC = yU
				}
			}
			for tz := tile.z - 1; tz <= tile.z+1; tz++ {
				tzPBC := tz
				if vPBCz == 1 {
					if tz > zU {
						tzPBC = zL
					}
					if tz < zL {
						tzPBC = zU
				    }
				}
				centers := t.centersInTile3d(txPBC, tyPBC, tzPBC)
				for _, c := range centers {
				    dx := x-c.x
				    dy := y-c.y
				    dz := z-c.z
				    if vPBCx == 1 {
						if dx > ws[0] * 0.5 {
							dx = dx - ws[0]
						}
						if dx <= -ws[0] * 0.5 {
							dx = dx + ws[0]
						}
					}
				    if vPBCy == 1 {
						if dy > ws[1] * 0.5 {
							dy = dy - ws[1]
						}
						if dy <= -ws[1] * 0.5 {
							dy = dy + ws[1]
						}
					}
				    if vPBCz == 1 {
						if dz > ws[2] * 0.5 {
							dz = dz - ws[2]
						}
						if dz <= -ws[1] * 0.5 {
							dz = dz + ws[2]
						}
					}
					dist := sqr(dx) + sqr(dy) + sqr(dz)
					if dist < mindist {
						nearest = c
						mindist = dist
					}
				}
			}
		}
	}
	//fmt.Println("nearest", x, y, ":", nearest)
	return int(nearest.region)
}

// Returns the list of Voronoi centers in tile(ix, iy), using only ix,iy to seed the random generator
func (t *tesselation3d) centersInTile3d(tx, ty, tz int) []center3d {
	pos := int3{tx, ty, tz}
	if c, ok := t.cache3d[pos]; ok {
		return c
	} else {
		// tile-specific seed that works for positive and negative tx, ty
		seed := (int64(tz)+(1<<24))*(1<<24)*(1<<24) + (int64(ty)+(1<<24))*(1<<24) + (int64(tx) + (1 << 24))
		t.rnd.Seed(seed ^ t.seed)
		N := t.poisson3d(LAMBDA3)
		c := make([]center3d, N)

		// absolute position of tile (m)
		x0, y0, z0 := float64(tx)*t.tilesize, float64(ty)*t.tilesize, float64(tz)*t.tilesize

		for i := range c {
			// random position inside tile
			c[i].x = x0 + t.rnd.Float64()*t.tilesize
			c[i].y = y0 + t.rnd.Float64()*t.tilesize
			c[i].z = z0 + t.rnd.Float64()*t.tilesize
			c[i].region = byte(t.rnd.Intn(t.maxRegion))
		}
		t.cache3d[pos] = c
		return c
	}
}

//func sqr(x float64) float64 { return x * x }

func (t *tesselation3d) tileOf3d(x, y, z float64) int3 {
	ix := int(math.Floor(x / t.tilesize))
	iy := int(math.Floor(y / t.tilesize))
	iz := int(math.Floor(z / t.tilesize))
	return int3{ix, iy, iz}
}

// Generate poisson distributed numbers (according to Knuth)

func (t *tesselation3d) poisson3d(lambda float64) int {
	L := math.Exp(-lambda)
	k := 1
	p := t.rnd.Float64()
	for p > L {
		k++
		p *= t.rnd.Float64()
	}
	return k - 1
}
