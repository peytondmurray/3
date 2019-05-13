#include "stencil.h"
#include "dwwidth.h"

/*
Get the indices of the region near the zero crossing of the Mz component; the region
near the domain wall.

Note the return value is a float
*/
extern "C" __global__ void dwPos(float* __restrict__ pos, float* __restrict__ mz, int Nx, int Ny, int Nz) {

    int _ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (_ix > 0 || iy >= Ny || iz >= Nz) return;

    // Find zero crossing along this row
    int ix = zcAlongX(mz, iy, iz, Nx, Ny, Nz);

    // Find the lower and upper bounds for the region near the zero crossing
    // s[toIndex(0, iy, iz, 2, Ny, Nz)] = float(ix-mw);
    // s[toIndex(1, iy, iz, 2, Ny, Nz)] = float(ix+mw+1);
    pos[idx2D(iy, iz)] = float(ix);
    // lowerBound[idx2D(iy, iz)] = float(ix-mw);
    // upperBound[idx2D(iy, iz)] = float(ix+mw+1);

    return;
}
