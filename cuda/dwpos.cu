#include "stencil.h"
#include "domainwall.h"

/*
Get the indices of the zero crossing of the Mz component, right near domain wall.
Note the return value is a float, and must be converted to int before being used elsewhere.
*/
extern "C" __global__ void dwPos(float* __restrict__ pos, float* __restrict__ mz, int Nx, int Ny, int Nz) {

    int _ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (_ix > 0 || iy >= Ny || iz >= Nz) return;

    // Find zero crossing along this row
    pos[idx2D(iy, iz)] = float(zcAlongX(mz, iy, iz, Nx, Ny, Nz));


    return;
}
