#include "stencil.h"
#include "subset.h"

extern "C" __global__ void nearDW(float* __restrict__ s, float* __restrict__ a, float* __restrict__ dwpos,
                                  int mw, int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= 2*mw+1 || iy >= Ny || iz >= Nz) return;

    s[idxSubset(ix, iy, iz, 2*mw+1, Ny, Nz)] = a[idx(int(dwpos[idx2D(iy, iz)] + ix - mw + 0.5), iy, iz)];
    // s[idxSubset(ix, iy, iz, 2*mw+1, Ny, Nz)] = int(dwpos[idx2D(iy, iz)] + ix - mw + 0.5); // Set s to be the index, for debugging purposes

    return;
}
