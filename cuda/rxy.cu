#include "stencil.h"

extern "C" __global__ void
setrxy(float* __restrict__ s,
       float* __restrict__ mx, float* __restrict__ my,float* __restrict__ mz,
       int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= Nx || iy >= Ny || iz >= Nz)
    {
        return;
    }

    int I = idx(ix, iy, iz);                      // central cell index
    s[I] = sqrtf(mx[I]*mx[I] + my[I]*my[I]);
}