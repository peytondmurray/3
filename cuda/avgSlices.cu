#include "stencil.h"

extern "C" __global__ void avgSlices(float* __restrict__ dst, float* __restrict__ src1, float* __restrict__ src2,
                                       int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= Nx || iy >= Ny || iz >= Nz) return;

    dst[idx(ix, iy, iz)] = 0.5*(src1[idx(ix, iy, iz)] + src2[idx(ix, iy, iz)]);

    return;
}
