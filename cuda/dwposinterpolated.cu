#include "stencil.h"

extern "C" __global__ void dwPosInterpolated(float* __restrict__ dst, float* __restrict__ zc, float* __restrict__ mz,
                                       int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix > 1 || iy >= Ny || iz >= Nz) return;

    int zcl = int(zc[idx2D(iy, iz)] + 0.5);  // Index of point to left of zero-crossing in 3D slice
    int izcl = idx(zcl, iy, iz);   // CUDA index of point to left of zc
    int izcr = idx(zcl+1, iy, iz); // CUDA index of point to right of zc

    dst[idx2D(iy, iz)] = float(zcl) - (mz[izcl]/(mz[izcr]-mz[izcl]));

    return;
}
