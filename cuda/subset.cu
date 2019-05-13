#include "stencil.h"
#include "subset.h"

extern "C" __global__ void subsetAlongX(float* __restrict__ s, float* __restrict__ a, float* __restrict__ lowerBound,
                                        float* __restrict__ upperBound, int NxSubset, int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= Nx || iy >= Ny || iz >= Nz) return;
    if (ix < int(lowerBound[idx2D(iy, iz)]+0.5) || ix > int(upperBound[idx2D(iy, iz)])+0.5) return;

    s[idxSubset(ix-int(lowerBound[idx2D(iy, iz)]+0.5), iy, iz, NxSubset, Ny, Nz)] = a[idx(ix, iy, iz)];

    return;
}

extern "C" __device__ int idxSubset(int ix, int iy, int iz, int NxSubset, int Ny, int Nz) {
    return ((iz*Ny) + iy)*NxSubset + ix;
}