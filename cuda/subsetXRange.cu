#include "stencil.h"
#include "subset.h"

/*
Extract a subset of a slice on the GPU which falls within some range of x-indices.
*/
extern "C" __global__ void subsetXRange(float* __restrict__ dst, float* __restrict__ src, int xmin, int xmax,
                                       int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= xmax-xmin || iy >= Ny || iz >= Nz) return;

    dst[idxSubset(ix, iy, iz, xmax-xmin, Ny, Nz)] = src[idx(ix+xmin, iy, iz)];

    return;
}
