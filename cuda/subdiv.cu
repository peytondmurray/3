#include "stencil.h"

extern "C" __global__ void
setangularvel(float* __restrict__ sx, float* __restrict__ sy, float* __restrict__ sz,
              float* __restrict__ ax, float* __restrict__ ay,float* __restrict__ az,
              float* __restrict__ bx, float* __restrict__ by, float* __restrict__ bz,
              float dt, int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= Nx || iy >= Ny || iz >= Nz)
    {
        return;
    }

    int I = idx(ix, iy, iz);                      // central cell index
    sx[I] = (ax[I] - bx[I])/dt;
    sy[I] = (ay[I] - by[I])/dt;
    sz[I] = (az[I] - bz[I])/dt;
    return;
}