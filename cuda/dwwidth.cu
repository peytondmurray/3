#include "stencil.h"
#include "dwwidth.h"

extern "C" __global__ void avgDWWidth(float* __restrict__ s, float* __restrict__ mz, int halfWidth,
                                      int Nx, int Ny, int Nz) {

    int _ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (_ix != 0 || iy >= Ny || iz >= Nz)
    {
        return;
    }

    int I = idx2D(iy, iz);

    // Find zero crossing along this row
    int ix = zcAlongY(mz, iy, iz, Nx, Ny, Nz);

    int iLo = ix-halfWidth;
    int iHi = ix+halfWidth+1;

    float N = float(2*halfWidth+1);                                     //# of points
    float maxN = N-1;

    float t4 = maxN*(maxN+1)/2;
    float t1 = N*indexMulSumAtanh(mz, iLo, iHi, iy, iz, Nx, Ny, Nz)/t4; // inf most of the time
    float t2 = sumAtanh(mz, iLo, iHi, iy, iz, Nx, Ny, Nz, s);           //Inf a bunch of the time
    float t3 = N*(maxN*(maxN+1)*(2*maxN+1)/6.0)/t4;                     // N*Sum(n^2)

    s[I] = (t3-t4)/(t1-t2);                                             // Inverse of the fit coefficient tanh(ax)
                                                                        // gives the domain wall width.

    return;
}

extern "C" __device__ float indexMulSumAtanh(float *arr, int ixLo, int ixHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    int k = 0;
    for (int i=ixLo; i<ixHi; i++) {
        ret += float(k)*atanhf(arr[idx(i, iy, iz)]);
        k ++;
    }
    return ret;
}

extern "C" __device__ float sumAtanh(float *arr, int ixLo, int ixHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    for (int i=ixLo; i<ixHi; i++) {
        ret += atanhf(arr[idx(i, iy, iz)]);
    }
    return ret;
}

extern "C" __device__ int zcAlongY(float *m, int iy, int iz, int Nx, int Ny, int Nz) {
    for (int ix=Nx-1; ix>1; ix--) {
        if (m[idx(ix, iy, iz)]*m[idx(ix-1, iy, iz)] < float(0)) return ix-1;
    }
    return -1;
}