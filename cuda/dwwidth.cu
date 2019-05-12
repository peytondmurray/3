#include "stencil.h"
#include "dwwidth.h"

extern "C" __global__ void avgDWWidth(float* __restrict__ s, float* __restrict__ mz, int halfWidth,
                                      int Nx, int Ny, int Nz) {

    int _ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (_ix > 1 || iy >= Ny || iz >= Nz)
    {
        return;
    }

    // Find zero crossing along this row
    int ix = zcAlongY(mz, iy, iz, Nx, Ny, Nz);

    // Find the lower and upper bounds for the region near the zero crossing
    int iLo = ix-halfWidth;
    int iHi = ix+halfWidth+1;

    float N = float(2*halfWidth+1); //# of cells running across the domain wall

    float sumXi = (N-1)*N/2;
    float t1 = N*atanhMulIndSum(mz, iLo, iHi, iy, iz, Nx, Ny, Nz);
    float t2 = atanhSum(mz, iLo, iHi, iy, iz, Nx, Ny, Nz)*sumXi;
    float t3 = N*(N-1)*N*(2*N-1)/6.0;
    float t4 = sumXi*sumXi;


    s[idx2D(iy, iz)] = (t3-t4)/(t1-t2);
    // s[idx2D(iy, iz)] = (t1-t2)/(t3-t4);

    // s[idx2D(iy, iz)] = (t3-t4)/(t1-t2);     // Inverse of the fit coefficient tanh(ax)
    //                                         // gives the domain wall width (here, in units of the cellsize)

    return;
}

extern "C" __device__ float atanhMulIndSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    int k = 0;
    for (int i=iLo; i<iHi; i++) {
        ret += float(k)*atanhf(arr[idx(i, iy, iz)]);
        k++;
    }
    return ret;
}

extern "C" __device__ float atanhSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    for (int i=iLo; i<iHi; i++) {
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