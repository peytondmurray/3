#ifndef _DWWIDTH_H_
#define _DWWIDTH_H_

extern "C" __device__ float atanhMulIndSum(float *, int, int, int, int, int, int, int);
extern "C" __device__ float atanhSum(float *, int, int, int, int, int, int, int);

// Find the zero crossing of Mz along the x-direction, searching in the direction of decreasing x. If no zero crossing
// is found, return -1.
extern "C" __device__ int zcAlongX(float *m, int iy, int iz, int Nx, int Ny, int Nz) {
    for (int ix=Nx-1; ix>1; ix--) {
        if (m[idx(ix, iy, iz)]*m[idx(ix-1, iy, iz)] < float(0)) return ix-1;
    }
    return -1;
}
#endif

