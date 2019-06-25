#ifndef _DWWIDTH_H_
#define _DWWIDTH_H_

// Convenience function for the domain wall width calculation. Sum(i*atanh(arr_i))
extern "C" __device__ float atanhMulIndSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    int k = 0;
    for (int i=iLo; i<iHi; i++) {
        ret += float(k)*atanhf(arr[idx(i, iy, iz)]);
        k++;
    }
    return ret;
}

// Convenience function for the domain wall width calculation. Sum(atanh(arr_i)).
extern "C" __device__ float atanhSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    for (int i=iLo; i<iHi; i++) {
        ret += atanhf(arr[idx(i, iy, iz)]);
    }
    return ret;
}

// Find the zero crossing of Mz along the x-direction, searching in the direction of decreasing x. If no zero crossing
// is found, return -1.
extern "C" __device__ int zcAlongX(float *m, int iy, int iz, int Nx, int Ny, int Nz) {
    for (int ix=Nx-1; ix>1; ix--) {
        if (m[idx(ix, iy, iz)]*m[idx(ix-1, iy, iz)] < float(0)) return ix-1;
    }
    return -1;
}
#endif

