#ifndef _SUBSET_H_
#define _SUBSET_H_

extern "C" __device__ int idxSubset(int ix, int iy, int iz, int NxSubset, int Ny, int Nz) {
    return ((iz*Ny) + iy)*NxSubset + ix;
}

#endif

