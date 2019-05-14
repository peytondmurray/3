#include "stencil.h"
#include "dwwidth.h"
#include "subset.h"

/*
Calculate the domain wall width along each x-coordinate of the simulation region. Takes the zero-crossing
of the z-component of the magnetization as the center of the wall, and fits a line

    y = ax + b

to atanh(Mz) using linear least squares. Then the inverse of the slope, 1/a, gives the domain wall width.
This function writes the domain wall width along each x-coordinate to the array s, which should
be of shape (1, n[Y], n[Z]) on the Go side of things.

Note the zero-crossing of Mz is searched in the direction of decreasing x. Also, the domain wall width is
given in units of the cell size in the x-direction, which allows for some optimizations. Finally, the
ansatz that the magnetization Mz = tanh(ax) is not always valid, esp. when the domain wall is highly curved.

Overall, I found this method (including the averaging over s) made calculating the domain wall width ~10x
faster than on CPU.
*/
extern "C" __global__ void DWWidth(float* __restrict__ s, float* __restrict__ mz, int halfWidth,
                                      int Nx, int Ny, int Nz) {

    int _ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (_ix >= 1 || iy >= Ny || iz >= Nz)
    {
        return;
    }

    // Find zero crossing along this row
    int ix = zcAlongX(mz, iy, iz, Nx, Ny, Nz);

    // if (ix == -1) {
    //     s[idx2D(iy, iz)] = 0;
    //     return;
    // }

    // Find the lower and upper bounds for the region near the zero crossing
    int iLo = ix-halfWidth;
    int iHi = ix+halfWidth+1;

    float N = float(2*halfWidth+1); //# of cells running across the domain wall

    float sumXi = (N-1)*N/2;
    float t1 = N*atanhMulIndSum(mz, iLo, iHi, iy, iz, Nx, Ny, Nz);
    float t2 = atanhSum(mz, iLo, iHi, iy, iz, Nx, Ny, Nz)*sumXi;
    float t3 = N*(N-1)*N*(2*N-1)/6.0;
    float t4 = sumXi*sumXi;

    s[idx2D(iy, iz)] = fabsf((t3-t4)/(t1-t2)); // This gives the inverse of the least squares fit slope to tanh(mz), i.e. the domain wall width (in units of the cell size).

    return;
}

// Convenience function for the calculation above. Sums(i*atanh(arr_i))
extern "C" __device__ float atanhMulIndSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    int k = 0;
    for (int i=iLo; i<iHi; i++) {
        ret += float(k)*atanhf(arr[idx(i, iy, iz)]);
        k++;
    }
    return ret;
}

// Convenience function for the calculation above. Sums atanh(arr_i).
extern "C" __device__ float atanhSum(float *arr, int iLo, int iHi, int iy, int iz, int Nx, int Ny, int Nz) {
    float ret = 0;
    for (int i=iLo; i<iHi; i++) {
        ret += atanhf(arr[idx(i, iy, iz)]);
    }
    return ret;
}
