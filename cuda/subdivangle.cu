#include "stencil.h"
#include "subdivangle.h"

// Don't use this, something the preprocessor is doing makes it different than what I have currently.
// #define PI 3.14159265358979323846264338327950288419716939937510582097494459

extern "C" __global__ void subDivAngle(float* __restrict__ s, float* __restrict__ a, float* __restrict__ b,
                                       int shift, float dt, int Nx, int Ny, int Nz) {

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    int iy = blockIdx.y * blockDim.y + threadIdx.y;
    int iz = blockIdx.z * blockDim.z + threadIdx.z;

    if (ix >= Nx || iy >= Ny || iz >= Nz) return;

    if (ix + shift < Nx && ix + shift >= 0) {
        s[idx(ix, iy, iz)] = deltaAngle(a[idx(ix, iy, iz)] - b[idx(ix+shift, iy, iz)])/dt;
    } else {
        s[idx(ix, iy, iz)] = 0;
    }

    return;
}

extern "C" __device__ float deltaAngle(float dA) {
    float pi = 3.14159265358979323846264338327950288419716939937510582097494459;
    if (dA < -1*pi) {
        return 2*pi + dA;
    } else if (dA > pi) {
        return 2*pi - dA;
    }
    return dA;
}
