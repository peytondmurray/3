#ifndef _DWWIDTH_H_
#define _DWWIDTH_H_

extern "C" __device__ float indexMulSumAtanh(float *, int, int, int, int, int, int, int);
extern "C" __device__ float sumAtanh(float *, int, int, int, int, int, int, int);
extern "C" __device__ int zcAlongY(float *, int, int, int, int, int);
#endif

