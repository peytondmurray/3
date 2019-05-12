#ifndef _DWWIDTH_H_
#define _DWWIDTH_H_

extern "C" __device__ float atanhMulIndSum(float *, int, int, int, int, int, int, int);
extern "C" __device__ float atanhSum(float *, int, int, int, int, int, int, int);
extern "C" __device__ int zcAlongY(float *, int, int, int, int, int);
#endif

