#pragma once                       // защита от двойного включения
#include <cuda_runtime.h>          // dim3, cudaMalloc, cudaMemcpy, …

// Индекс в объёмном массиве nx×ny×nz (строчн. порядок)
#define IDX(i,j,k,nx,ny)  (((k)*(ny) + (j))*(nx) + (i))

// сюда же можно поместить inline-helpers
__device__ __forceinline__
int clampi(int v, int lo, int hi) { return max(lo, min(v, hi)); }