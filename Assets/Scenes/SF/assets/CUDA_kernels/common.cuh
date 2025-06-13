#pragma once
#include <cuda_runtime.h>

// Индекс в объёме (строчное хранение: x – fastest)
#define IDX(i,j,k,nx,ny)  (((k)*(ny) + (j))*(nx) + (i))

// Безопасное обращение с обрезкой по краям
__device__ __forceinline__ int safeIdx(int ii,int jj,int kk,
                                       int nx,int ny,int nz)
{
    ii = max(0, min(ii, nx-1));
    jj = max(0, min(jj, ny-1));
    kk = max(0, min(kk, nz-1));
    return IDX(ii,jj,kk,nx,ny);
}