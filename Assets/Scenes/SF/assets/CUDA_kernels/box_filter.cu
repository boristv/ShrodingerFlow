//#include "common.cuh"

#define IDX(i,j,k,nx,ny)  (((k)*(ny) + (j))*(nx) + (i))

/* ───────────────────────────────────────────────────────── *
1. 3×3×3 box filter                                     *
 * ───────────────────────────────────────────────────────── */
extern "C" __global__
void box_filter(float *dst, const float *src,
                int nx, int ny, int nz)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
    if (i>=nx || j>=ny || k>=nz) return;
    
    float sum=0.f;
    #pragma unroll
    for (int dz=-1; dz<=1; ++dz)
        for (int dy=-1; dy<=1; ++dy)
            for (int dx=-1; dx<=1; ++dx)
            {
                int ii = min(max(i+dx,0), nx-1);
                int jj = min(max(j+dy,0), ny-1);
                int kk = min(max(k+dz,0), nz-1);
                sum += src[IDX(ii,jj,kk,nx,ny)];
            }
    dst[IDX(i,j,k,nx,ny)] = sum / 27.f;
}