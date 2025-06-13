#include "cuda_runtime.h"

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

/* ───────────────────────────────────────────────────────── *
   2. Compute |S| and ν_t                                   *
 * ───────────────────────────────────────────────────────── */
extern "C" __global__
void compute_strain(float *nuT,
                    const float *u, const float *v, const float *w,
                    float dx, float dy, float dz,
                    float Cs,  float delta,
                    int nx, int ny, int nz)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
    if (i>=nx || j>=ny || k>=nz) return;

    // centred diffs (clamped)
    auto idx=[&](int ii,int jj,int kk){return IDX(min(max(ii,0),nx-1), min(max(jj,0),ny-1), min(max(kk,0),nz-1), nx, ny);} ;

    float du_dx = (u[idx(i+1,j,k)] - u[idx(i-1,j,k)])/(2.f*dx);
    float du_dy = (u[idx(i,j+1,k)] - u[idx(i,j-1,k)])/(2.f*dy);
    float du_dz = (u[idx(i,j,k+1)] - u[idx(i,j,k-1)])/(2.f*dz);

    float dv_dx = (v[idx(i+1,j,k)] - v[idx(i-1,j,k)])/(2.f*dx);
    float dv_dy = (v[idx(i,j+1,k)] - v[idx(i,j-1,k)])/(2.f*dy);
    float dv_dz = (v[idx(i,j,k+1)] - v[idx(i,j,k-1)])/(2.f*dz);

    float dw_dx = (w[idx(i+1,j,k)] - w[idx(i-1,j,k)])/(2.f*dx);
    float dw_dy = (w[idx(i,j+1,k)] - w[idx(i,j-1,k)])/(2.f*dy);
    float dw_dz = (w[idx(i,j,k+1)] - w[idx(i,j,k-1)])/(2.f*dz);

    float S11 = du_dx;
    float S22 = dv_dy;
    float S33 = dw_dz;
    float S12 = 0.5f*(du_dy + dv_dx);
    float S13 = 0.5f*(du_dz + dw_dx);
    float S23 = 0.5f*(dv_dz + dw_dy);

    float magS = sqrtf(2.f*(S11*S11 + S22*S22 + S33*S33)
                       + 4.f*(S12*S12 + S13*S13 + S23*S23));

    float delta2 = Cs*delta;
    delta2 *= delta2;
    nuT[IDX(i,j,k,nx,ny)] = delta2 * magS;
}

/* ───────────────────────────────────────────────────────── *
   3. Add SGS stress divergence                              *
 * ───────────────────────────────────────────────────────── */
extern "C" __global__
void add_turb_viscosity(float *u, float *v, float *w,
                        const float *nuT,
                        float dx,float dy,float dz,
                        int nx,int ny,int nz,
                        float dt)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
    if (i>=nx || j>=ny || k>=nz) return;

    auto idx=[&](int ii,int jj,int kk){return IDX(min(max(ii,0),nx-1), min(max(jj,0),ny-1), min(max(kk,0),nz-1), nx, ny);} ;

    // compute ν_t gradient (6‑point Laplacian‑like) & apply to velocity
    float nuE = nuT[idx(i+1,j,k)];
    float nuW = nuT[idx(i-1,j,k)];
    float nuN = nuT[idx(i,j+1,k)];
    float nuS = nuT[idx(i,j-1,k)];
    float nuTz= nuT[idx(i,j,k+1)];
    float nuB = nuT[idx(i,j,k-1)];
    float nuC = nuT[idx(i,j,k)  ];

    float lap_nu = (nuE + nuW - 2.f*nuC)/(dx*dx)
                 + (nuN + nuS - 2.f*nuC)/(dy*dy)
                 + (nuTz+ nuB - 2.f*nuC)/(dz*dz);

    u[IDX(i,j,k,nx,ny)] += dt * lap_nu;
    v[IDX(i,j,k,nx,ny)] += dt * lap_nu;
    w[IDX(i,j,k,nx,ny)] += dt * lap_nu;
}