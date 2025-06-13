#include "common.cuh"

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