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

// Вычисляет |S| и SGS‑ν_t по Smagorinsky
extern "C" __global__ void computeStrain(float *nuT,
                              const float *u,const float *v,const float *w,
                              float dx,float dy,float dz,
                              float Cs,float Delta,
                              int nx,int ny,int nz)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z*blockDim.z + threadIdx.z;
    if(i>=nx || j>=ny || k>=nz) return;

    // центральные разности 1‑го порядка (для CUDA10.2 — экономим регистры)
    float du_dx = (u[safeIdx(i+1,j,k,nx,ny,nz)] - u[safeIdx(i-1,j,k,nx,ny,nz)])/(2*dx);
    float du_dy = (u[safeIdx(i,j+1,k,nx,ny,nz)] - u[safeIdx(i,j-1,k,nx,ny,nz)])/(2*dy);
    float du_dz = (u[safeIdx(i,j,k+1,nx,ny,nz)] - u[safeIdx(i,j,k-1,nx,ny,nz)])/(2*dz);

    float dv_dx = (v[safeIdx(i+1,j,k,nx,ny,nz)] - v[safeIdx(i-1,j,k,nx,ny,nz)])/(2*dx);
    float dv_dy = (v[safeIdx(i,j+1,k,nx,ny,nz)] - v[safeIdx(i,j-1,k,nx,ny,nz)])/(2*dy);
    float dv_dz = (v[safeIdx(i,j,k+1,nx,ny,nz)] - v[safeIdx(i,j,k-1,nx,ny,nz)])/(2*dz);

    float dw_dx = (w[safeIdx(i+1,j,k,nx,ny,nz)] - w[safeIdx(i-1,j,k,nx,ny,nz)])/(2*dx);
    float dw_dy = (w[safeIdx(i,j+1,k,nx,ny,nz)] - w[safeIdx(i,j-1,k,nx,ny,nz)])/(2*dy);
    float dw_dz = (w[safeIdx(i,j,k+1,nx,ny,nz)] - w[safeIdx(i,j,k-1,nx,ny,nz)])/(2*dz);

    // Симметризированный градиент (тензор деформации)
    float S11 = du_dx;
    float S22 = dv_dy;
    float S33 = dw_dz;

    float S12 = 0.5f*(du_dy + dv_dx);
    float S13 = 0.5f*(du_dz + dw_dx);
    float S23 = 0.5f*(dv_dz + dw_dy);

    float magS = sqrtf(2.0f*(S11*S11 + S22*S22 + S33*S33) +
                       4.0f*(S12*S12 + S13*S13 + S23*S23));

    float nu_t = (Cs * Delta)*(Cs * Delta) * magS;
    nuT[IDX(i,j,k,nx,ny)] = nu_t;
}