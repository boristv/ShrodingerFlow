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
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int k = blockIdx.z * blockDim.z + threadIdx.z;
    if (i >= nx || j >= ny || k >= nz) return;

    // νt в шести соседних ячейках + центр
    float nuE = nuT[safeIdx(i+1, j  , k  , nx, ny, nz)];
    float nuW = nuT[safeIdx(i-1, j  , k  , nx, ny, nz)];
    float nuN = nuT[safeIdx(i  , j+1, k  , nx, ny, nz)];
    float nuS = nuT[safeIdx(i  , j-1, k  , nx, ny, nz)];
    float nuTz= nuT[safeIdx(i  , j  , k+1, nx, ny, nz)];
    float nuB = nuT[safeIdx(i  , j  , k-1, nx, ny, nz)];
    float nuC = nuT[IDX(i, j, k, nx, ny)];   // центр

    // «Лапласиан-подобный» вклад
    float lap_nu = (nuE + nuW - 2.0f * nuC) / (dx * dx)
                 + (nuN + nuS - 2.0f * nuC) / (dy * dy)
                 + (nuTz + nuB - 2.0f * nuC) / (dz * dz);

    int id = IDX(i, j, k, nx, ny);
    u[id] += dt * lap_nu;
    v[id] += dt * lap_nu;
    w[id] += dt * lap_nu;
}