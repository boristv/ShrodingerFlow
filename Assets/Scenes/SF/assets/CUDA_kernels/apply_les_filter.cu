extern "C" __global__
void apply_les_filter(cuFloatComplex* psi, int Nx, int Ny, int Nz,
                      float sizex, float sizey, float sizez, float k_cutoff)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    if (i >= Nx || j >= Ny || k >= Nz) return;

    int index = i * Ny * Nz + j * Nz + k;

    // Индексы в спектре: сдвигаем относительно центра
    int ix = i - Nx / 2;
    int iy = j - Ny / 2;
    int iz = k - Nz / 2;

    // Волновые числа
    float kx = ix / sizex;
    float ky = iy / sizey;
    float kz = iz / sizez;

    float k_mag2 = kx * kx + ky * ky + kz * kz;
    float filter = expf(-k_mag2 / (k_cutoff * k_cutoff));

    cuFloatComplex val = psi[index];
    psi[index].x = val.x * filter;
    psi[index].y = val.y * filter;
}
