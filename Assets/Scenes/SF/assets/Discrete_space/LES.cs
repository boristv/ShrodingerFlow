using System;
using ManagedCuda;
using ManagedCuda.VectorTypes;
using source.assets.CUDA_kernels;
using source.assets.Discrete_space;
using source.assets.Discrete_space.utils; // Velocity, SpaceProperties

namespace source.assets.Les
{
    /// <summary>
    /// Статический класс, реализующий Smagorinsky‑LES поверх квантового солвера.
    /// Интерфейс стабилизирован: 4 публичных метода, которые дергаются из
    ///     ISF.Init()  →  LES.Init(...)
    ///     LES_step() →  Filter → ComputeNuT → ApplySGS
    /// Внутри держим только буфер ν_t; velocity буферы приходят извне.
    /// </summary>
    public static class LES
    {
        /* ----------------  device buffers  ---------------- */
        private static CudaDeviceVariable<float> nu_t;      // SGS‑вязкость

        /* ------------------  CUDA ядра  ------------------- */
        public static CudaKernel boxFilter;
        public static CudaKernel computeStrain;
        public static CudaKernel addTurbViscosity;

        /* сохранённые размеры блока/сетки (общие для всех ядёр) */
        private static dim3 block, grid;

        /* ============================================================= */
        #region  Init
        /// <summary>
        /// Вызывается однократно из ISF.Init().
        /// Загружает ядра, настраивает grid/block и выделяет ν_t.
        /// </summary>
        public static void Init(SpaceProperties p)
        {
            // 1. Буфер SGS‑вязкости
            nu_t = new CudaDeviceVariable<float>(p.num);

            // 2. Загрузка PTX / CUDA ядёр
            boxFilter        = KernelLoader.load_kernel("boxFilter");
            computeStrain    = KernelLoader.load_kernel("computeStrain");
            addTurbViscosity = KernelLoader.load_kernel("addTurbViscosity");

            // 3. Настройка сетки
            block = new dim3(8, 8, 8);
            grid  = new dim3((p.resx + 7) / 8,
                             (p.resy + 7) / 8,
                             (p.resz + 7) / 8);

            // — boxFilter использует shared‑memory (10³ float)
            boxFilter.BlockDimensions      = block;
            boxFilter.GridDimensions       = grid;
            boxFilter.DynamicSharedMemory  = (uint)((8 + 2) * (8 + 2) * (8 + 2) * sizeof(float));

            computeStrain.BlockDimensions   = block;
            computeStrain.GridDimensions    = grid;

            addTurbViscosity.BlockDimensions = block;
            addTurbViscosity.GridDimensions  = grid;
        }
        #endregion

        /* ============================================================= */
        #region  Public API
        /// <summary>
        /// 3×3×3 box‑filter: src → dst.  Вызываем по компонентам.
        /// </summary>
        public static void Filter(Velocity src, Velocity dst, SpaceProperties p)
        {
            boxFilter.Run(dst.vx.DevicePointer, src.vx.DevicePointer,
                          p.resx, p.resy, p.resz);
            boxFilter.Run(dst.vy.DevicePointer, src.vy.DevicePointer,
                          p.resx, p.resy, p.resz);
            boxFilter.Run(dst.vz.DevicePointer, src.vz.DevicePointer,
                          p.resx, p.resy, p.resz);
        }

        /// <summary>
        /// Вычисляем Smagorinsky ν_t = (C_s Δ)^2 |S|  на всём объёме.
        /// </summary>
        public static void ComputeNuT(Velocity filtered, SpaceProperties p)
        {
            float delta = p.filterFac * MathF.Min(p.dx, MathF.Min(p.dy, p.dz));

            computeStrain.Run(nu_t.DevicePointer,
                              filtered.vx.DevicePointer,
                              filtered.vy.DevicePointer,
                              filtered.vz.DevicePointer,
                              p.dx, p.dy, p.dz,
                              p.Cs, delta,
                              p.resx, p.resy, p.resz);
        }

        /// <summary>
        /// Добавляет дивергенцию SGS‑стресса к полю скоростей.
        /// </summary>
        public static void ApplySGS(Velocity vel, SpaceProperties p)
        {
            addTurbViscosity.Run(vel.vx.DevicePointer,
                                 vel.vy.DevicePointer,
                                 vel.vz.DevicePointer,
                                 nu_t.DevicePointer,
                                 p.dx, p.dy, p.dz,
                                 p.resx, p.resy, p.resz,
                                 p.dt);
        }
        #endregion
    }
}