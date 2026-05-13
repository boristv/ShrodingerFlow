using ComputeShaderSF;
using UnityEngine;

namespace ShrodingerFlow.AnalyticalTests
{
    /// <summary>
    /// Один «снимок» физических диагностик на текущем шаге.
    /// Все интегралы умножены на dV (для сетки 64³, sizeX=4 — это (4/64)³ ≈ 2.4e-4).
    /// </summary>
    public struct IsfDiagnosticsSnapshot
    {
        /// <summary>E_kin = 0.5 · Σ |u|² · dV.</summary>
        public float KineticEnergy;
        /// <summary>Σ ρ · dV, где ρ = |ψ₁|² + |ψ₂|². После Normalize ≈ sizeX·sizeY·sizeZ.</summary>
        public float DensityIntegral;
        /// <summary>max |ω| на сетке (центральные разности с периодикой).</summary>
        public float MaxVorticity;
        /// <summary>RMS(|ω|).</summary>
        public float RmsVorticity;
        /// <summary>RMS(|u|).</summary>
        public float RmsVelocity;
        /// <summary>max |u|.</summary>
        public float MaxVelocity;
        /// <summary>Σ |ω|² · dV — энстрофия (≈ интегральная мера активности вихря).</summary>
        public float Enstrophy;
        /// <summary>Центроид «вихревого облака», взвешенный по |ω|² (для трекинга кольца).</summary>
        public Vector3 VorticityCentroid;
    }

    /// <summary>
    /// Консервативные характеристики ISF: гидродинамический импульс, спиральность, дивергенция.
    /// Вычисляются через единственный GetData() per буфер скорости.
    /// </summary>
    public struct IsfConservationSnapshot
    {
        /// <summary>P = ½ ∫ x × ω dV — гидродинамический импульс. Сохраняется при свободной эволюции.</summary>
        public Vector3 Impulse;
        /// <summary>|P| — модуль импульса.</summary>
        public float ImpulseMag;
        /// <summary>H = ∫ v·ω dV — спиральность. Сохраняется при ν=0 (идеальный ISF).</summary>
        public float Helicity;
        /// <summary>RMS(∇·v) — мера нарушения несжимаемости. Должна быть ≈ 0 после pressure projection.</summary>
        public float DivergenceRms;
    }

    /// <summary>
    /// Утилиты считывают буферы CSISF/CSVelocity на CPU и считают интегральные характеристики.
    /// Это медленный путь — он осознанный, тесты не критичны к FPS, важна цифровая точность.
    /// </summary>
    public static class IsfDiagnostics
    {
        /// <summary>
        /// Вычисляет гидродинамический импульс P, спиральность H и RMS дивергенции в одном проходе.
        /// Все производные — центральные разности с периодическими ГУ.
        /// </summary>
        public static IsfConservationSnapshot SampleConservation(CSISF isf, CSVelocity vel)
        {
            int rx = isf.resX, ry = isf.resY, rz = isf.resZ;
            int num = rx * ry * rz;
            float dx = isf.dx, dy = isf.dy, dz = isf.dz;
            float dV = dx * dy * dz;
            float twoDx = 2f * dx, twoDy = 2f * dy, twoDz = 2f * dz;

            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            vel.vx.GetData(vx);
            vel.vy.GetData(vy);
            vel.vz.GetData(vz);

            double impX = 0, impY = 0, impZ = 0;
            double helicity = 0;
            double divSqSum = 0;

            for (int i = 0; i < rx; i++)
            {
                int ip = (i + 1) % rx, im = (i - 1 + rx) % rx;
                for (int j = 0; j < ry; j++)
                {
                    int jp = (j + 1) % ry, jm = (j - 1 + ry) % ry;
                    for (int k = 0; k < rz; k++)
                    {
                        int kp = (k + 1) % rz, km = (k - 1 + rz) % rz;
                        int idx    = i  * ry * rz + j  * rz + k;
                        int idx_ip = ip * ry * rz + j  * rz + k;
                        int idx_im = im * ry * rz + j  * rz + k;
                        int idx_jp = i  * ry * rz + jp * rz + k;
                        int idx_jm = i  * ry * rz + jm * rz + k;
                        int idx_kp = i  * ry * rz + j  * rz + kp;
                        int idx_km = i  * ry * rz + j  * rz + km;

                        // Завихрённость: ω = ∇ × v (центральные разности)
                        float dvz_dy = (vz[idx_jp] - vz[idx_jm]) / twoDy;
                        float dvy_dz = (vy[idx_kp] - vy[idx_km]) / twoDz;
                        float dvx_dz = (vx[idx_kp] - vx[idx_km]) / twoDz;
                        float dvz_dx = (vz[idx_ip] - vz[idx_im]) / twoDx;
                        float dvy_dx = (vy[idx_ip] - vy[idx_im]) / twoDx;
                        float dvx_dy = (vx[idx_jp] - vx[idx_jm]) / twoDy;

                        float wx = dvz_dy - dvy_dz;
                        float wy = dvx_dz - dvz_dx;
                        float wz = dvy_dx - dvx_dy;

                        float xp = isf.pxCPU[idx];
                        float yp = isf.pyCPU[idx];
                        float zp = isf.pzCPU[idx];

                        // Импульс: P = ½ ∫ x × ω dV
                        impX += 0.5 * (yp * wz - zp * wy) * dV;
                        impY += 0.5 * (zp * wx - xp * wz) * dV;
                        impZ += 0.5 * (xp * wy - yp * wx) * dV;

                        // Спиральность: H = ∫ v·ω dV
                        helicity += (vx[idx] * wx + vy[idx] * wy + vz[idx] * wz) * dV;

                        // Дивергенция: ∇·v
                        float divV = (vx[idx_ip] - vx[idx_im]) / twoDx
                                   + (vy[idx_jp] - vy[idx_jm]) / twoDy
                                   + (vz[idx_kp] - vz[idx_km]) / twoDz;
                        divSqSum += divV * divV;
                    }
                }
            }

            var impulse = new Vector3((float)impX, (float)impY, (float)impZ);
            return new IsfConservationSnapshot
            {
                Impulse      = impulse,
                ImpulseMag   = impulse.magnitude,
                Helicity     = (float)helicity,
                DivergenceRms = (float)System.Math.Sqrt(divSqSum / num)
            };
        }

        /// <summary>
        /// Читает буферы скорости на CPU и возвращает массив |v|² = vx²+vy²+vz² для каждой ячейки.
        /// Используется тестом оптической согласованности вместо ρ=|ψ|², которая однородна в ISF.
        /// </summary>
        public static float[] SampleVelocitySquared(CSVelocity vel, int num)
        {
            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            vel.vx.GetData(vx);
            vel.vy.GetData(vy);
            vel.vz.GetData(vz);

            var vel2 = new float[num];
            for (int i = 0; i < num; i++)
                vel2[i] = vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
            return vel2;
        }

        /// <summary>Запросить «полный» снимок: единственный GetData() per buffer, всё остальное — на CPU.</summary>
        public static IsfDiagnosticsSnapshot SampleAll(CSISF isf, CSVelocity vel)
        {
            int rx = isf.resX, ry = isf.resY, rz = isf.resZ;
            int num = rx * ry * rz;

            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            vel.vx.GetData(vx);
            vel.vy.GetData(vy);
            vel.vz.GetData(vz);

            var psi1 = new Vector2[num];
            var psi2 = new Vector2[num];
            isf.psi1.GetData(psi1);
            isf.psi2.GetData(psi2);

            float dx = isf.dx, dy = isf.dy, dz = isf.dz;
            float dV = dx * dy * dz;

            double kineticSum = 0.0;
            double densitySum = 0.0;
            double velSqSum = 0.0;
            float maxVel2 = 0f;

            for (int i = 0; i < num; i++)
            {
                float u = vx[i], v = vy[i], w = vz[i];
                float v2 = u * u + v * v + w * w;
                velSqSum += v2;
                kineticSum += v2;
                if (v2 > maxVel2) maxVel2 = v2;

                float p1 = psi1[i].x * psi1[i].x + psi1[i].y * psi1[i].y;
                float p2 = psi2[i].x * psi2[i].x + psi2[i].y * psi2[i].y;
                densitySum += p1 + p2;
            }

            double enstrophySum = 0.0;
            double omegaSqSum = 0.0;
            float maxOmega2 = 0f;

            double cX = 0.0, cY = 0.0, cZ = 0.0, wSum = 0.0;

            float twoDx = 2f * dx;
            float twoDy = 2f * dy;
            float twoDz = 2f * dz;

            for (int i = 0; i < rx; i++)
            {
                int ip = (i + 1) % rx;
                int im = (i - 1 + rx) % rx;
                for (int j = 0; j < ry; j++)
                {
                    int jp = (j + 1) % ry;
                    int jm = (j - 1 + ry) % ry;
                    for (int k = 0; k < rz; k++)
                    {
                        int kp = (k + 1) % rz;
                        int km = (k - 1 + rz) % rz;

                        int idx = i * ry * rz + j * rz + k;

                        int idx_ip = ip * ry * rz + j * rz + k;
                        int idx_im = im * ry * rz + j * rz + k;
                        int idx_jp = i * ry * rz + jp * rz + k;
                        int idx_jm = i * ry * rz + jm * rz + k;
                        int idx_kp = i * ry * rz + j * rz + kp;
                        int idx_km = i * ry * rz + j * rz + km;

                        float dvz_dy = (vz[idx_jp] - vz[idx_jm]) / twoDy;
                        float dvy_dz = (vy[idx_kp] - vy[idx_km]) / twoDz;
                        float dvx_dz = (vx[idx_kp] - vx[idx_km]) / twoDz;
                        float dvz_dx = (vz[idx_ip] - vz[idx_im]) / twoDx;
                        float dvy_dx = (vy[idx_ip] - vy[idx_im]) / twoDx;
                        float dvx_dy = (vx[idx_jp] - vx[idx_jm]) / twoDy;

                        float omegaX = dvz_dy - dvy_dz;
                        float omegaY = dvx_dz - dvz_dx;
                        float omegaZ = dvy_dx - dvx_dy;

                        float w2 = omegaX * omegaX + omegaY * omegaY + omegaZ * omegaZ;
                        omegaSqSum += w2;
                        enstrophySum += w2;
                        if (w2 > maxOmega2) maxOmega2 = w2;

                        double weight = w2;
                        wSum += weight;
                        cX += weight * isf.pxCPU[idx];
                        cY += weight * isf.pyCPU[idx];
                        cZ += weight * isf.pzCPU[idx];
                    }
                }
            }

            var s = new IsfDiagnosticsSnapshot
            {
                KineticEnergy = (float)(0.5 * kineticSum * dV),
                DensityIntegral = (float)(densitySum * dV),
                MaxVorticity = Mathf.Sqrt(maxOmega2),
                RmsVorticity = (float)System.Math.Sqrt(omegaSqSum / num),
                RmsVelocity = (float)System.Math.Sqrt(velSqSum / num),
                MaxVelocity = Mathf.Sqrt(maxVel2),
                Enstrophy = (float)(enstrophySum * dV),
                VorticityCentroid = wSum > 0.0
                    ? new Vector3((float)(cX / wSum), (float)(cY / wSum), (float)(cZ / wSum))
                    : Vector3.zero
            };
            return s;
        }

        /// <summary>
        /// Рендерит срез |ω| (двумерное сечение через указанную ось/индекс) в Texture2D с применением viridis-подобной палитры.
        /// Возвращаемую текстуру можно сохранить через <see cref="TestArtifactWriter.Png"/>.
        /// </summary>
        public static Texture2D RenderVorticitySlice(CSISF isf, CSVelocity vel, SliceAxis axis, int sliceIndex)
        {
            int rx = isf.resX, ry = isf.resY, rz = isf.resZ;
            int num = rx * ry * rz;

            var vx = new float[num];
            var vy = new float[num];
            var vz = new float[num];
            vel.vx.GetData(vx);
            vel.vy.GetData(vy);
            vel.vz.GetData(vz);

            float twoDx = 2f * isf.dx;
            float twoDy = 2f * isf.dy;
            float twoDz = 2f * isf.dz;

            float OmegaMag(int i, int j, int k)
            {
                int ip = (i + 1) % rx, im = (i - 1 + rx) % rx;
                int jp = (j + 1) % ry, jm = (j - 1 + ry) % ry;
                int kp = (k + 1) % rz, km = (k - 1 + rz) % rz;

                int idx_ip = ip * ry * rz + j * rz + k;
                int idx_im = im * ry * rz + j * rz + k;
                int idx_jp = i * ry * rz + jp * rz + k;
                int idx_jm = i * ry * rz + jm * rz + k;
                int idx_kp = i * ry * rz + j * rz + kp;
                int idx_km = i * ry * rz + j * rz + km;

                float dvz_dy = (vz[idx_jp] - vz[idx_jm]) / twoDy;
                float dvy_dz = (vy[idx_kp] - vy[idx_km]) / twoDz;
                float dvx_dz = (vx[idx_kp] - vx[idx_km]) / twoDz;
                float dvz_dx = (vz[idx_ip] - vz[idx_im]) / twoDx;
                float dvy_dx = (vy[idx_ip] - vy[idx_im]) / twoDx;
                float dvx_dy = (vx[idx_jp] - vx[idx_jm]) / twoDy;

                float wx = dvz_dy - dvy_dz;
                float wy = dvx_dz - dvz_dx;
                float wz = dvy_dx - dvx_dy;
                return Mathf.Sqrt(wx * wx + wy * wy + wz * wz);
            }

            int w, h;
            switch (axis)
            {
                case SliceAxis.X: w = ry; h = rz; break;
                case SliceAxis.Y: w = rx; h = rz; break;
                default: w = rx; h = ry; break;
            }

            var values = new float[w * h];
            float maxV = 0f;
            for (int a = 0; a < w; a++)
            {
                for (int b = 0; b < h; b++)
                {
                    float v;
                    switch (axis)
                    {
                        case SliceAxis.X:
                            v = OmegaMag(Mathf.Clamp(sliceIndex, 0, rx - 1), a, b);
                            break;
                        case SliceAxis.Y:
                            v = OmegaMag(a, Mathf.Clamp(sliceIndex, 0, ry - 1), b);
                            break;
                        default:
                            v = OmegaMag(a, b, Mathf.Clamp(sliceIndex, 0, rz - 1));
                            break;
                    }
                    values[a + b * w] = v;
                    if (v > maxV) maxV = v;
                }
            }

            var tex = new Texture2D(w, h, TextureFormat.RGB24, false);
            var pixels = new Color32[w * h];
            float inv = maxV > 1e-12f ? 1f / maxV : 0f;
            for (int p = 0; p < pixels.Length; p++)
            {
                float t = values[p] * inv;
                pixels[p] = Colormap.Viridis(t);
            }
            tex.SetPixels32(pixels);
            tex.Apply(false);
            return tex;
        }
    }

    public enum SliceAxis { X, Y, Z }

    /// <summary>Простая аппроксимация viridis-палитры. Для дополнительной выразительности срезов.</summary>
    internal static class Colormap
    {
        private static readonly Vector3[] _viridisAnchors =
        {
            new Vector3(0.267f, 0.005f, 0.329f),
            new Vector3(0.282f, 0.140f, 0.458f),
            new Vector3(0.254f, 0.265f, 0.530f),
            new Vector3(0.207f, 0.372f, 0.553f),
            new Vector3(0.164f, 0.471f, 0.558f),
            new Vector3(0.128f, 0.567f, 0.551f),
            new Vector3(0.135f, 0.659f, 0.518f),
            new Vector3(0.267f, 0.749f, 0.441f),
            new Vector3(0.478f, 0.821f, 0.318f),
            new Vector3(0.741f, 0.873f, 0.150f),
            new Vector3(0.993f, 0.906f, 0.144f),
        };

        public static Color32 Viridis(float t)
        {
            t = Mathf.Clamp01(t);
            float scaled = t * (_viridisAnchors.Length - 1);
            int i = Mathf.Clamp(Mathf.FloorToInt(scaled), 0, _viridisAnchors.Length - 2);
            float f = scaled - i;
            Vector3 a = _viridisAnchors[i];
            Vector3 b = _viridisAnchors[i + 1];
            Vector3 c = Vector3.Lerp(a, b, f);
            return new Color32(
                (byte)Mathf.Clamp(Mathf.RoundToInt(c.x * 255f), 0, 255),
                (byte)Mathf.Clamp(Mathf.RoundToInt(c.y * 255f), 0, 255),
                (byte)Mathf.Clamp(Mathf.RoundToInt(c.z * 255f), 0, 255),
                255);
        }
    }
}
