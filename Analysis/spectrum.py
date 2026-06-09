#!/usr/bin/env python3
"""
Энергетический спектр E(k) поля скорости при разных ℏ — для раздела «проверяемые результаты».
Читает дампы из ../VelocityDumps/ (включи _dumpVelocity в SFHybridCS, прогони сцену
при нескольких значениях hbar), строит E(k) в log-log и сохраняет PNG.

Запуск:  python3 Analysis/spectrum.py
"""
import os, glob, struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP_DIR = os.path.join(HERE, "..", "VelocityDumps")
OUT = os.path.join(HERE, "spectrum_vs_hbar.png")


def load(path):
    with open(path, "rb") as f:
        rx, ry, rz = struct.unpack("<iii", f.read(12))
        dx, dy, dz, hbar = struct.unpack("<ffff", f.read(16))
        n = rx * ry * rz
        buf = np.frombuffer(f.read(n * 4 * 3), dtype="<f4")
    vx = buf[0:n].reshape((rx, ry, rz))
    vy = buf[n:2 * n].reshape((rx, ry, rz))
    vz = buf[2 * n:3 * n].reshape((rx, ry, rz))
    return dict(rx=rx, ry=ry, rz=rz, dx=dx, dy=dy, dz=dz, hbar=hbar, vx=vx, vy=vy, vz=vz)


def spectrum(d):
    """Радиально усреднённый энергетический спектр E(k) = 1/2 <|u_hat(k)|^2>."""
    rx, ry, rz = d["rx"], d["ry"], d["rz"]
    # 3D FFT каждой компоненты (нормировка 1/N)
    N = rx * ry * rz
    Ux = np.fft.fftn(d["vx"]) / N
    Uy = np.fft.fftn(d["vy"]) / N
    Uz = np.fft.fftn(d["vz"]) / N
    E = 0.5 * (np.abs(Ux) ** 2 + np.abs(Uy) ** 2 + np.abs(Uz) ** 2)
    # волновые числа (в единицах 1/длина), изотропная сетка по индексам
    kx = np.fft.fftfreq(rx, d=d["dx"]) * 2 * np.pi
    ky = np.fft.fftfreq(ry, d=d["dy"]) * 2 * np.pi
    kz = np.fft.fftfreq(rz, d=d["dz"]) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    Kmag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2).ravel()
    Ef = E.ravel()
    kmax = Kmag.max()
    nb = max(rx, ry, rz) // 2
    edges = np.linspace(0, kmax, nb + 1)
    idx = np.digitize(Kmag, edges) - 1
    Ek = np.zeros(nb)
    for b in range(nb):
        m = idx == b
        if m.any():
            Ek[b] = Ef[m].sum()
    kc = 0.5 * (edges[1:] + edges[:-1])
    return kc, Ek


def main():
    files = sorted(glob.glob(os.path.join(DUMP_DIR, "vel_*.bin")))
    if not files:
        print("Нет дампов в", os.path.abspath(DUMP_DIR),
              "\nВключи _dumpVelocity на SFHybridCS и прогони сцену (желательно при разных hbar).")
        return
    # последний дамп для каждого hbar
    by_hbar = {}
    for p in files:
        d = load(p)
        by_hbar[round(d["hbar"], 4)] = p  # перезапись → останется последний (наибольший step)

    plt.figure(figsize=(7, 5))
    for hbar in sorted(by_hbar):
        d = load(by_hbar[hbar])
        kc, Ek = spectrum(d)
        mask = (kc > 0) & (Ek > 0)
        plt.loglog(kc[mask], Ek[mask], lw=1.8, label=f"ℏ = {hbar:g}")

    # опорный наклон Колмогорова k^(-5/3)
    kref = np.array([kc[mask][1], kc[mask][len(kc[mask]) // 2]])
    if len(kref) == 2 and kref[0] > 0:
        c = Ek[mask][1] * kref[0] ** (5 / 3)
        plt.loglog(kref, c * kref ** (-5 / 3), "k--", lw=1, label="k$^{-5/3}$ (Колмогоров)")

    plt.xlabel("волновое число k")
    plt.ylabel("E(k)")
    plt.title("Энергетический спектр поля скорости при разных ℏ")
    plt.legend()
    plt.grid(True, which="both", ls=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(OUT, dpi=150)
    print("сохранено:", OUT)


if __name__ == "__main__":
    main()
