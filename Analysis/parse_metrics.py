#!/usr/bin/env python3
"""
Парсер журнала [Hybrid3D] из Unity Editor.log → CSV + графики для раздела «проверяемые результаты».
Метрики: alphaSum (сохранение/баланс α), occupied% (охват), |u|mean (устойчивость потока),
Ekin (энергия), repErr/divRMS (ошибка представимости/несжимаемость), drift% (тест сохранения объёма).

Запуск:  python3 Analysis/parse_metrics.py [путь_к_логу]
По умолчанию берёт ~/Library/Logs/Unity/Editor.log
"""
import os, re, sys, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Library/Logs/Unity/Editor.log")
CSV_OUT = os.path.join(HERE, "metrics.csv")
PNG_OUT = os.path.join(HERE, "metrics.png")

FIELDS = ["step", "hbar", "alphaSum", "alphaMax", "occupied_pct", "u_mean", "u_max",
          "Ekin", "repErr", "divRMS", "drift_pct"]

PATTERNS = {
    "step": r"step=(\d+)",
    "hbar": r"hbar=([\d.]+)",
    "alphaSum": r"alphaSum=([\d.eE+-]+)",
    "alphaMax": r"alphaMax=([\d.eE+-]+)",
    "occupied_pct": r"\(([\d.]+)%\)",
    "u_mean": r"\|u\|mean=([\d.eE+-]+)",
    "u_max": r"\|u\|max=([\d.eE+-]+)",
    "Ekin": r"Ekin=([\d.eE+-]+)",
    "repErr": r"repErr=([\d.eE+-]+)",
    "divRMS": r"divRMS=([\d.eE+-]+)",
    "drift_pct": r"drift=([+\-]?[\d.]+)%",
}


def parse():
    rows = []
    with open(LOG, "r", errors="ignore") as f:
        for line in f:
            if "[Hybrid3D]" not in line or "step=" not in line:
                continue
            row = {}
            for k, pat in PATTERNS.items():
                m = re.search(pat, line)
                row[k] = float(m.group(1)) if m else ""
            rows.append(row)
    return rows


def main():
    if not os.path.exists(LOG):
        print("лог не найден:", LOG); return
    rows = parse()
    if not rows:
        print("Нет строк [Hybrid3D] в логе. Запусти сцену в Play (метрики печатаются каждые _debugEverySteps).")
        return
    with open(CSV_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"сохранено: {CSV_OUT} ({len(rows)} строк)")

    step = [r["step"] for r in rows]
    def col(name):
        return [r[name] if r[name] != "" else float("nan") for r in rows]

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].plot(step, col("u_mean")); ax[0, 0].set_title("|u|mean — устойчивость потока"); ax[0, 0].set_xlabel("шаг")
    ax[0, 1].plot(step, col("Ekin"), color="tab:orange"); ax[0, 1].set_title("Ekin — кинетическая энергия"); ax[0, 1].set_xlabel("шаг")
    ax[1, 0].semilogy(step, col("repErr"), color="tab:red"); ax[1, 0].set_title("repErr — ошибка представимости (норм. дивергенция)"); ax[1, 0].set_xlabel("шаг")
    drift = col("drift_pct")
    if any(d == d for d in drift):  # есть not-NaN
        ax[1, 1].plot(step, drift, color="tab:green"); ax[1, 1].set_title("drift% — сохранение объёма ∫α")
    else:
        ax[1, 1].plot(step, col("occupied_pct"), color="tab:green"); ax[1, 1].set_title("occupied% — охват объёма")
    ax[1, 1].set_xlabel("шаг")
    for a in ax.ravel():
        a.grid(True, ls=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(PNG_OUT, dpi=150)
    print("сохранено:", PNG_OUT)


if __name__ == "__main__":
    main()
