import numpy as np
import matplotlib.pyplot as plt

# ---------- исходные данные ----------
g = 9.81
dm = np.array([245.2, 245.5, 245.8, 245.3, 245.5, 245.8, 245.2, 245.3, 245.9, 245.3])  # г
P = np.cumsum(dm) / 1000 * g                                                           # Н

n_1_down = np.array([27.4, 25.6, 23.9, 22.5, 20.8, 19.7, 17.8, 16.2, 14.6, 13.2])      # см
n_1_up   = np.array([31.0, 28.8, 26.8, 25.0, 23.2, 21.5, 19.7, 18.0, 16.8, 14.7])
n_2_down = np.array([28.2, 26.4, 24.7, 22.9, 21.7, 19.6, 17.8, 16.7, 14.5, 12.7])
n_2_up_raw = np.array([14.4, 16.2, 17.7, 19.4, 21.1, 22.9, 24.7, 26.5, 28.4, 30.4])


n_2_up = n_2_up_raw[::-1][1:]           

d, sigma_d = 0.51, 0.01                  
r, sigma_r = 20.0, 1.0                  
l, sigma_l = 1330.0, 5.0                 
h, sigma_h = 1760.0, 5.0                 
sigma_n = 0.1                           

c = r / (2 * h) * 10                     # мм удлинения на 1 см шкалы: Δl = n·r/(2h)
sigma_dl = np.sqrt(2) * sigma_n * c      

series = {
    "серия 1, нагружение": (P,     n_1_down, "o", "tab:blue"),
    "серия 1, разгрузка":  (P,     n_1_up,   "s", "tab:orange"),
    "серия 2, нагружение": (P,     n_2_down, "^", "tab:green"),
    "серия 2, разгрузка":  (P[:9], n_2_up,   "v", "tab:red"),
}


data = {name: (x, c * (n[0] - n)) for name, (x, n, _, _) in series.items()}
Pall = np.concatenate([x for x, _ in data.values()])
Lall = np.concatenate([dl for _, dl in data.values()])
N = len(Pall)

# ---------- МНК по всем точкам: P = k·Δl + P0 ----------
mean = np.mean
k = (mean(Pall * Lall) - mean(Pall) * mean(Lall)) / (mean(Lall**2) - mean(Lall)**2)   # Н/мм
P0 = mean(Pall) - k * mean(Lall)
eps_rand = 1 / (k * np.sqrt(N - 2)) * np.sqrt(
    (mean(Pall**2) - mean(Pall)**2) / (mean(Lall**2) - mean(Lall)**2) - k**2)
eps_dl = np.hypot(sigma_r / r, sigma_h / h)       # систематическая, eps_P пренебрежимо мала
eps_k = np.hypot(eps_rand, eps_dl)

S = np.pi * d**2 / 4
E = k * l / S                                      # Н/мм² = МПа
eps_S, eps_l = 2 * sigma_d / d, sigma_l / l
eps_E = np.sqrt(eps_S**2 + eps_k**2 + eps_l**2)

print(f"N = {N}")
print(f"k = {k:.2f} Н/мм")
print(f"eps_случ = {eps_rand*100:.1f}%  eps_сист = {eps_dl*100:.1f}%  eps_k = {eps_k*100:.1f}%")
print(f"S = {S:.3f} мм²")
print(f"E = {E/1000:.0f} ГПа  eps_E = {eps_E*100:.1f}%  sigma_E = {E/1000*eps_E:.0f} ГПа")

# ---------- график ----------
fig, ax = plt.subplots(figsize=(7.5, 5))
for name, (x, dl) in data.items():
    _, _, marker, color = series[name]
    yerr = np.where(np.arange(len(x)) == 0, 0, sigma_dl)  
    ax.errorbar(x, dl, yerr=yerr, fmt=marker, color=color, capsize=3,
                elinewidth=1, markersize=4, label=name)

xx = np.array([P.min(), P.max()])
ax.plot(xx, (xx - P0) / k, "k-", lw=1.2, label=f"МНК по всем точкам ($k={k:.1f}$ Н/мм)")

ax.set_xlabel("$P$, Н")
ax.set_ylabel(r"$\Delta l$, мм")
ax.set_xlim(0, 25)
ax.set_ylim(-0.05, None)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("graph.png", dpi=170)
plt.show()