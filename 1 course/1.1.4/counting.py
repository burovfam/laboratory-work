import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson, norm
import sys

base = r"C:\Users\basil\Visual Studio Code Projects\laboratory-work\1 course\1.1.4"

output_file = open(base + r"\results.txt", "w", encoding="utf-8")
sys.stdout = output_file

# ============================================================
# п.14: Группировка данных по интервалам τ
# ============================================================

f = open(base + r"\search.txt").readlines()

def res(f, p):
    a = []
    for i in range(0, len(f), p):
        k = sm = 0
        while k != p and i + k < len(f):
            sm += int(f[i + k])
            k += 1
        if k == p:
            a.append(sm)

    with open(fr"{base}\result_{p}.txt", "w") as file:
        for x in a:
            file.write(str(x) + "\n")

res(f, 10)
res(f, 20)
res(f, 40)
res(f, 80)

def load(p):
    with open(fr"{base}\result_{p}.txt", "r") as file:
        return [int(line.strip()) for line in file.readlines()]

a10 = load(10)
a20 = load(20)
a40 = load(40)
a80 = load(80)


# ============================================================
# п.16: Расчёт статистик <n>, σ_n, σ_<n>, j, σ_j для каждого τ
# ============================================================

def stats(data, tau):
    n_mean = np.mean(data)
    sigma_n = np.std(data, ddof=0)
    N = len(data)
    sigma_mean = sigma_n / np.sqrt(N)
    j = n_mean / tau
    sigma_j = sigma_mean / tau

    print(f"--- τ = {tau} с ---")
    print(f"N = {N}")
    print(f"<n> = {n_mean:.2f}")
    print(f"σ_n = {sigma_n:.2f}   (теория: √<n> = {np.sqrt(n_mean):.2f})")
    print(f"σ_<n> = {sigma_mean:.3f}")
    print(f"j = {j:.3f} ± {sigma_j:.3f}  частиц/с")
    print()

    return n_mean, sigma_n, sigma_mean, j, sigma_j

n_mean_10, sigma_n_10, sigma_mean_10, j_10, sigma_j_10 = stats(a10, 10)
n_mean_20, sigma_n_20, sigma_mean_20, j_20, sigma_j_20 = stats(a20, 20)
n_mean_40, sigma_n_40, sigma_mean_40, j_40, sigma_j_40 = stats(a40, 40)
n_mean_80, sigma_n_80, sigma_mean_80, j_80, sigma_j_80 = stats(a80, 80)


# ============================================================
# п.15: Гистограммы для каждого τ (4 subplot'а)
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

def draw_hist(ax, data, tau):
    bins = range(min(data), max(data) + 2)
    ax.hist(data, bins=bins, density=True, align='left',
            edgecolor='black', alpha=0.7)
    ax.set_title(f"tau = {tau} с")
    ax.set_xlabel("n (число отсчётов)")
    ax.set_ylabel("$w_n$")

draw_hist(axes[0, 0], a10, 10)
draw_hist(axes[0, 1], a20, 20)
draw_hist(axes[1, 0], a40, 40)
draw_hist(axes[1, 1], a80, 80)

plt.tight_layout()
plt.savefig(base + r"\histograms.png", dpi=150)
plt.close(fig)


# ============================================================
# п.18: Проверка свойства Пуассона σ_n ≈ √<n>
# ============================================================

def check_poisson_property(tau, n_mean, sigma_n):
    theory_sigma = np.sqrt(n_mean)
    relative_diff = abs(sigma_n - theory_sigma) / theory_sigma * 100

    print(f"--- τ = {tau} с: проверка свойства Пуассона ---")
    print(f"σ_n (эксперимент) = {sigma_n:.2f}")
    print(f"√<n> (теория)     = {theory_sigma:.2f}")
    print(f"расхождение       = {relative_diff:.1f}%")
    print()

check_poisson_property(10, n_mean_10, sigma_n_10)
check_poisson_property(20, n_mean_20, sigma_n_20)
check_poisson_property(40, n_mean_40, sigma_n_40)
check_poisson_property(80, n_mean_80, sigma_n_80)


# ============================================================
# п.19: Доли случаев в пределах σ, 2σ, 3σ
# ============================================================

def check_sigma_intervals(data, tau, n_mean, sigma_n):
    N = len(data)
    data_np = np.array(data)
    deviations = np.abs(data_np - n_mean)

    count_1s = np.sum(deviations <= sigma_n)
    count_2s = np.sum(deviations <= 2 * sigma_n)
    count_3s = np.sum(deviations <= 3 * sigma_n)

    frac_1s = count_1s / N * 100
    frac_2s = count_2s / N * 100
    frac_3s = count_3s / N * 100

    print(f"--- τ = {tau} с: доли в пределах σ, 2σ, 3σ ---")
    print(f"|n-<n>| <= σ_n  : {count_1s}/{N} = {frac_1s:.1f}%   (теория Гаусса: 68.3%)")
    print(f"|n-<n>| <= 2σ_n : {count_2s}/{N} = {frac_2s:.1f}%   (теория Гаусса: 95.4%)")
    print(f"|n-<n>| <= 3σ_n : {count_3s}/{N} = {frac_3s:.1f}%   (теория Гаусса: 99.7%)")
    print()

check_sigma_intervals(a10, 10, n_mean_10, sigma_n_10)
check_sigma_intervals(a20, 20, n_mean_20, sigma_n_20)
check_sigma_intervals(a40, 40, n_mean_40, sigma_n_40)
check_sigma_intervals(a80, 80, n_mean_80, sigma_n_80)


# ============================================================
# п.17: Наложение теоретических распределений (Пуассон, Гаусс)
# ============================================================

def plot_with_theory(data, tau, n_mean, sigma_n):
    fig, ax = plt.subplots(figsize=(8, 6))

    bins = np.arange(min(data), max(data) + 2)
    ax.hist(data, bins=bins, density=True, align='left',
            edgecolor='black', alpha=0.6, label="эксперимент")

    n_range = np.arange(min(data), max(data) + 1)

    poisson_pmf = poisson.pmf(n_range, mu=n_mean)
    ax.plot(n_range, poisson_pmf, 'ro-', label="теория (Пуассон)", markersize=4)

    x_smooth = np.linspace(min(data), max(data), 300)
    gauss_pdf = norm.pdf(x_smooth, loc=n_mean, scale=sigma_n)
    ax.plot(x_smooth, gauss_pdf, 'g--', label="теория (Гаусс)")

    ax.set_title(f"tau = {tau} с, наложение теоретических распределений")
    ax.set_xlabel("n")
    ax.set_ylabel("$w_n$")
    ax.legend()

    plt.tight_layout()
    plt.savefig(base + fr"\histogram_theory_tau{tau}.png", dpi=150)
    plt.close(fig)

plot_with_theory(a10, 10, n_mean_10, sigma_n_10)


sys.stdout = sys.__stdout__
output_file.close()

print("finish")