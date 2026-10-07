import numpy as np
import matplotlib.pyplot as plt
from math import factorial


file = r"C:\Users\basil\Visual Studio Code Projects\laboratory-work\1 course\1.1.4\work_file.txt"

f = [x for x in open(file).readlines() if not x.startswith("#")]


def sum_values(f, p):

    a = []

    for i in range(0, len(f), p):

        k = sm = 0

        while k != p and i + k < len(f):

            sm += int(f[i+k])
            k += 1

        if k == p:
            a.append(sm)

    n = sorted(set(a))
    wn = []
    N = len(a)

    for x in n:
        wn.append(a.count(x) / N)

    plt.bar(n, wn, width=1, alpha=0.5,
            label=f"τ = {p} с",
            linewidth=2, edgecolor="black")

    return a, n, wn


def stats(a, p):

    n_mean = np.mean(a)
    sigma_n = np.std(a, ddof=0)
    N = len(a)

    sigma_mean = sigma_n / np.sqrt(N)

    j = n_mean / p
    sigma_j = sigma_mean / p

    sigma_theory = np.sqrt(n_mean)

    print(f"\nРасчёт статистик для τ = {p} с")
    print(f"N = {N}")
    print(f"<n> = {n_mean}")
    print(f"σ_n = {sigma_n}")
    print(f"σ_<n> = {sigma_mean}")
    print(f"j = {j}")
    print(f"σ_j = {sigma_j}")

    print(f"σ_n (эксперимент) = {sigma_n}")
    print(f"√<n> (теория Пуассона) = {sigma_theory}")
    print(f"Разность = {abs(sigma_n - sigma_theory)}")

    return n_mean, sigma_n, sigma_mean, j, sigma_j


def poisson(n, n_mean):

    return n_mean**n / factorial(n) * np.exp(-n_mean)


def gaussian(n, n_mean, sigma_n):

    return 1 / (sigma_n * np.sqrt(2 * np.pi)) * \
           np.exp(-(n - n_mean)**2 / (2 * sigma_n**2))


def fractions(a, n_mean, sigma_n):

    N = len(a)

    f1 = sum(abs(x - n_mean) <= sigma_n for x in a) / N
    f2 = sum(abs(x - n_mean) <= 2 * sigma_n for x in a) / N
    f3 = sum(abs(x - n_mean) <= 3 * sigma_n for x in a) / N

    print(f"\nДоли попаданий для <n> = {n_mean}")
    print(f"±σ:  {f1}  ({f1 * 100:.1f}%)")
    print(f"±2σ: {f2}  ({f2 * 100:.1f}%)")
    print(f"±3σ: {f3}  ({f3 * 100:.1f}%)")

    return f1, f2, f3



a10, n10, wn10 = sum_values(f, 10)
a20, n20, wn20 = sum_values(f, 20)
a40, n40, wn40 = sum_values(f, 40)
a80, n80, wn80 = sum_values(f, 80)




plt.xlabel("n")
plt.ylabel("$w_n$")
plt.title("Гистограммы распределения")
plt.legend()
plt.grid(axis="y")

plt.show()




s10 = stats(a10, 10)
s20 = stats(a20, 20)
s40 = stats(a40, 40)
s80 = stats(a80, 80)



fractions(a10, s10[0], s10[1])
fractions(a20, s20[0], s20[1])
fractions(a40, s40[0], s40[1])
fractions(a80, s80[0], s80[1])


# Эксперимент + распределение Пуассона

plt.figure()

plt.bar(n20, wn20, width=1, alpha=0.5,
        label="Эксперимент",
        linewidth=2, edgecolor="black")

x = np.arange(min(n20), max(n20) + 1)

w_poisson = [poisson(n, s20[0]) for n in x]

plt.plot(x, w_poisson, "o-",
         label="Распределение Пуассона")

plt.xlabel("n")
plt.ylabel("$w_n$")
plt.title("Экспериментальное распределение и распределение Пуассона")
plt.legend()
plt.grid(axis="y")

plt.show()


# Эксперимент + распределение Гаусса


plt.figure()

plt.bar(n20, wn20, width=1, alpha=0.5,
        label="Эксперимент",
        linewidth=2, edgecolor="black")

x = np.linspace(min(n20), max(n20), 500)

w_gaussian = gaussian(x, s20[0], s20[1])

plt.plot(x, w_gaussian,
         label="Нормальное распределение")

plt.xlabel("n")
plt.ylabel("$w_n$")
plt.title("Экспериментальное распределение и распределение Гаусса")
plt.legend()
plt.grid(axis="y")

plt.show()


# Пуассон + Гаусс


plt.figure()

plt.bar(n20, wn20, width=1, alpha=0.5,
        label="Эксперимент",
        linewidth=2, edgecolor="black")

x = np.arange(min(n20), max(n20) + 1)

w_poisson = [poisson(n, s20[0]) for n in x]

plt.plot(x, w_poisson, "o-",
         label="Пуассон")

x_gauss = np.linspace(min(n20), max(n20), 500)

w_gaussian = gaussian(x_gauss, s20[0], s20[1])

plt.plot(x_gauss, w_gaussian,
         label="Гаусс")

plt.xlabel("n")
plt.ylabel("$w_n$")
plt.title("Экспериментальное распределение, Пуассон и Гаусс")
plt.legend()
plt.grid(axis="y")

plt.show()