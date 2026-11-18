import numpy as np
import matplotlib.pyplot as plt

M = np.array([0.0660409, 0.1078315, 0.1668092,
              0.2575714, 0.4007189])

Omega = np.array([0.0356791, 0.0580776, 0.0903639,
                  0.1392735, 0.2190943])

# МНК: Omega = a*M + b
a, b = np.polyfit(M, Omega, 1)

# Теоретическая зависимость
I0 = 8.3475e-4
nu = 395.3
omega0 = 2 * np.pi * nu
a_theory = 1 / (I0 * omega0)

M_line = np.linspace(0, 0.42, 300)

plt.figure(figsize=(8, 5.5))
plt.scatter(M, Omega, label="Эксперимент")
plt.plot(M_line, a * M_line + b,
         label=f"МНК: Ω = {a:.4f}M {b:+.4f}")
plt.plot(M_line, a_theory * M_line, "--",
         label="Теория")

plt.xlabel(r"$M$, Н$\cdot$м")
plt.ylabel(r"$\Omega$, рад/с")
plt.title("Зависимость скорости прецессии от момента силы")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("grafik_precessiya_mnk.png", dpi=300)
plt.show()