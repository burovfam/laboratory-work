import math
import matplotlib.pyplot as plt

g = 9.81
R = 0.1146
r = 0.0302
z0 = 2.20
sigma_R = 0.0005
sigma_r = 0.0003
sigma_z0 = 0.005

k = g * R * r / (4 * math.pi ** 2 * z0)
rel_k = math.sqrt((sigma_R / R) ** 2 + (sigma_r / r) ** 2 + (sigma_z0 / z0) ** 2)
sigma_k = k * rel_k
print("k =", k, " sigma_k/k =", rel_k, " sigma_k =", sigma_k)

t0 = [65.26, 65.05, 70.50, 86.79, 43.82]
n0 = [15, 15, 16, 20, 10]
T0_list = []
for i in range(5):
    T0_list.append(t0[i] / n0[i])
T0 = sum(T0_list) / 5
summa = 0
for x in T0_list:
    summa = summa + (x - T0) ** 2
sigma_rand = math.sqrt(summa / (5 * 4))
sigma_syst = 0.01 / 15
sigma_T0 = math.sqrt(sigma_rand ** 2 + sigma_syst ** 2)
print("T0 =", T0, " sigma_rand =", sigma_rand, " sigma_syst =", sigma_syst, " sigma_T0 =", sigma_T0)

m0 = 1.0668
sigma_m = 0.0005
I0 = k * m0 * T0 ** 2
rel_I0 = math.sqrt(rel_k ** 2 + (sigma_m / m0) ** 2 + (2 * sigma_T0 / T0) ** 2)
sigma_I0 = I0 * rel_I0
print("I0 =", I0, " rel =", rel_I0, " sigma =", sigma_I0)


def telo(imya, m, T, sigma_T):
    M = m0 + m
    I_sys = k * M * T ** 2
    rel = math.sqrt(rel_k ** 2 + (2 * sigma_T / T) ** 2)
    sigma_sys = I_sys * rel
    I_t = I_sys - I0
    sigma_t = math.sqrt(sigma_sys ** 2 + sigma_I0 ** 2)
    print(imya, ": M =", M, " I_sys =", I_sys, " sigma_sys =", sigma_sys, " I =", I_t, " sigma =", sigma_t, " T rel =", sigma_T / T)
    return I_t, sigma_t


I1, s1 = telo("ободок", 0.7476, 4.195, 0.035)
I2, s2 = telo("крышка", 0.590, 3.936, 0.002)

M12 = m0 + 0.7476 + 0.590
I12 = k * M12 * 3.914 ** 2 - I0
I12_teor = I1 + I2
s12 = math.sqrt(s1 ** 2 + s2 ** 2)
print("M12 =", M12, " I12_eksp =", I12, " I12_teor =", I12_teor, " delta =", abs(I12 - I12_teor), " sigma_sum =", s12)

md = 1.4418
M = m0 + md
rows = [(9, 62.39), (10, 62.60), (11, 62.81), (12, 63.21), (13, 64.14),
        (14, 65.20), (15, 66.49), (16, 67.56), (17, 68.96), (18, 71.34)]
hs = []
h2 = []
T_d = []
I_d = []
for nn, tt in rows:
    h = (nn - 9) * 0.005
    T = tt / 20
    I = k * M * T ** 2 - I0
    hs.append(h)
    h2.append(h * h)
    T_d.append(T)
    I_d.append(I)
    print(nn, round(h, 3), round(h * h, 6), round(tt, 2), round(T, 4), round(I, 5))

n = len(h2)
mx = sum(h2) / n
my = sum(I_d) / n
sxx = 0
sxy = 0
for i in range(n):
    sxx = sxx + (h2[i] - mx) ** 2
    sxy = sxy + (h2[i] - mx) * (I_d[i] - my)
a = sxy / sxx
b = my - a * mx
ostatki = 0
syy = 0
for i in range(n):
    ostatki = ostatki + (I_d[i] - a * h2[i] - b) ** 2
    syy = syy + (I_d[i] - my) ** 2
s2_ost = ostatki / (n - 2)
sigma_a = math.sqrt(s2_ost / sxx)
sigma_b = math.sqrt(s2_ost * (1 / n + mx ** 2 / sxx))
r2 = 1 - ostatki / syy
print("МНК: a =", a, "+-", sigma_a, " b =", b, "+-", sigma_b, " R^2 =", r2)
print("a / m_d =", a / md, " (по теории a = m_d =", md, ")")

R_diska = 9 * 0.005
I_teor0 = md * R_diska ** 2 / 2
print("теория: I_д0 = m R^2/2 =", I_teor0, " (R =", R_diska, "м)")
print("a - m_d =", a - md, " в сигмах:", (a - md) / sigma_a)

x_line = [0, 0.0022]
y_line = []
y_teor = []
for x in x_line:
    y_line.append((a * x + b) * 1000)
    y_teor.append((md * x + I_teor0) * 1000)
y_points = []
for v in I_d:
    y_points.append(v * 1000)

plt.rcParams["font.size"] = 11
plt.figure(figsize=(7.5, 5))
plt.plot(h2, y_points, "o", color="tab:blue", markersize=6, markeredgecolor="black", markeredgewidth=0.6, label="Эксперимент")
plt.plot(x_line, y_line, "-", color="tab:red", linewidth=1.8, label="Линейная аппроксимация (МНК)")
plt.plot(x_line, y_teor, "--", color="gray", linewidth=1.5, label="Теория: $I = MR^2/2 + Mh^2$")
tekst = "$I = ah^2 + b$\n$a = %.2f \\pm %.2f$ кг\n$b = (%.2f \\pm %.2f)\\cdot10^{-3}$ кг·м$^2$\n$R^2 = %.3f$" % (a, sigma_a, b * 1000, sigma_b * 1000, r2)
plt.text(0.97, 0.05, tekst, transform=plt.gca().transAxes, ha="right", va="bottom", fontsize=10,
         bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.95))
plt.xlim(0, 0.0022)
plt.ylim(bottom=0)
plt.xlabel("$h^2$, м$^2$")
plt.ylabel("Момент инерции диска $I$, $10^{-3}$ кг·м$^2$")
plt.title("Разрезанный диск: зависимость $I$ от $h^2$")
plt.grid(alpha=0.35)
plt.legend(loc="upper left", fontsize=10)
plt.tight_layout()
plt.savefig("graph_disk.png", dpi=200)
print("graph_disk.png сохранён")
