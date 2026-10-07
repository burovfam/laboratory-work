import numpy as np
import matplotlib.pyplot as plt

# данные момента сил и угла поворота для 3 серий
M1 = np.array([0.0513, 0.1025, 0.1538, 0.2050, 0.2563, 0.3075, 0.4100, 0.5125,
               0.4100, 0.3075, 0.2563, 0.2050, 0.1538, 0.1025, 0.0513])
phi1 = np.array([0.0004, 0.0085, 0.0165, 0.0266, 0.0346, 0.0441, 0.0606, 0.0788,
                 0.0610, 0.0441, 0.0330, 0.0266, 0.0172, 0.0085, 0.0007])

M2 = np.array([0.0513, 0.1025, 0.1538, 0.2050, 0.2563, 0.3075, 0.4100, 0.5125,
               0.5125, 0.4100, 0.3075, 0.2563, 0.2050, 0.1538, 0.1025, 0.0513])
phi2 = np.array([-0.0003, 0.0077, 0.0175, 0.0251, 0.0340, 0.0434, 0.0593, 0.0801,
                 0.0801, 0.0606, 0.0431, 0.0347, 0.0249, 0.0175, 0.0088, -0.0002])

M3 = np.array([0.0513, 0.1025, 0.1538, 0.2050, 0.2563, 0.3075, 0.4100, 0.5125,
               0.4100, 0.3075, 0.2563, 0.2050, 0.1538, 0.1025, 0.0513])
phi3 = np.array([0.0003, 0.0077, 0.0172, 0.0259, 0.0320, 0.0438, 0.0623, 0.0791,
                 0.0606, 0.0431, 0.0350, 0.0256, 0.0172, 0.0077, -0.0003])

# работа функции МНК, определение углового коэффициента, который потом перейдет в модуль кручения
f1 = (np.mean(M1*phi1) - np.mean(M1)*np.mean(phi1)) / (np.mean(phi1**2) - np.mean(phi1)**2)
f2 = (np.mean(M2*phi2) - np.mean(M2)*np.mean(phi2)) / (np.mean(phi2**2) - np.mean(phi2)**2)
f3 = (np.mean(M3*phi3) - np.mean(M3)*np.mean(phi3)) / (np.mean(phi3**2) - np.mean(phi3)**2)

b1 = np.mean(M1) - f1*np.mean(phi1)
b2 = np.mean(M2) - f2*np.mean(phi2)
b3 = np.mean(M3) - f3*np.mean(phi3)

# теор данные для модуля сдвига
l = 1.485
R = 0.003
# поиск модуля сдвига
G1 = 2*l*f1 / (np.pi*R**4)
G2 = 2*l*f2 / (np.pi*R**4)
G3 = 2*l*f3 / (np.pi*R**4)

# погрешности измерений (в СИ)
dM = 0.001
dphi = 0.0001
dR = 0.0000025
dl = 0.0005

df1_rel = np.sqrt((dM/np.mean(M1))**2 + (dphi/np.mean(phi1))**2)
df2_rel = np.sqrt((dM/np.mean(M2))**2 + (dphi/np.mean(phi2))**2)
df3_rel = np.sqrt((dM/np.mean(M3))**2 + (dphi/np.mean(phi3))**2)

df1 = f1*df1_rel
df2 = f2*df2_rel
df3 = f3*df3_rel

dG1_rel = np.sqrt(df1_rel**2 + (4*dR/R)**2 + (dl/l)**2)
dG2_rel = np.sqrt(df2_rel**2 + (4*dR/R)**2 + (dl/l)**2)
dG3_rel = np.sqrt(df3_rel**2 + (4*dR/R)**2 + (dl/l)**2)

dG1 = G1*dG1_rel
dG2 = G2*dG2_rel
dG3 = G3*dG3_rel

# средние значения f и G, средние значения погрешностей
f_sr = (f1+f2+f3)/3
G_sr = (G1+G2+G3)/3
df_sr = (df1+df2+df3)/3
dG_sr = (dG1+dG2+dG3)/3

# построение графиков (оси перевёрнуты: M по X, φ по Y)
fig, ax = plt.subplots(1, 3, figsize=(15,5))

for i, (M, phi, f, b, title) in enumerate(zip(
    [M1, M2, M3], [phi1, phi2, phi3], [f1, f2, f3], [b1, b2, b3],
    ["Серия 1", "Серия 2", "Серия 3"]
)):
    # точки: по X — момент, по Y — угол
    ax[i].scatter(M, phi, color='red', s=50, label='Точки', zorder=5)

    # МНК-прямая в новых координатах: φ = (1/f)·M - b/f
    x = np.linspace(min(M), max(M), 100)
    ax[i].plot(x, (x - b)/f, 'b-', linewidth=2,
               label=f'МНК: φ = {1/f:.3f}·M {-b/f:+.4f}')

    ax[i].set_xlabel('M, Н·м')
    ax[i].set_ylabel('φ, рад')
    ax[i].set_title(title)
    ax[i].legend()
    ax[i].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# вывод необходимых значений для самой таблицы
print("Серия 1: f =", round(f1,2), "±", round(df1,2), "Н·м/рад, G =", round(G1/1e9,1), "±", round(dG1/1e9,1), "ГПа")
print("Серия 2: f =", round(f2,2), "±", round(df2,2), "Н·м/рад, G =", round(G2/1e9,1), "±", round(dG2/1e9,1), "ГПа")
print("Серия 3: f =", round(f3,2), "±", round(df3,2), "Н·м/рад, G =", round(G3/1e9,1), "±", round(dG3/1e9,1), "ГПа")
print("Среднее: f =", round(f_sr,2), "±", round(df_sr,2), "Н·м/рад")
print("Среднее: G =", round(G_sr/1e9,1), "±", round(dG_sr/1e9,1), "ГПа")