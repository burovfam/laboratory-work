import matplotlib.pyplot as plt
import numpy as np
import os

r = np.array([3, 4, 6, 8, 9, 10])
T = np.array([2.05, 2.25, 2.80, 3.40, 3.65, 3.90])

x = T**2
y = r**2

coefficients = np.polyfit(x, y, 1)
k = coefficients[0]
b = coefficients[1]

x_trend = np.linspace(0, 18, 100)
y_trend = k * x_trend + b

plt.figure(figsize=(8, 6))

plt.scatter(x, y, color='red', s=50, zorder=3, label='Экспериментальные точки')
plt.plot(x_trend, y_trend, color='blue', linewidth=2, label='Линия тренда (МНК)')

plt.title('График зависимости $r^2$ от $T^2$', fontsize=14)
plt.xlabel('$T^2$, с$^2$', fontsize=12)
plt.ylabel('$r^2$, см$^2$', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

equation_text = f'$y = {k:.2f}x + {b:.2f}$'
plt.text(0.05, 0.95, equation_text, transform=plt.gca().transAxes, fontsize=12,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

script_dir = os.path.dirname(os.path.abspath(__file__))
save_path = os.path.join(script_dir, 'graph.png')

plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"График сохранен как: {save_path}")
print(f"Уравнение тренда: y = {k:.4f} * x + {b:.4f}")
print(f"Коэффициент наклона k = {k:.4f}")