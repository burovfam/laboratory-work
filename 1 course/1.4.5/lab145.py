import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L, dL = 0.5, 0.001         
rho_tab = 568.4e-6         
g = 9.81
dm = 0.01e-3                

m_g = [110.6+461.4+492.7]
for add in (480.7, 500.3, 499.1, 481.7):
    m_g.append(m_g[-1] + add)
m = np.array(m_g) * 1e-3
T = m * g
dT = T * dm / m             

nu = [
 [135.3, 271.6, 402.1, 542.2, 677.0, 812.9, 948.9],
 [161.4, 323.2, 485.9, 648.4, 808.3, 969.6, 1130.9],
 [186.7, 372.7, 560.5, 746.6, 931.8],
 [208.4, 416.7, 624.2, 833.1, 1041.3],
 [230.8, 455.0, 682.4, 909.1, 1137.1],
]

# --- МНК-прямая nu = k*n через начало координат ---
k, sk, u, su = [], [], [], []
for y in nu:
    y = np.array(y); n = np.arange(1, len(y)+1)
    kk = np.sum(n*y)/np.sum(n**2)
    res = y - kk*n
    s_k = np.sqrt(np.sum(res**2)/((len(y)-1)*np.sum(n**2)))
    uu = 2*L*kk
    s_u = uu*np.sqrt((dL/L)**2 + (s_k/kk)**2)
    k.append(kk); sk.append(s_k); u.append(uu); su.append(s_u)
    print(f"k={kk:.3f}±{s_k:.3f} Гц  u={uu:.2f}±{s_u:.2f} м/с  sigma_nu(остатки)={np.sqrt(np.sum(res**2)/(len(y)-1)):.2f} Гц")
k, sk, u, su = map(np.array, (k, sk, u, su))

# --- u^2 = T/rho: прямая через начало координат ---
u2 = u**2; su2 = 2*u*su
w = 1/su2**2
a = np.sum(w*T*u2)/np.sum(w*T**2)       # a = 1/rho
sa = 1/np.sqrt(np.sum(w*T**2))
chi2 = np.sum(w*(u2-a*T)**2)/(len(T)-1)
print("chi2/ndf =", chi2)
rho = 1/a
srho_fit = sa/a**2
srho = rho*np.sqrt((srho_fit/rho)**2 + (2*dL/L)**2)
print(f"1/rho = {a:.1f}±{sa:.1f} ;  rho = {rho*1e3:.4f}±{srho*1e3:.4f} г/м  (fit {srho_fit*1e3:.4f})")
print(f"табл. {rho_tab*1e3:.4f} г/м, расхождение {(rho-rho_tab)/rho_tab*100:+.2f}%")

# теория
print("\nТеория по точным массам:")
for i in range(5):
    ut = np.sqrt(T[i]/rho_tab)
    print(f"m{i+1}: m={m[i]*1e3:.1f} г, T={T[i]:.3f} Н, u_теор={ut:.2f}, u_эксп={u[i]:.2f}±{su[i]:.2f}, откл {(u[i]-ut)/ut*100:+.2f}%")

# --- графики ---
fig, ax = plt.subplots(figsize=(8,6))
for i, y in enumerate(nu):
    n = np.arange(1, len(y)+1)
    p = ax.plot(n, y, "o", label=f"T = {T[i]:.2f} Н")
    nn = np.array([0, 7.3])
    ax.plot(nn, k[i]*nn, "-", color=p[0].get_color(), lw=1)
ax.set_xlabel("Номер гармоники $n$"); ax.set_ylabel(r"$\nu_n$, Гц")
ax.set_xlim(0, 7.5); ax.set_ylim(0, None); ax.grid(alpha=.3); ax.legend()
plt.tight_layout(); plt.savefig("nu_vs_n.png", dpi=200); plt.close()

fig, ax = plt.subplots(figsize=(8,6))
ax.errorbar(T, u2, yerr=su2, fmt="o", capsize=3, label="эксперимент")
tt = np.linspace(0, T.max()*1.05, 50)
ax.plot(tt, a*tt, "-", label=fr"МНК: $u^2 = T/\rho_l$, $\rho_l$ = {rho*1e3:.3f} г/м")
ax.plot(tt, tt/rho_tab, "--", label=fr"табличное $\rho_l$ = {rho_tab*1e3:.3f} г/м")
ax.set_xlabel("$T$, Н"); ax.set_ylabel("$u^2$, м$^2$/с$^2$")
ax.set_xlim(0, None); ax.set_ylim(0, None); ax.grid(alpha=.3); ax.legend()
plt.tight_layout(); plt.savefig("u2_vs_T.png", dpi=200); plt.close()

