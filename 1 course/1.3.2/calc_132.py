import numpy as np
import matplotlib.pyplot as plt

g = 9.81
mean = np.mean

def mnk(x, y):
    b = (mean(x*y) - mean(x)*mean(y)) / (mean(x**2) - mean(x)**2)
    N = len(x)
    eps = 1/(b*np.sqrt(N-2)) * np.sqrt((mean(y**2)-mean(y)**2)/(mean(x**2)-mean(x)**2) - b**2)
    a = mean(y) - b*mean(x)
    return b, a, eps

# =====================================================================
# ЧАСТЬ I. Статический метод
# =====================================================================
L   = 148.5
sL  = 0.1
D   = np.mean([10.52, 10.40]); sD = abs(10.52-10.40)/2
d_rod, sd_rod = 0.60, 0.005
l_rod, sl_rod = 175.0, 0.2
s_dl = 0.1

m_load = [50,100,150,200,250,300,400,500]
m_unl1 = [400,300,250,200,150,100,50]
m_unl2 = [400,300,250,200,150,100,50]
S1_load=[-0.1,2.3,5.2,7.45,10.1,12.9,17.6,23.8]; S1_unl=[18.0,12.8,10.3,7.4,5.2,2.6,-0.05]
S2_load=[0.1,2.5,4.9,7.9,10.3,13.1,18,23.4];     S2_unl=[18.1,13.1,9.8,7.9,5.1,2.5,0.2]
S3_load=[0.1,2.3,5.1,7.7,9.9,13,18.5,23.9];      S3_unl=[18,12.8,10.4,7.6,5.1,2.3,-0.1]
runs = [("серия 1, нагружение",m_load,S1_load,"o",True),("серия 1, разгрузка",m_unl1,S1_unl,"o",False),
        ("серия 2, нагружение",m_load,S2_load,"s",True),("серия 2, разгрузка",m_unl2,S2_unl,"s",False),
        ("серия 3, нагружение",m_load,S3_load,"^",True),("серия 3, разгрузка",m_unl2,S3_unl,"^",False)]
phi_all=[]; M_all=[]
for name,m,dl,_,_ in runs:
    m=np.array(m)/1000; dl=np.array(dl)
    M_all += list(m*g*D/100)
    phi_all += list(dl/(2*L))
phi_all=np.array(phi_all); M_all=np.array(M_all); N1=len(phi_all)
f, M0, eps_r1 = mnk(phi_all, M_all)
eps_s1 = np.hypot(sD/D, sL/L)
eps_f = np.hypot(eps_r1, eps_s1)
d_m, l_m = d_rod/100, l_rod/100
G1 = 32*l_m*f/(np.pi*d_m**4)
eps_G1 = np.sqrt(eps_f**2 + (sl_rod/l_rod)**2 + (4*sd_rod/d_rod)**2)
print("== ЧАСТЬ I ==")
print(f"D = {D:.2f} ± {sD:.2f} см; N = {N1}")
print(f"f = {f:.2f} Н·м/рад; eps_случ = {eps_r1*100:.2f}%; eps_сист = {eps_s1*100:.2f}%; eps_f = {eps_f*100:.2f}%; sigma_f = {f*eps_f:.2f}")
print(f"G = {G1/1e9:.1f} ГПа; eps_G = {eps_G1*100:.1f}%; sigma_G = {G1/1e9*eps_G1:.1f} ГПа")
print(f"phi_max = {phi_all.max():.4f} рад, M_max = {M_all.max():.3f} Н·м")

# =====================================================================
# ЧАСТЬ II. Динамический метод
# =====================================================================
m1,m2 = 0.3775, 0.3731
l_w, sl_w = 175.0, 0.2
d_meas = np.array([1.47,1.46,1.47,1.48,1.47,1.48,1.48,1.47,1.48,1.47])
d_w = d_meas.mean(); sd_w = 0.01
sT, sr = 0.02, 0.1

r_cm = np.array([3,4,6,8,9,10.0])
T    = np.array([2.05,2.25,2.80,3.40,3.65,3.90])
x = T**2; y = (r_cm/100)**2
N2 = len(x)
k, a, eps_r2 = mnk(x, y)
eps_s2 = 2*np.sqrt((sT/T.mean())**2 + (sr/r_cm.mean())**2)
eps_k = np.hypot(eps_r2, eps_s2)
f2 = 4*np.pi**2*(m1+m2)*k
dw_m, lw_m = d_w/1000, l_w/100
G2 = 32*lw_m*f2/(np.pi*dw_m**4)
eps_G2 = np.sqrt(eps_k**2 + (sl_w/l_w)**2 + (4*sd_w/d_w)**2)
print("\n== ЧАСТЬ II ==")
print(f"d = {d_w:.3f} мм")
print(f"k = {k*1e3:.3f}e-3 м²/с²; eps_случ = {eps_r2*100:.1f}%; eps_сист = {eps_s2*100:.1f}%; eps_k = {eps_k*100:.1f}%")
print(f"f = {f2*1e2:.3f}e-2 ± {f2*eps_k*1e2:.3f}e-2 Н·м/рад")
print(f"G = {G2/1e9:.1f} ГПа; eps_G = {eps_G2*100:.1f}%; sigma_G = {G2/1e9*eps_G2:.1f} ГПа")

# =====================================================================
# ГРАФИКИ (только вывод на экран)
# =====================================================================
cols = {"серия 1":"tab:blue","серия 2":"tab:green","серия 3":"tab:red"}

fig1, ax1 = plt.subplots(figsize=(7.5,5))
for name,m,dl,mk,filled in runs:
    s_ = name.split(",")[0]
    ph = np.array(dl)/(2*L); MM = np.array(m)/1000*g*D/100
    ax1.errorbar(ph, MM, xerr=s_dl/(2*L), fmt=mk, color=cols[s_], mfc=cols[s_] if filled else "none",
                capsize=2, elinewidth=1, markersize=5, label=name)
xx1 = np.array([phi_all.min(), phi_all.max()])
ax1.plot(xx1, M0+f*xx1, "k-", lw=1.2, label=f"МНК ($f={f:.2f}$ Н·м/рад)")
ax1.set_xlabel(r"$\varphi$, рад"); ax1.set_ylabel(r"$M$, Н$\cdot$м"); ax1.grid(alpha=.3); ax1.legend(fontsize=8)
fig1.tight_layout()

fig2, ax2 = plt.subplots(figsize=(7.5,5))
ax2.errorbar(x, y*1e4, xerr=2*T*sT, yerr=2*r_cm*sr, fmt="o", color="tab:blue", capsize=3, elinewidth=1,
            markersize=5, label="экспериментальные точки")
xx2 = np.array([0, x.max()*1.05]); ax2.plot(xx2, (a+k*xx2)*1e4, "k-", lw=1.2,
            label=f"МНК ($k={k*1e3:.2f}\\cdot10^{{-3}}$ м$^2$/с$^2$)")
ax2.set_xlabel(r"$T^2$, с$^2$"); ax2.set_ylabel(r"$r^2$, см$^2$"); ax2.set_xlim(0,None)
ax2.grid(alpha=.3); ax2.legend(fontsize=9); fig2.tight_layout()

plt.show()