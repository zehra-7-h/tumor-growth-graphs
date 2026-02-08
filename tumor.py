# -*- coding: utf-8 -*-
"""
TÜMÖR BÜYÜMESİ MODELİ - LOJİSTİK BÜYÜME DENKLEMİ
"""

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

# ===================== PARAMETRELER =====================
P0 = 1.0          # Başlangıç tümör boyutu (mm^3)
r = 0.5           # Büyüme hızı (1/gün)
K = 100.0         # Maksimum tümör boyutu (taşıma kapasitesi) (mm^3)
t_span = (0, 20)  # Zaman aralığı (gün)
t_eval = np.linspace(0, 20, 200)  # Grafik için zaman noktaları

# ===================== ANALİTİK ÇÖZÜM =====================
def analytic_solution(t, P0, r, K):
    """
    Lojistik denklemin analitik çözümü
    P(t) = K / [1 + ((K-P0)/P0) * exp(-r*t)]
    """
    return K / (1 + ((K - P0) / P0) * np.exp(-r * t))

# ===================== SAYISAL ÇÖZÜM (ODE Çözücü) =====================
def logistic_growth(t, P, r, K):
    """
    Lojistik büyüme diferansiyel denklemi: dP/dt = r*P*(1 - P/K)
    solve_ivp için uygun formatta
    """
    return r * P * (1 - P / K)

# Diferansiyel denklemi çöz
sol = solve_ivp(logistic_growth, t_span, [P0], args=(r, K),
                t_eval=t_eval, method='RK45', rtol=1e-6)

# ===================== GRAFİK 1: Temel Büyüme Eğrisi =====================
fig = plt.figure(figsize=(14, 10))

plt.subplot(2, 2, 1)
# Analitik ve sayısal çözümü karşılaştır
P_analytic = analytic_solution(t_eval, P0, r, K)
plt.plot(t_eval, P_analytic, 'b-', linewidth=2, label='Analitik Çözüm')
plt.plot(sol.t, sol.y[0], 'ro', markersize=4, label='Sayısal Çözüm (RK45)', alpha=0.6)
plt.axhline(y=K, color='g', linestyle='--', label=f'Taşıma Kapasitesi (K={K})')
plt.xlabel('Zaman (gün)', fontsize=11)
plt.ylabel('Tümör Boyutu (mm³)', fontsize=11)
plt.title('Lojistik Tümör Büyümesi - Analitik vs Sayısal Çözüm', fontsize=12, fontweight='bold')
plt.legend(loc='lower right')
plt.grid(True, alpha=0.3)

# ===================== GRAFİK 2: Farklı Büyüme Hızları (r) =====================
plt.subplot(2, 2, 2)
r_values = [0.2, 0.5, 0.8, 1.2]  # Farklı büyüme hızları
for r_val in r_values:
    P_temp = analytic_solution(t_eval, P0, r_val, K)
    plt.plot(t_eval, P_temp, label=f'r = {r_val}', linewidth=2)
plt.xlabel('Zaman (gün)', fontsize=11)
plt.ylabel('Tümör Boyutu (mm³)', fontsize=11)
plt.title('Farklı Büyüme Hızlarının Etkisi', fontsize=12, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)

# ===================== GRAFİK 3: Farklı Başlangıç Boyutları =====================
plt.subplot(2, 2, 3)
P0_values = [0.5, 1.0, 5.0, 20.0]  # Farklı başlangıç boyutları
for P0_val in P0_values:
    P_temp = analytic_solution(t_eval, P0_val, r, K)
    plt.plot(t_eval, P_temp, label=f'P0 = {P0_val}', linewidth=2)
plt.xlabel('Zaman (gün)', fontsize=11)
plt.ylabel('Tümör Boyutu (mm³)', fontsize=11)
plt.title('Farklı Başlangıç Boyutlarının Etkisi', fontsize=12, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)

# ===================== GRAFİK 4: Tedavi Simülasyonu (r değişimi) =====================
plt.subplot(2, 2, 4)
# Tedavi öncesi (yüksek r)
r_before = 0.7
P_before = analytic_solution(t_eval, P0, r_before, K)

# Tedavi sonrası (düşük r) - 10. günde tedavi başlıyor
t_before = t_eval[t_eval <= 10]
t_after = t_eval[t_eval >= 10]
r_after = 0.1  # Tedaviyle büyüme hızı azalıyor

# Tedavi sonrası için yeni başlangıç koşulu
P_at_10 = analytic_solution(10, P0, r_before, K)
P_after = analytic_solution(t_after - 10, P_at_10, r_after, K)

plt.plot(t_before, analytic_solution(t_before, P0, r_before, K),
         'b-', linewidth=3, label='Tedavi Öncesi (r=0.7)')
plt.plot(t_after, P_after, 'r-', linewidth=3, label='Tedavi Sonrası (r=0.1)')
plt.axvline(x=10, color='k', linestyle='--', linewidth=2, label='Tedavi Başlangıcı')
plt.xlabel('Zaman (gün)', fontsize=11)
plt.ylabel('Tümör Boyutu (mm³)', fontsize=11)
plt.title('Kemoterapi Etkisinin Simülasyonu', fontsize=12, fontweight='bold')
plt.legend()
plt.grid(True, alpha=0.3)

plt.suptitle('TÜMÖR BÜYÜMESİ MODELİ - LOJİSTİK DENKLEM ANALİZİ',
             fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()


output_file = "tumor_buyumesi_grafikler.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"✓ Grafikler '{output_file}' dosyasına kaydedildi.")
plt.close(fig)

# ===================== EKRAN ÇIKTILARI =====================
print("\n" + "="*70)
print("TÜMÖR BÜYÜMESİ MODELİ - LOJİSTİK DENKLEM SONUÇLARI")
print("="*70)
print(f"\nPARAMETRELER:")
print(f"  • Başlangıç tümör boyutu (P0): {P0} mm³")
print(f"  • Büyüme hızı (r): {r} 1/gün")
print(f"  • Taşıma kapasitesi (K): {K} mm³")

print(f"\nANALİTİK SONUÇLAR:")
days = [0, 5, 10, 15, 20]
for day in days:
    P = analytic_solution(day, P0, r, K)
    growth_rate = r * P * (1 - P / K)
    percent_of_K = (P / K) * 100
    print(f"  {day:2d}. gün: {P:7.2f} mm³ | Büyüme hızı: {growth_rate:6.3f} mm³/gün | K'nin %{percent_of_K:5.1f}'i")

print(f"\nÖNEMLİ ANALİZ NOKTALARI:")
# En hızlı büyüme anı (P = K/2)
if P0 < K/2:
    # DOĞRU FORMÜL: t = (1/r) * ln((K - P0)/P0)
    t_half = (1/r) * np.log((K - P0)/P0)
    P_half = K/2  # En hızlı büyüme K/2 noktasında
    print(f"  • En hızlı büyüme (P=K/2): {t_half:.2f}. günde, boyut: {P_half:.2f} mm³")
    print(f"  • Maksimum büyüme hızı: {r * P_half * (1 - P_half/K):.3f} mm³/gün")
elif P0 == K/2:
    print(f"  • En hızlı büyüme: Başlangıçta (P0 = K/2 = {K/2} mm³)")
    print(f"  • Maksimum büyüme hızı: {r * P0 * (1 - P0/K):.3f} mm³/gün")
else:
    print(f"  • En hızlı büyüme: Başlangıçta (P0={P0} > K/2={K/2})")
    print(f"  • Büyüme hızı sürekli azalıyor")

# 95% K'ye ulaşma zamanı
if P0 < K * 0.95:
    # DOĞRU FORMÜL: t = (1/r) * ln((K/0.95 - P0)/(P0 * (1/0.95 - 1)))
    t_95 = (1/r) * np.log((K/0.95 - P0)/(P0 * (1/0.95 - 1)))
    print(f"  • %95 kapasiteye ulaşma: {t_95:.2f}. gün")
else:
    print(f"  • Zaten %{P0/K*100:.1f} kapasitede")

print(f"\nTEDAVİ SİMÜLASYONU (10. günde kemoterapi):")
P_no_treatment = analytic_solution(20, P0, 0.7, K)
P_with_treatment = analytic_solution(10, analytic_solution(10, P0, 0.7, K), 0.1, K)
reduction = ((P_no_treatment - P_with_treatment) / P_no_treatment) * 100
print(f"  • Tedavisiz 20. gün: {P_no_treatment:.1f} mm³")
print(f"  • Tedavili 20. gün: {P_with_treatment:.1f} mm³")
print(f"  • Boyut azalması: %{reduction:.1f}")

print(f"\n" + "="*70)
print("GRAFİK DOSYASI: 'tumor_buyumesi_grafikler.png'")
print("Dosyayı görüntülemek için bu dizine bakın:")
print(f"  {os.path.abspath(output_file)}")
print("="*70)

# Ek: Model doğrulama
print(f"\nMODEL DOĞRULAMA (t=10 gün için):")
t_test = 10.0
P_test = analytic_solution(t_test, P0, r, K)
dP_dt_analytic = r * P_test * (1 - P_test / K)
# Sayısal türev
h = 1e-5
P_plus = analytic_solution(t_test + h, P0, r, K)
P_minus = analytic_solution(t_test - h, P0, r, K)
dP_dt_numeric = (P_plus - P_minus) / (2 * h)
print(f"  • Tümör boyutu: {P_test:.4f} mm³")
print(f"  • Analitik dP/dt: {dP_dt_analytic:.6f}")
print(f"  • Sayısal dP/dt: {dP_dt_numeric:.6f}")
print(f"  • Hata: {abs(dP_dt_analytic - dP_dt_numeric):.2e}")