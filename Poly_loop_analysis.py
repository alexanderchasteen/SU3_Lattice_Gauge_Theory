import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N_T = 6.0

# --- Load data ---
data = np.loadtxt("Poly_Loop_Data.csv", delimiter=",")
r = data[0, :]
configs = data[1:, :]  # shape (Ncfg, Nr)

# --- Compute mean Polyakov loop ---
means = np.mean(configs, axis=0)

# --- Propagate standard deviation through log ---
P_std = np.std(configs, axis=0, ddof=1)
V = -np.log(means) / N_T
V_err = P_std / (N_T * means)

# --- Cornell potential model ---
def potential(r, A, B, sigma):
    return A + B / r + sigma * r

# --- Fit using propagated errors ---
params, cov = curve_fit(potential, r, V, sigma=V_err, absolute_sigma=True)
A_fit, B_fit, sigma_fit = params
errors = np.sqrt(np.diag(cov))
A_err, B_err, sigma_err = errors

print("\nFitted parameters (propagated std):")
print(f"  A      = {A_fit:.4f} ± {A_err:.4f}")
print(f"  B      = {B_fit:.4f} ± {B_err:.4f}")
print(f"  sigma  = {sigma_fit:.4f} ± {sigma_err:.4f}")

# --- Plot 1: With propagated error bars ---
plt.figure(figsize=(7,4))
plt.title("Polyakov Loop Potential (with propagated error)")
plt.errorbar(r, V, yerr=V_err, fmt='o', markersize=3, capsize=2,
             alpha=0.7, label="Data (propagated std)")
plt.plot(r, potential(r, *params), "-", lw=2, label=f"Fit σ={sigma_fit:.3f}")
plt.xlabel("Distance r (lattice units)")
plt.ylabel("Potential V(r)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --- Plot 2: Without error bars ---
plt.figure(figsize=(7,4))
plt.title("Polyakov Loop Potential (mean only)")
plt.plot(r, V, "o", markersize=3, alpha=0.7, label="Mean potential")
plt.plot(r, potential(r, *params), "-", lw=2, label=f"Fit σ={sigma_fit:.3f}")
plt.xlabel("Distance r (lattice units)")
plt.ylabel("Potential V(r)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
