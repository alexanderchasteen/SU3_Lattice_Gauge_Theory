import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

a=1
beta=2

samples = pd.read_csv('output.csv').to_numpy()

counts, bin_edges = np.histogram(samples, bins=100, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.bar(bin_centers, counts, width=bin_edges[1]-bin_edges[0], alpha=0.6, label='Sampled')

# Plot theoretical distribution (normalized)
x_theory = np.linspace(-1, 1, 500)
p_theory = np.sqrt(1 - x_theory**2) * np.exp(a*beta*x_theory)

print( np.trapz(p_theory, x_theory))
p_theory /= np.trapz(p_theory, x_theory)  # normalize area to 1
plt.plot(x_theory, p_theory, 'r', lw=2, label='Theory')

plt.xlabel('x0')
plt.ylabel('Probability density')
plt.legend()
plt.show()
plt.show()