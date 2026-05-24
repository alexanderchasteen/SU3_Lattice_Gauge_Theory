import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("Poly_Loop_Real_Data.csv", delimiter=",")
data2 = np.loadtxt("Poly_Loop_Complex_Data.csv", delimiter=",")
couplings = data1[0, :]


for i in range(len(couplings)):
    real_configs = data1[1:, i] 
    imaginary_configs = data2[1:, i] 
    # normalize = np.sqrt(real_configs*real_configs+imaginary_configs*imaginary_configs)
    plt.scatter(real_configs,imaginary_configs)
    plt.show()





