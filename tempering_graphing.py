import pandas as pd
import numpy as np  
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
from scipy.special import erf
from sklearn.metrics import r2_score


SPAT_SIZE=4
TEMP_SIZE=4

df=pd.read_csv('MC'+str(SPAT_SIZE)+'^3*'+str(TEMP_SIZE)+'_analysis.csv')
beta_array=df['Beta'].to_numpy()
avg_plaq=df['Avg Plaq'].to_numpy()
plaq_sus=df['Plaq Sus'].to_numpy()
SD_avg_plaq=df['Jackknife Action SD'].to_numpy()
SD_plaq_sus=df['Jackknife Sus SD'].to_numpy()

df2=pd.read_csv('MC'+str(SPAT_SIZE)+'^3*'+str(TEMP_SIZE)+'_svendsen_analysis.csv')
svendsen_beta_array=df2['Beta'].to_numpy()
svendsen_plaq_array=df2['Avg Plaq'].to_numpy()          
svendsen_sus_array=df2['Plaq Sus'].to_numpy()
svendsen_plaq_SD_array=df2['Jackknife Action SD'].to_numpy()
svendsen_sus_SD_array=df2['Jackknife Sus SD'].to_numpy()


jacknife_data_array = np.loadtxt('MC'+str(SPAT_SIZE)+'^3*'+str(TEMP_SIZE)+'_jackknife_samples.csv', delimiter=',')
beta_array_jacknife=jacknife_data_array[0,:]


plt.rc('font', size=15)
plt.errorbar(beta_array,avg_plaq,yerr=SD_avg_plaq,fmt='o',label='Parallel Tempering Data',color='blue',ecolor='red',elinewidth=3)
plt.errorbar(svendsen_beta_array,svendsen_plaq_array,yerr=svendsen_plaq_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.title('Average Plaquette vs Beta')
plt.xlabel('Beta')
plt.ylabel('Average Plaquette')
plt.show()

plt.errorbar(beta_array,plaq_sus,yerr=SD_plaq_sus,fmt='o',label='Parallel Tempering Data',color='blue',ecolor='red',elinewidth=3)
plt.errorbar(svendsen_beta_array,svendsen_sus_array,yerr=svendsen_sus_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.show()


# # # #Optimal paramters for 6
# beta_max=0.440
# beta_min=0.438
# initial_guess=[0.439, 0.005, 0.04]


# # # Optimal paramters for 5
beta_max=0.4394
beta_min=0.4388
initial_guess=[0.4392, 0.005, 165]

# # #Optimal paramters for 4
# beta_max=0.438
# beta_min=0.436
# initial_guess=[0.438, 0.005, 0.04]

# # # Optimal paramters for 3
# beta_max=0.436
# beta_min=0.428
# initial_guess=[0.430, 0.005, 0.04]

mask = (svendsen_beta_array >= beta_min) & (svendsen_beta_array <= beta_max)

svendsen_beta_array     = svendsen_beta_array[mask]
svendsen_sus_array      = svendsen_sus_array[mask]
svendsen_sus_SD_array  = svendsen_sus_SD_array[mask]
beta_array_jacknife = beta_array_jacknife[mask]
jacknife_data_array = jacknife_data_array[:, mask]


def lorentzian(x, x0, gamma, a):
    return a * gamma**2 / ((x - x0)**2 + gamma**2) 

param,_=curve_fit(lorentzian,svendsen_beta_array, svendsen_sus_array,p0 = initial_guess )
fitA,fitB,fitC=param
R2=r2_score(svendsen_sus_array,lorentzian(svendsen_beta_array,fitA,fitB,fitC))
R2=round(R2,4)
x0=fitA
gamma=fitB
a=fitC




initial_guess=[x0, gamma, a]
cts_beta=np.arange(beta_min,beta_max,(beta_max-beta_min)/1000)
fitY=lorentzian(cts_beta,fitA,fitB,fitC)

plt.rc('font', size=15)
plt.errorbar(svendsen_beta_array,svendsen_sus_array,yerr=svendsen_sus_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.errorbar(beta_array,plaq_sus,yerr=SD_plaq_sus,fmt='o',label='Parallel Tempering Data',color='blue',ecolor='red',elinewidth=3)
plt.plot(cts_beta, fitY, '-', label='Lorentzian Fit R^2='+str(R2),color='blue',alpha=0.5)
plt.legend()
plt.title('Plaquette Susceptibility vs Beta')
plt.xlabel('Beta')
plt.ylabel('Plaquette Susceptibility')  
plt.show()

critical_beta_array = []
FWHM_array = []
peak_height_array = []
for i in range(1, jacknife_data_array.shape[0]):
    # plt.plot(beta_array_jacknife, jacknife_data_array[i,:], color='gray', alpha=0.3)
    # plt.show()
    param,_=curve_fit(lorentzian,beta_array_jacknife,jacknife_data_array[i,:],p0=initial_guess,maxfev=10000 )
    fitA,fitB,fitC=param
    critical_betas = fitA
    FWHM = 2*fitB
    peak_height = fitC
    critical_beta_array.append(critical_betas)
    FWHM_array.append(FWHM)
    peak_height_array.append(peak_height)

def jackknife_variance(data):
    N = len(data)
    mean = np.mean(data)
    return (N - 1) * np.mean((data - mean)**2)


SD_critical_beta=np.sqrt(jackknife_variance(critical_beta_array))
SD_FWHM=np.sqrt(jackknife_variance(FWHM_array))
SD_peak_height=np.sqrt(jackknife_variance(peak_height_array))

print('Critical Beta:', x0, '+/-', SD_critical_beta)
print('FWHM:', 2*gamma, '+/-', SD_FWHM)
print('Peak Height:', a, '+/-', SD_peak_height)
print(str(SPAT_SIZE),',',x0,',',2*gamma,',',a,',', SD_critical_beta,',', SD_FWHM,',', SD_peak_height)