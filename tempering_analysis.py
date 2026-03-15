import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

SIZE=3
Nplaq=SIZE**4*6
df=pd.read_csv('Raw_MC_Data_M3.csv')


# df=pd.read_csv('/home/alexa/Z2_Lattice/Z2_Lattice_PT/PT'+str(SIZE)+'_RawMC.csv')
beta_array=df.columns.to_numpy()
beta_array=beta_array.astype(float)
CONFIGS=len(beta_array)
beta_nbhd_radius=(beta_array[1]-beta_array[0])/2


# Used on the blocked values to compute variances on quantities. Can resample wrt the mean the blocks which have the data in the right form (i.e computed expectations on subpopulations and doing statistics on the subpopulations of measurments)
def jacknife_var_on_mean(array):
    n=len(array)
    mean=np.mean(array)
    jacknife_means=[]
    for i in range(n):
        jacknife_sample=np.delete(array,i)
        jacknife_mean=np.mean(jacknife_sample)
        jacknife_means.append(jacknife_mean)
    jacknife_means=np.array(jacknife_means)
    variance=(n-1)/n * np.sum((jacknife_means - mean)**2)
    return variance

# block the data and compute statistical functions on each block, like mean, variance, stuff...
def block_via_function(array,blocks,function):
    size_of_array=len(array)
    block_size=int(size_of_array/blocks)
    blocked_array=[]
    for id in range(blocks):
        start=id*block_size
        end=start+block_size
        block=array[start:end]
        block_stat=function(block)
        blocked_array.append(block_stat)
    return np.array(blocked_array)

# This function will let us block our data set into N blocks and each block we will do the reweighting procedure. 
# The function will output the expectation value of the array. Its a pain in the ass to generalize it for any functon of observable (kind of, just not worth the thought as of now)
# So for our purposeses we want to apply this to array^2 and array
def svendsen_blocking_SPECIAL(array,blocks,weights):
    size_of_array=len(array)
    block_size=int(size_of_array/blocks)
    blocked_array=[]
    for id in range(blocks):
        start=id*block_size
        end=start+block_size
        block=np.array(array[start:end])
        block_weights=np.array(weights[start:end])
        block = block*block_weights
        normalize=np.sum(block_weights)
        block_stat=np.sum(block)/normalize
        blocked_array.append(block_stat)
    return np.array(blocked_array)


# returns the jacknife samples, for sus, so we can curve fit and then do statistics on each fit
def jacknife_samples(array):
    jacknife_samples=[]
    for i in range(len(array)):
        v=0
        for j in range(len(array)):
            if i==j:
                v+=0
            else:
                v+=array[j]
        v=v/(len(array)-1)
        jacknife_samples.append(v)
    return jacknife_samples

def fourth_moment(array):
    m = np.mean(array)
    s=0.0
    for x in array:
        diff = x-m
        s += diff*diff*diff*diff
    s /= len(array)
    return s


def binder_cum(array):
    fourth_mom = fourth_moment(array)
    second_mom = np.var(array)
    binder_cum = 1-(fourth_mom)/(3*second_mom*second_mom)
    return binder_cum


def svendsen_binder_block(array, blocks, weights):
    eps=1e-14
    size = len(array)
    block_size = size // blocks
    binder_blocks = []
    for i in range(blocks):
        start = i * block_size
        end   = start + block_size
        x = array[start:end]
        w = weights[start:end]
        wsum = np.sum(w)
        if wsum == 0:
            binder_blocks.append(np.nan)
            continue
        mean = np.sum(w * x) / wsum
        diff = x - mean
        mu2 = np.sum(w * diff**2) / wsum
        if mu2 < eps:
            binder_blocks.append(np.nan)
            continue

        mu4 = np.sum(w * diff**4) / wsum
        binder_blocks.append(1.0 - mu4 / (3.0 * mu2 * mu2))
    return np.array(binder_blocks)


blocks=50
avg_plaq=[]
SD_avg_plaq=[]
svendsen_plaq_array=[]
svendsen_plaq_SD_array=[]

plaq_sus=[]
SD_plaq_sus=[]
svendsen_sus_array=[]
svendsen_sus_SD_array=[]

binder = []
SD_Binder = []
svendsen_binder = []
svensen_binder_SD_array =[]

svendsen_beta_array=[]

all_jackknife_samples = []

for i in range(CONFIGS):
    array=df.iloc[:,i].to_numpy()

    # block the data according to the obsverable we want to measure
    blocked_plaq = block_via_function(array,blocks,function = np.mean)
    avg_plaq.append(np.mean(blocked_plaq))
    SD_avg_plaq.append(np.sqrt(jacknife_var_on_mean(blocked_plaq)))

    blocked_plaq_sus = Nplaq*block_via_function(array,blocks,function = np.var)
    plaq_sus.append(np.mean(blocked_plaq_sus))
    SD_plaq_sus.append(np.sqrt(jacknife_var_on_mean(blocked_plaq_sus)))

    blocked_binder = block_via_function(array,blocks,function = binder_cum)
    binder.append(np.mean(blocked_binder))
    SD_Binder.append(np.sqrt(jacknife_var_on_mean(blocked_binder)))

    svendsen_length=500
    for k in np.linspace(-beta_nbhd_radius, beta_nbhd_radius, 1000):
        new_beta=beta_array[i]+k
        svendsen_beta_array.append(new_beta)
        # compute weights
        delta_beta=k
        logw_array = delta_beta * Nplaq * array
        maxlogw = np.max(logw_array)
        w_array = np.exp(np.clip(logw_array - maxlogw, -700, 700))  # prevent overflow
       
        # block the data corresponing to the weights. Compute statistics on new weighted data, representing the informaton at beta'
        svendsen_block_of_P = svendsen_blocking_SPECIAL(array,blocks,w_array)
        avg_plq_svendsen=np.mean(svendsen_block_of_P)
        svendsen_plaq_array.append(avg_plq_svendsen)
        svendsen_plaq_SD_array.append(np.sqrt(jacknife_var_on_mean(svendsen_block_of_P)))

        
        svendsen_block_of_P2 = svendsen_blocking_SPECIAL(array**2, blocks, w_array)
        svendsen_block_of_sus =  Nplaq * (svendsen_block_of_P2 - svendsen_block_of_P**2)
        plaq_sus_svendsen = np.mean(svendsen_block_of_sus)
        svendsen_sus_array.append(plaq_sus_svendsen)  
        svendsen_sus_SD_array.append(np.sqrt(jacknife_var_on_mean(svendsen_block_of_sus)))

        # svendsen_block_of_P3 = svendsen_blocking_SPECIAL(array**3, blocks, w_array)
        # svendsen_block_of_P4 = svendsen_blocking_SPECIAL(array**4, blocks, w_array)
        # svendsen_fourth_mom = svendsen_block_of_P4-4*svendsen_block_of_P*svendsen_block_of_P3+6*svendsen_block_of_P2*svendsen_block_of_P**2-3*svendsen_block_of_P**4
        # denom = (svendsen_block_of_P2 - svendsen_block_of_P**2)
        # binder_block_svendsen =1-svendsen_fourth_mom/(3*denom*denom)
        # svendsen_binder_cum = np.mean(binder_block_svendsen)
        # svendsen_binder.append(svendsen_binder_cum)
        # svensen_binder_SD_array.append(np.sqrt(jacknife_var_on_mean(binder_block_svendsen)))

        binder_block_svendsen = svendsen_binder_block(array, blocks, w_array)
        svendsen_binder.append(np.nanmean(binder_block_svendsen))
        svensen_binder_SD_array.append(np.sqrt(jacknife_var_on_mean(binder_block_svendsen[~np.isnan(binder_block_svendsen)])))

       

        jacknife_samples_at_svendesen_beta=jacknife_samples(svendsen_block_of_sus)
        all_jackknife_samples.append(jacknife_samples_at_svendesen_beta)
      
plt.errorbar(beta_array,avg_plaq,yerr=SD_avg_plaq,fmt='o',label='Parallel Tempering Data',color='blue',ecolor='red',elinewidth=3)
plt.errorbar(svendsen_beta_array,svendsen_plaq_array,yerr=svendsen_plaq_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.show()
plt.errorbar(beta_array,plaq_sus,yerr=SD_plaq_sus,fmt='o',label='Parallel Tempering Data',color='blue',ecolor='red',elinewidth=3)
plt.errorbar(svendsen_beta_array,svendsen_sus_array,yerr=svendsen_sus_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.show()
plt.errorbar(svendsen_beta_array,svendsen_binder,yerr=svensen_binder_SD_array,fmt='o',label='Svendsen Reweighted Data',color='blue',ecolor='yellow',elinewidth=3,capsize=0,markersize=1)
plt.show()


folder = './Analyzed_Data/'
montecarlo_data=np.transpose(np.array([beta_array,avg_plaq,plaq_sus,SD_avg_plaq,SD_plaq_sus]))
filename=folder + 'PT'+str(SIZE)+'_analysis.csv'

header = ["Beta", "Avg Plaq", "Plaq Sus", "Jackknife Action SD", "Jackknife Sus SD"]
with open(filename, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csv_writer = csv.writer(csvfile)

    # Write all rows at once
    csv_writer.writerow(header)
    csv_writer.writerows(montecarlo_data)


print(f"CSV file '{filename}' created successfully.")


svendsen_data=np.transpose(np.array([svendsen_beta_array,svendsen_plaq_array,svendsen_sus_array,svendsen_binder,svendsen_plaq_SD_array,svendsen_sus_SD_array,svensen_binder_SD_array]))
filename2=folder + 'PT'+str(SIZE)+'_svendsen_analysis.csv'

header = ["Beta", "Avg Plaq", "Plaq Sus","Binder Cum", "Jackknife Action SD", "Jackknife Sus SD","Binder Cum SD"]
with open(filename2, 'w', newline='') as csvfile:
    # Create a CSV writer object
    csv_writer = csv.writer(csvfile)

    # Write all rows at once
    csv_writer.writerow(header)
    csv_writer.writerows(svendsen_data)


print(f"CSV file '{filename2}' created successfully.")



jackknife_filename =folder + 'PT'+str(SIZE)+'_jackknife_samples.csv'


jackknife_matrix = np.transpose(np.array(all_jackknife_samples))

# Create header row (beta values)
header = [f"{beta:.6f}" for beta in svendsen_beta_array]

with open(jackknife_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    writer.writerows(jackknife_matrix)

print(f"CSV file '{jackknife_filename}' created successfully with betas as columns.")


