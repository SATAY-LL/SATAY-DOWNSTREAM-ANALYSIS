# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 12:31:52 2021

@author: floor
"""
## import modules
from script_STEP1_fitness_calculation_from_reads import fitness1
import pandas as pd
import numpy as np
import statistics 
import matplotlib as plt
from collections import Counter
from scipy.stats import norm
import scipy

#function for fitness per domain.
#input = wig file
#output= dataframe with fitness values per domain

def fitness_perdomain (filepath_and_name):
    #filepath_and_name = r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam.wig"
    domain_start_stop = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\AGR_Yeast_Genes (2).tsv", sep = "\t", names = ["gene", "start bp gene", "stop bp gene", "start bp domain", "stop bp domain" , 'domain', "chromosome"] )
    
    print('begin step I')
    DataTabel, n = fitness1(filepath_and_name)
    print('end step I')
    
    # sort on chromosomes and add collumn to inject the fitness values
    domain_start_stop['fitness'] = ''
    domain_start_stop['fitness']= domain_start_stop['fitness'].astype('object')
    domain_start_stop = domain_start_stop.sort_values('chromosome')
    domain_start_stop = domain_start_stop.reset_index()
    
    
    ## CHROMOSOME I
    for i in range(0,len(DataTabel)):
        print(i)
        
        if DataTabel.at[i,'chromosome'] == 'chrI':        
            print('chr1')
            
            for j in range (0, 720):           
                if DataTabel.at[i, 'tn start position'] <= (domain_start_stop.at[j, 'stop bp gene']) and DataTabel.at[i, 'tn start position'] >= (domain_start_stop.at[j, 'start bp gene']):
                    TNposition  =DataTabel.at[i, 'tn start position'] - (domain_start_stop.at[j, 'start bp gene'])
                    print(j)
                    
                    if TNposition <= domain_start_stop.at[j, 'stop bp domain'] and TNposition >= domain_start_stop.at[j, 'start bp domain']:
                        
                        a = DataTabel.at[i, 'fitness']
                        a.astype(np.float64)
                        print(a)
                        domain_start_stop.at[j, 'fitness'] = a
                    
                        break
                    break
                   
        
    ## CHROMOSOME II
        if DataTabel.at[i,'chromosome'] == 'chrII':        
            print('chr2')
            
            for j in range (721, 3886):           
                if DataTabel.at[i, 'tn start position'] <= (domain_start_stop.at[j, 'stop bp gene']) and DataTabel.at[i, 'tn start position'] >= (domain_start_stop.at[j, 'start bp gene']):
                    TNposition  =DataTabel.at[i, 'tn start position'] - (domain_start_stop.at[j, 'start bp gene'])
                    print(j)
                    
                    if TNposition <= domain_start_stop.at[j, 'stop bp domain'] and TNposition >= domain_start_stop.at[j, 'start bp domain']:
                        
                        a = DataTabel.at[i, 'fitness']
                        a.astype(np.float64)
                        print(a)
                        domain_start_stop.at[j, 'fitness'] = a
                    
                        break
                    break
        
        ## STOP after two chromosomes (this doesnt work?)
        if DataTabel.at[i,'chromosome'] == 'chrIII':
            i=len(DataTabel) + 1
    return(domain_start_stop)

##compute fitness per domain for NRP1_1, both replicates
filepath_and_name1A = r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-1_merged-DpnII-NlaIII-a_trimmed.sorted.bam.wig"
filepath_and_name1B = r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-1_merged-DpnII-NlaIII-b_trimmed.sorted.bam.wig"

DataTabel1a = fitness_perdomain(filepath_and_name1A)
DataTabel1b = fitness_perdomain(filepath_and_name1B)

## average fitness per domain
DataTabel1a['fitness_avr1a'] =''
DataTabel1b['fitness_avr1b'] =''

for I in range(0, len(DataTabel1a)-1):
    if DataTabel1a.at[I, 'fitness'] != []:
        DataTabel1a.at[I, 'fitness_avr1a'] =  statistics.mean(DataTabel1a.at[I, 'fitness1a'])


for J in range(0, len(DataTabel1b)-1):
    if DataTabel1b.at[J, 'fitness1b'] != []:
     DataTabel1b.at[J,'fitness_avr1b'] =  statistics.mean(DataTabel1b.at[J, 'fitness1b'])


## select for domains that have fitness values in both datatabel
fitness1a =[]
fitness1b = []
difference = []
difference_notrounded = []
#genes10 = []
chromosome  = 'chrI'
for i in range(0,len(DataTabel1a)):
    if  type(DataTabel1a.at[i,'fitness_avr1a']) == float and  type(DataTabel1b.at[i,'fitness_avr1b']) == float  and len(DataTabel1a.at[i,'fitness1a']) >5 and len(DataTabel1b.at[i,'fitness1b']) >5: #and tabel.at[i, 'chromosome'] == chromosome:
       
        fitness1a.append(DataTabel1a.at[i,'fitness_avr1a'])
        fitness1b.append(DataTabel1b.at[i,'fitness_avr1b'])
        #genes10.append(tabel.at[i, 'gene'])
        
## plot to have a first look
plt.figure()
plt.style.use('Solarize_Light2')
plt.scatter(fitness1a, fitness1b,  c='black', alpha = 0.6, marker="x")


## compute difference between two datasets
difference = []
difference_notrounded = []
for k in range(0,len(fitness1a)-1):
    x= (fitness1a[k]-fitness1b[k])
    difference.append(round( x, 2 ) )
    difference_notrounded.append(x)
abundance =Counter(difference)
mu, std = norm.fit(difference_notrounded) 


## plot differences
lists = sorted(abundance.items()) # sorted by key, return a list of tuples
xdata,ydata = zip(*lists) # unpack a list of pairs into two tuples

plt.figure()
plt.title('Deviation between domains of Nrp1-1a and Nrp1-1b')
plt.xlabel('Absolute difference between datasets (per domain)')
plt.ylabel("Abundance")
plt.scatter(xdata, ydata, label= 'Difference in fitness per domain, STD='+ str(round(std, 3)),alpha=0.8, color='orange')
spline = scipy.interpolate.UnivariateSpline(xdata, ydata)
spline.set_smoothing_factor(11500)
#plt.plot(xdata, spline(xdata), 'green', lw=3)
halvemax= (spline(0)/2)
FWHM_points = scipy.interpolate.UnivariateSpline(xdata, ydata - halvemax, s=11500).roots()
x_valuesHWHM = [FWHM_points[1], FWHM_points[0]]
y_valueHWHM = [halvemax, halvemax]
FWHM = FWHM_points[1] -FWHM_points[0]
#plt.plot(x_valuesHWHM, y_valueHWHM, label = 'FWHM (=' +str(round(FWHM,3)) + ')', color = 'black')
# plt.ylim(0, 200)
# plt.xlim(-1, 1)
plt.legend()
plt.savefig('domains', bbox_inches='tight', dpi=1000)
