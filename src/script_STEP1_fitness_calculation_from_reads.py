# -*- coding: utf-8 -*-
"""
Created on Thu May  6 16:15:48 2021

@author: floor
"""
import numpy as np
import math 
import pandas as pd
import wiggelen
import matplotlib.pyplot as plt
import statistics
import statsmodels.stats.proportion as sm
from def_portion_total_reads import portion_of_allreads


# STEP 1:
# Calculate fitness for each insertion using reads and the HO locus

#the function to calculate a fitness value. This function is used in step II. 
#filepath_and_name = r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam.wig'
def fitness1 (filepath_and_name): 
    data = portion_of_allreads(filepath_and_name, 0.4)
                                      
    print('Start STEP1')                                 
      # #HO locus is neutral, compare to GF of this locus>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      # #first find avarage reads of HO locus, then compute growthfactor. Can also be done other way around.
    # dataHO =  data[ 
    #   (data['chromosome']  == 'chrIV') &
    #   (46270 < data['tn start position']) & (data['tn start position'] <  48031) 
    #   ]
    # readsHO = dataHO['#reads'].mean()  
    # print(readsHO)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    N=117601914802*10000000 #number of transposons before sequencing (())
    
      #use imput to find n and m
    n =data['50% reads'].sum()      #total number of reads#m = len(ni) #number of different type transposons
    m =len(data) #number of different transposons
    
    ni = data['50% reads'].to_numpy() #array with reads for every insertion
    Ni = np.zeros(m)
    Ci = np.zeros(m)
    gf = np.zeros(m)
    fitness = np.zeros(m)
    gfHO = []
    pcr = 35  #number of PCR cycles
    a = 2**-pcr
    data['confidence_levels'] = ''
    
    #forloop to compute growthfactor and fitness of insertions, relative to HO locus 
    for i in range (0,m):
        
        Ni[i] = np.longdouble(ni[i]) *N/n #calculate number of each tn present before sequencing
        Ci[i] = Ni[i] * a  # calculate number of tn present before PCR
        
        if Ci[i] == 0:
            gf[i] = 0
        else: 
            gf[i] = math.log(Ci[i],(2)) #calculate the growth factor 
    #delete ADE2 from data:
        if  (data.at[i,'chromosome']  == 'chrXV') & (564476 < data.at[i,'tn start position']) & (data.at[i,'tn start position'] <  566191):
             data = data.drop([i], axis=0) 
    gfHO = statistics.mean(gf) #use this line if average of HO locus is taken after calculations of gf 
    #gfHO = math.log(np.longdouble(readsHO)*N/n*a, (2)) #use this line if average of reads HO locus is calculated above
    data = data.reset_index()
    for i in range (0,m):
        fitness [i] = (gf[i]/gfHO) #define fitness relative to HO locus
    # for i in range(0, len(data)):
    #     data.at[i,'confidence_levels'] = sm.proportion_confint(data.at[i,'#reads'], n, alpha=0.05, method='agresti_coull')
    data = data.join(pd.DataFrame(fitness))
    #join(pd.DataFrame(fitness))
    data.rename(columns = {0:'fitness'}, inplace = True)
    #for i in range (0,len(data)):
    #    if data.at[i,'fitness']>3.5:
    #         print(data.at[i,'chromosome'], data.at[i,'tn start position'])
    return(data)
#data = pd.DataFrame(fitness1(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-1_merged-techrep-a_techrep-b_trimmed.sorted.bam.wig'))

#plot >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# plt.plot(data['fitness'], 'o')
# plt.title('Fitness relative to the HO locus')
# plt.style.use('seaborn') 
# plt.xlabel('start position transposon')
# plt.ylabel('relative fitness')
# plt.rcParams.update({'font.size': 22})
# plt.show()
# #plt.errorbar(data['index'], (data['#reads']), yerr = data['confidence_levels'])
# plt.savefig('plots/afb', dpi=800)

   
