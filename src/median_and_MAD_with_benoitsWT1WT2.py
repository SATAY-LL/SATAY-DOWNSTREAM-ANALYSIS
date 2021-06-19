# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 15:57:01 2021

@author: floor
"""

import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
import pylab
import statistics
from statsmodels import robust

#in this code, a method called Median Absolute Deviation is used to remove outliers from the reads per gene.
#The input is data from yeastmine and from earlier analysis steps:
# - ChromosomeRegion_AllGenes.tsv
# - data_step2_readsWT1.csv
# - data_step2_readsWT2.csv
# - data_step2WT1.csv
# - data_step2WT2.csv

dataframe = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\data\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
dataframe['readsWT1'] = ''
dataframe['reads_avr1'] = ''
dataframe['readsWT2'] = ''
dataframe['reads_avr2'] = ''
#dataframe['error1'] = ''
#dataframe['error2'] = ''
dataframe['fitnessWT1'] = ''
dataframe['fitnessWT2'] = ''
dataframe['fitness_avr1'] =  ''
dataframe['fitness_avr2'] =  ''
dataframe['fitness_avr2median'] = ''
dataframe['fitness_avr1median'] = ''
dataframe['reads_avr2MAD'] = ''
dataframe['reads_avr1MAD']= ''

# import readsWT1
dataframe['readsWT1']= dataframe['readsWT1'].astype('object')

with open('data\data_step2_readsWT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['readsWT1'] = data

#import reads WT2
dataframe['readsWT2']= dataframe['readsWT2'].astype('object')

with open('data\data_step2_readsWT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data2 = list(reader)
dataframe['readsWT2'] = data2

#import fitness WT1
dataframe['fitnessWT1']= dataframe['fitnessWT1'].astype('object')

with open('data\data_step2WT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data3 = list(reader)
dataframe['fitnessWT1'] = data3

#import fitness WT2
dataframe['fitnessWT2']= dataframe['fitnessWT2'].astype('object')

with open('data\data_step2WT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data4 = list(reader)
dataframe['fitnessWT2'] = data4

del(data, data2,data3,data4)
z = 1.960 #for 95 percent confidence interval
for i in range(0,len(dataframe)):


    if (dataframe.at[i,'fitness_avr1']) == [] or len(dataframe.at[i, 'fitnessWT1']) < 6:
        dataframe.at[i,'reads_avr1'] =  ''
        dataframe.at[i,'fitness_avr1'] =  ''
        #dataframe.at[i, 'error1'] = 0

    else:
        dataframe.at[i,'reads_avr1'] = statistics.median(dataframe.at[i,'readsWT1'])
        
        #using the median absolute deviation, remove outliers
        MAD = robust.mad(dataframe.at[i,'readsWT1']) # MAD with b = 1.4826, so its assumed to be normal
        minimal = dataframe.at[i, 'reads_avr1']-2.5*MAD
        maximal = dataframe.at[i, 'reads_avr1']+2.5*MAD
        dataframe.at[i, 'reads_avr1'] =  statistics.median(dataframe.at[i,'readsWT1'])
        readspergene = np.array(dataframe.at[i,'readsWT1'])
        readspergene = readspergene[(readspergene>=minimal) & (readspergene<= maximal)]   
        dataframe.at[i, 'readsWT1'] = readspergene  
        dataframe.at[i,'reads_avr1MAD'] = statistics.median(dataframe.at[i,'readsWT1'])

        dataframe.at[i,'fitness_avr1'] = statistics.mean(dataframe.at[i,'fitnessWT1'])
        dataframe.at[i,'fitness_avr1median'] = statistics.median(dataframe.at[i,'fitnessWT1'])
        #dataframe.at[i, 'error1'] =  math.log(np.longdouble(2*(z*(statistics.stdev(dataframe.at[i, 'readsWT1']))/(np.sqrt(len(dataframe.at[i,'readsWT1'])))))*N/n1*a, (2))/gfHO1
    
    if (dataframe.at[i,'fitness_avr2']) == [] or len(dataframe.at[i, 'fitnessWT2']) < 6:
        dataframe.at[i,'reads_avr2'] =  ''
        dataframe.at[i,'fitness_avr2'] =  ''
        #dataframe.at[i, 'error2'] = 0

    else:
        
        dataframe.at[i,'reads_avr2'] = statistics.median(dataframe.at[i,'readsWT2'])
        
        #using the median absolute deviation, remove outliers
        MAD = robust.mad(dataframe.at[i,'readsWT2']) # MAD with b = 1.4826, so its assumed to be normal
        minimal = dataframe.at[i, 'reads_avr2']-2.5*MAD
        maximal = dataframe.at[i, 'reads_avr2']+2.5*MAD
        dataframe.at[i, 'reads_avr2'] =  statistics.median(dataframe.at[i,'readsWT2'])
        readspergene = np.array(dataframe.at[i,'readsWT2'])
        readspergene = readspergene[(readspergene>=minimal) & (readspergene<= maximal)]   
        dataframe.at[i, 'readsWT2'] = readspergene  
        dataframe.at[i,'reads_avr2MAD'] = statistics.median(dataframe.at[i,'readsWT2'])

        dataframe.at[i,'fitness_avr2'] = statistics.mean(dataframe.at[i,'fitnessWT2'])
        dataframe.at[i,'fitness_avr2median'] = statistics.median(dataframe.at[i,'fitnessWT2'])
        #dataframe.at[i, 'error2'] =  math.log(np.longdouble(2*(z*(statistics.stdev(dataframe.at[i, 'readsWT2']))/(np.sqrt(len(dataframe.at[i,'readsWT2'])))))*N/n2*a, (2))/gfHO2

##filter on values that you want to see 
readsWT1 = []
readsWT2 = []
readsWT1MAD = []
readsWT2MAD =[] 
#err1 = []
#err2 = []
WT1 =[]
WT2 = []
WT1med = []
WT2med = []
chromosome  = 'chrXII'
for i in range(0,len(dataframe)):
    if  type(dataframe.at[i,'fitness_avr1']) == float and  type(dataframe.at[i,'fitness_avr2']) == float  and len(dataframe.at[i,'fitnessWT1']) >5 and len(dataframe.at[i,'fitnessWT2']) >5 and dataframe.at[i,'reads_avr1']<500 and dataframe.at[i,'reads_avr2'] <500 :
        readsWT1.append(dataframe.at[i,'reads_avr1'])
        readsWT2.append((dataframe.at[i,'reads_avr2']))
        readsWT2MAD.append(dataframe.at[i,'reads_avr2MAD'])
        readsWT1MAD.append(dataframe.at[i,'reads_avr1MAD'])
#        err1.append((dataframe.at[i,'error1']))
#        err2.append((dataframe.at[i,'error2']))
        WT1.append(dataframe.at[i,'fitness_avr1'])
        WT2.append(dataframe.at[i,'fitness_avr2'])
        WT1med.append(dataframe.at[i,'fitness_avr1median'])
        WT2med.append(dataframe.at[i,'fitness_avr2median'])
        
#and dataframe.at[i, 'chromosome'] == chromosome 
        
#compute mean square error
a = np.array(WT1) 
b = np.array(WT2) 
mses = np.sqrt(((a-b)**2).mean())

aa = np.array(WT1med) 
bb = np.array(WT2med) 
msesmed = np.sqrt(((aa-bb)**2).mean())

A = np.array(readsWT1) 
B = np.array(readsWT2) 
msesreads = np.sqrt(((A-B)**2).mean())

AA = np.array(readsWT1MAD) 
BB = np.array(readsWT2MAD) 
msesreadsMAD = np.sqrt(((AA-BB)**2).mean())

print('mean square error of fitness using the mean is', mses)
print('MSE of fitness using median is', msesmed)
print('MSE of reads (using median) is', msesreads)
print('MSE of reads with median absolute deviation is', msesreadsMAD)

 ## PLOTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#plot reads wt1 wt2     
plt.style.use('Solarize_Light2')
plt.plot(readsWT1, readsWT2,  'x', color='black')
#plt.errorbar(WT1, WT2,  xerr=err1, yerr=err2)
plt.xlabel('reads WT1')
plt.ylabel('reads WT2')  
plt.title('Reads WT1 and WT2, ' )#(+str(chromosome))
plt.axline([0, 0], [1, 1], label = 'x=y')

# calculate the trendline
z = np.polyfit(readsWT2, readsWT1, 1)
p = np.poly1d(z)
pylab.plot(readsWT1,p(readsWT1),"r--", label = 'fit of scatter')

plt.legend()
plt.savefig('afb.jpg', dpi=1200)
plt.show()

#plot fitness WT1 and WT2
plt.figure ()   
plt.plot(WT1, WT2, 'x', color = 'black', label = 'mean, MSE = ' +str(mses))
plt.plot(WT1med, WT2med, 'x', color = 'red', label = 'median, MSE=' +str(msesmed))
plt.xlabel('fitness WT1')
plt.ylabel('fitness WT2')  
plt.title('Fitness WT1 and WT2 ' )#(+str(chromosome))
plt.axline([0, 0], [1, 1], label = 'x=y')

# calculate the trendline
z = np.polyfit(WT2, WT1, 1)
p = np.poly1d(z)
pylab.plot(WT1,p(WT1),"r--", label = 'fit of scatter')
plt.legend()
plt.show()
