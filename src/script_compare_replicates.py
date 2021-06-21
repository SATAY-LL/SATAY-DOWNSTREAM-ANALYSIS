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
import math
from tkinter import Tk
from tkinter import filedialog

        
        
#%% Selecting all datasets from your local files 

root = Tk()
root.filename =  filedialog.askopenfilename(title = "Chromosome regions file",filetypes = (("tsv files","*.tsv"),("all files","*.*")))
filename_chrom=root.filename
root.withdraw()        
#%%

#STEP III:
#compare our fitness values to cellmap data
#get data from step II and put it in a dataframe

##### you should change the variable name to something more appropiate for this datafile, otherwise you have to comment
## what it does 
dataframe = pd.read_csv(filename_chrom, sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )

dataframe['readsWT1'] = ''
dataframe['reads_avr1'] = ''
dataframe['readsWT2'] = ''
dataframe['reads_avr2'] = ''
dataframe['error1'] = ''
dataframe['error2'] = ''
dataframe['fitnessWT1'] = ''
dataframe['fitnessWT2'] = ''
dataframe['fitness_avr1'] =  ''
dataframe['fitness_avr2'] =  ''
dataframe['fitness_avr2median'] = ''
dataframe['fitness_avr1median'] = ''
# import readsWT1
dataframe['readsWT1']= dataframe['readsWT1'].astype('object')

## You should do all the imports before to have at the beginning a clear overview of all the data you will be needing 
with open('data_step2_readsWT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['readsWT1'] = data

#import reads WT2
dataframe['readsWT2']= dataframe['readsWT2'].astype('object')

with open('data_step2_readsWT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data2 = list(reader)
dataframe['readsWT2'] = data2

#import fitness WT1
dataframe['fitnessWT1']= dataframe['fitnessWT1'].astype('object')

with open('data_step2WT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data3 = list(reader)
dataframe['fitnessWT1'] = data3

#import fitness WT2
dataframe['fitnessWT2']= dataframe['fitnessWT2'].astype('object')

with open('data_step2WT2.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data4 = list(reader)
dataframe['fitnessWT2'] = data4

del(data,data2,data3,data4)

#%%% what all these variables mean ??? 
z = 1.960 #for 95 percent confidence interval
N=117601914802*10000000 
pcr = 35
a = 2**pcr
n1 = 31793901
n2 = 15302907
gfHO1 = 3.88665
gfHO2 = 4.10502

#%% For loop for? 

for i in range(0,len(dataframe)):


    if (dataframe.at[i,'readsWT1']) == [] or len(dataframe.at[i, 'readsWT1']) < 6 or dataframe.at[i, 'readsWT1'] == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]:
        dataframe.at[i,'reads_avr1'] =  ''
        dataframe.at[i,'fitness_avr1'] =  ''
        dataframe.at[i, 'error1'] = 0

        
    else:
        dataframe.at[i,'reads_avr1'] = statistics.median(dataframe.at[i,'readsWT1'])
        dataframe.at[i,'fitness_avr1'] = statistics.mean(dataframe.at[i,'fitnessWT1'])
        dataframe.at[i,'fitness_avr1median'] = statistics.median(dataframe.at[i,'fitnessWT1'])
        dataframe.at[i, 'error1'] =  math.log(np.longdouble(2*(z*(statistics.stdev(dataframe.at[i, 'readsWT1']))/(np.sqrt(len(dataframe.at[i,'readsWT1'])))))*N/n1*a, (2))/gfHO1
    
    if (dataframe.at[i,'readsWT2']) == [] or len(dataframe.at[i, 'readsWT2']) < 6:
        dataframe.at[i,'reads_avr2'] =  ''
        dataframe.at[i,'fitness_avr2'] =  ''
        dataframe.at[i, 'error2'] = 0

    else:
        dataframe.at[i,'reads_avr2'] = statistics.median(dataframe.at[i,'readsWT2'])
        dataframe.at[i,'fitness_avr2'] = statistics.mean(dataframe.at[i,'fitnessWT2'])
        dataframe.at[i,'fitness_avr2median'] = statistics.median(dataframe.at[i,'fitnessWT2'])
        dataframe.at[i, 'error2'] =  math.log(np.longdouble(2*(z*(statistics.stdev(dataframe.at[i, 'readsWT2']))/(np.sqrt(len(dataframe.at[i,'readsWT2'])))))*N/n2*a, (2))/gfHO2
#%%
##filter on values that you want to see 
readsWT1 = []
readsWT2 = [] 
err1 = []
err2 = []
WT1 =[]
WT2 = []
WT1med = []
WT2med = []
chromosome  = 'chrXII'
for i in range(0,len(dataframe)):
    if  type(dataframe.at[i,'reads_avr1']) == float and  type(dataframe.at[i,'reads_avr2']) == float and dataframe.at[i, 'reads_avr1']<500 and len(dataframe.at[i,'readsWT1']) >5 and len(dataframe.at[i,'readsWT2']) >5 :
        readsWT1.append(dataframe.at[i,'reads_avr1'])
        readsWT2.append((dataframe.at[i,'reads_avr2']))
        err1.append((dataframe.at[i,'error1']))
        err2.append((dataframe.at[i,'error2']))
        WT1.append(dataframe.at[i,'fitness_avr1'])
        WT2.append(dataframe.at[i,'fitness_avr2'])
        WT1med.append(dataframe.at[i,'fitness_avr1median'])
        WT2med.append(dataframe.at[i,'fitness_avr2median'])
        
#and dataframe.at[i, 'chromosome'] == chromosome 
#%%        
#compute mean square error
a = np.array(WT1) 
b = np.array(WT2) 
mses = np.sqrt(((a-b)**2).mean())

aa = np.array(WT1med) 
bb = np.array(WT2med) 
msesmed = np.sqrt(((aa-bb)**2).mean())
#%%
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
print ("y=%.6fx+(%.6f)"%(z[0],z[1]))

plt.legend()
plt.savefig('afb.jpg', dpi=1200)
plt.show()

#plot fitness WT1 and WT2
plt.figure ()   
#plt.plot(WT1, WT2, 'x', color = 'black', label = 'mean, MSE = ' +str(mses))
plt.plot(WT1med, WT2med, 'x', color = 'red', label = 'median, MSE=' +str(msesmed))
plt.xlabel('fitness WT1')
plt.ylabel('fitness WT2')  
plt.title('Fitness WT1 and WT2 (median),' )#(+str(chromosome))
plt.axline([0, 0], [1, 1], label = 'x=y')

# calculate the trendline
z = np.polyfit(WT2, WT1, 1)
p = np.poly1d(z)
pylab.plot(WT1,p(WT1),"r--", label = 'fit of scatter')
plt.legend()
print ("y=%.6fx+(%.6f)"%(z[0],z[1]))
plt.show()

plt.figure()
plt.plot(readsWT1, 'o')