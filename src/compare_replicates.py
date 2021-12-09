# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:17:53 2021

@author: floor
"""
import pandas as pd
import numpy as np
import csv
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import pylab
import statistics
from scipy.stats import norm
import math
from statsmodels import robust
from sklearn.metrics import r2_score
from collections import Counter
import scipy
from scipy.optimize import curve_fit
import scipy.stats as ss


#STEP III:
#compare our fitness of 1 file to another
#get data from step II and put it in a tabel
tabel = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\data\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
tabel['fitness1merged'] = ''
tabel['fitness2merged'] = ''
tabel['fitness_avr1merged'] =  ''
tabel['fitness_avr2merged'] =  ''
tabel['fitness1a'] = ''
tabel['fitness2a'] = ''
tabel['fitness_avr1a'] =  ''
tabel['fitness_avr2a'] =  ''


#import fitness TR 1
tabel['fitness1a']= tabel['fitness1a'].astype('object')

with open(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\processed\adding random reads\data_step2_fitnessdNrp1-2-merged_10percent_EXC.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data3 = list(reader)
tabel['fitness1a'] = data3

#import fitness tech rep 2
tabel['fitness2a']= tabel['fitness2a'].astype('object')

with open(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\processed\adding random reads\data_step2_fitnessdNrp1-1-merged_10percent_EXC.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data4 = list(reader)
tabel['fitness2a'] = data4

del(data3,data4)


#find median of fitness values
for i in range(0, len(tabel)-1):
    if tabel.at[i, 'fitness1a'] != []:
        tabel.at[i, 'fitness_avr1a'] =  statistics.mean(tabel.at[i, 'fitness1a'])


for i in range(0, len(tabel)-1):
    if tabel.at[i, 'fitness2a'] != []:
        tabel.at[i,'fitness_avr2a'] =  statistics.mean(tabel.at[i, 'fitness2a'])

fitness1a =[]
fitness2a = []
difference = []
difference_notrounded = []
genes10 = []
chromosome  = 'chrXII'
for i in range(0,len(tabel)):
    if  type(tabel.at[i,'fitness_avr1a']) == float and  type(tabel.at[i,'fitness_avr2a']) == float  and len(tabel.at[i,'fitness2a']) >5 and len(tabel.at[i,'fitness1a']) >5: #and tabel.at[i, 'chromosome'] == chromosome:
       
        fitness1a.append(tabel.at[i,'fitness_avr1a'])
        fitness2a.append(tabel.at[i,'fitness_avr2a'])
        genes10.append(tabel.at[i, 'gene'])

#compute difference between TR1 and TR2
for j in range(0,len(fitness1a)-1):
    xsingle= (fitness1a[j]-fitness2a[j])
    difference.append(round( xsingle, 2 ) )
    difference_notrounded.append(xsingle)
abundance =Counter(difference)


#import fitness merged 1
tabel['fitness1merged']= tabel['fitness1merged'].astype('object')

with open(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\processed\adding random reads\data_step2_fitnessdNrp1-1-merged_60percent_EXC.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data1 = list(reader)
tabel['fitness1merged'] = data1

#import fitness merged 2
tabel['fitness2merged']= tabel['fitness2merged'].astype('object')

with open(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\processed\adding random reads\data_step2_fitnessdNrp1-2-merged_60percent_EXC.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data2 = list(reader)
tabel['fitness2merged'] = data2

del(data1,data2)


transposon = []
#median of fitness values per gene
for i in range(0, len(tabel)-1):
    if tabel.at[i, 'fitness1merged'] != []:
        tabel.at[i, 'fitness_avr1merged'] =  statistics.mean(tabel.at[i, 'fitness1merged'])       
        
for i in range(0, len(tabel)-1):
    if tabel.at[i, 'fitness2merged'] != []:
        tabel.at[i,'fitness_avr2merged'] =  statistics.mean(tabel.at[i, 'fitness2merged'])

fitness1merged =[]
fitness2merged = []
differencemerged = []
differencemerged_notrounded = []
genes60= []
chromosome  = 'chrXII'
## calculate diffenences:
for i in range(0,len(tabel)):
    if  type(tabel.at[i,'fitness_avr1merged']) == float and tabel.at[i, 'gene'] in genes10 and  type(tabel.at[i,'fitness_avr2merged']) == float and len(tabel.at[i,'fitness2merged']) >5 and len(tabel.at[i,'fitness1merged']) >5  : #and tabel.at[i, 'chromosome'] == chromosome:
       
        fitness1merged.append(tabel.at[i,'fitness_avr1merged'])
        fitness2merged.append(tabel.at[i,'fitness_avr2merged'])
        genes60.append(tabel.at[i, 'gene'])
for j in range(0,len(fitness1merged)-1):
    xmerged= (fitness1merged[j]-fitness2merged[j])
    
    
    differencemerged.append( round(xmerged,2))
    differencemerged_notrounded.append(xmerged)
abundancemerged =Counter(differencemerged)

varianceBR = statistics.pvariance(differencemerged_notrounded)
varianceTR = statistics.pvariance(difference_notrounded)


### PLOT OF HISTOGRAM DIFFERENCE + GAUSSIAN::
plt.style.use('Solarize_Light2')    
binwidth  = 0.01
_,bins,_ = plt.hist(differencemerged_notrounded, bins=np.arange(min(differencemerged_notrounded), max(differencemerged_notrounded) + binwidth, binwidth), label='Difference between replicates')

muM, stdM = norm.fit(differencemerged_notrounded) 
p = np.multiply(norm.pdf(bins, muM, stdM), 305)
#plt.plot(bins, p,  linewidth=2, label = 'Normal distr. tech. rep.', color = 'indigo') #, label = 'normal distribution of single, STD=' +  str(round(stdS,3)))
plt.xlabel('Difference between datasets')
plt.title('Histogram 99% samples of -Nrp1')
plt.ylabel('Abundance')
plt.legend(bbox_to_anchor=(0.55, 0.85))
plt.savefig('histogram differences merged with normal distribution', bbox_inches='tight', dpi=1000)
## SCATTER PLOTS OF 2REPLICATES::
#plot reads wt1 wt2     
plt.figure()
plt.style.use('Solarize_Light2')
plt.scatter(fitness1merged, fitness2merged,  c='black', alpha = 0.1, marker="x")


plt.axline([0, 0], [1, 1], label = 'x=y')

# calculate the trendline
z = np.polyfit(fitness1merged, fitness2merged, 1)
p = np.poly1d(z)
y_hat = np.poly1d(z)(fitness1merged)

plt.plot(fitness1a, y_hat, "r--", lw=1)
text = f"$fitness2a={z[0]:0.3f}\;fitness1a{z[1]:+0.3f}$\n$R^2 = {r2_score(fitness2merged,y_hat):0.3f}$"
plt.gca().text(0.05, 0.95, text,transform=plt.gca().transAxes,
        fontsize=9, verticalalignment='top', horizontalalignment= 'right',color = 'blue')
plt.plot(fitness1merged,p(fitness1merged),label="trendline", color = 'red')

plt.xlabel('Fitness Nrp1 1 merged 60%')
plt.ylabel('Fitness Nrp1 2 merged 60%')  
plt.title('Fitness Nrp1 1 merged 60% vs. 2 merged 60% ' )#(+str(chromosome))
plt.legend()
plt.savefig('afb.jpg', dpi=1200)
plt.show()     



## PLOTS of DIFFERENCES
muS, stdS = norm.fit(difference)

lists = sorted(abundance.items()) # sorted by key, return a list of tuples

xdata,ydata = zip(*lists) # unpack a list of pairs into two tuples


plt.figure()
plt.title('Deviation between random samples of Nrp1-2 merged \n vs Nrp1-1 merged 10%')
plt.xlabel('Absolute difference between datasets')
plt.ylabel("Abundance")
plt.scatter(xdata, ydata, label= 'diff. 10% sample and \n 10% sample, STD='+ str(round(stdS, 3)),alpha=0.8, color='orange')
spline = scipy.interpolate.UnivariateSpline(xdata, ydata)
spline.set_smoothing_factor(11500)
#plt.plot(xdata, spline(xdata), 'green', lw=3)
halvemax= (spline(0)/2)
FWHM_points = scipy.interpolate.UnivariateSpline(xdata, ydata - halvemax, s=11500).roots()
x_valuesHWHM = [FWHM_points[1], FWHM_points[0]]
y_valueHWHM = [halvemax, halvemax]
FWHM = FWHM_points[1] -FWHM_points[0]
#plt.plot(x_valuesHWHM, y_valueHWHM, label = 'FWHM (=' +str(round(FWHM,3)) + ')', color = 'black')
plt.ylim(0, 200)
plt.xlim(-1, 1)
plt.legend()
plt.savefig('60%', bbox_inches='tight', dpi=1000)

#fit exponential functions::
# #fit to two exponential functions
# x_1halve = np.asarray(xdata[3:107])
# x_2halve = np.asarray(xdata[107:-3])
# y_1halve = np.asarray(ydata[3:107])
# y_2halve = np.asarray(ydata[107:-3])

# #exponential fit to 1st halve
# log_x_dataS_1halve = np.log(x_1halve)
# log_y_dataS_1halve = np.log(y_1halve)

# power1_1HS, power2_1HS = np.polyfit(x_1halve, log_y_dataS_1halve, 1)
# yS = np.exp(power1_1HS) * np.exp(power2_1HS*x_1halve)


# # plt.plot(x_1halve, yS, color ='indigo')

# #exponential fit to 2st halve
# log_x_dataS_2halve = np.log(x_2halve)
# log_y_dataS_2halve = np.log(y_2halve)

# power1_2HS, power2_2HS = np.polyfit(x_2halve, log_y_dataS_2halve, 1)
# y_2HS = np.exp(-power1_2HS) * np.exp(-(power2_2HS*x_2halve))


# #Find half WIDHT HALVE MAXIMUM of exponential fits
# peak_top = 0
# for i in ydata:
#     if i > peak_top:
#         peak_top = i
# half_max = peak_top/2
# left = (np.log(half_max)-power1_1HS)/power2_1HS
# right = (np.log(half_max)+power1_2HS)/(-power2_2HS)
# FWHM = right- left
# x_values = [left, right]
# y_value = [half_max, half_max]
# # plt.plot(x_values, y_value, label = 'FWHM (=' +str(round(FWHM,3)) + ')', color = 'black')
# #plt.plot(x_2halve, y_2HS, color = 'indigo' , label = 'fit 1a & 2a to two exponentials \n outer 3 values removed ')
# del (x_values, y_value, left, right, FWHM, half_max, i)


# #MERGED / second comparison plot
lists2 = sorted(abundancemerged.items()) # sorted by key, return a list of tuples


a, b = zip(*lists2) # unpack a list of pairs into two tuples
plt.figure()
plt.scatter(a, b, c='black', marker='x', label = 'diff. 60% sample and \n 60% sample, STD='+ str(round(stdM, 3)), alpha=0.8)
plt.title('Deviation between random samples of Nrp1-2 merged \n vs. Nrp1-1 merged 60% ')
plt.xlabel('Absolute difference between datasets')
plt.ylabel("Abundance")
plt.legend(loc='center left', bbox_to_anchor=(0.55, 0.8))

#plt.ylim(0, 800)
plt.xlim(-1, 1)
spline = scipy.interpolate.UnivariateSpline(a, b)
spline.set_smoothing_factor(15000)
#plt.plot(a, spline(a), 'red', lw=3, alpha=0.6)
halvemax= (spline(0)/2)
FWHM_merg = scipy.interpolate.UnivariateSpline(a, b - halvemax, s=15000).roots()
x_valuesHWHM = [FWHM_merg[1], FWHM_merg[0]]
y_valueHWHM = [halvemax, halvemax]
FWHM_M = FWHM_merg[1]-FWHM_merg[0]
#plt.plot(x_valuesHWHM, y_valueHWHM, label = 'FWHM (=' +str(round(FWHM_M,3)) + ')', color = 'black')
plt.ylim(0, 200)
plt.xlim(-1, 1)
plt.legend()
plt.savefig('50%', bbox_inches='tight', dpi=1000)

# #also fit to two exponential functions to the second comparison
# a_1halve = np.asarray(a[:106])
# a_2halve = np.asarray(a[104:-7])
# b_1halve = np.asarray(b[:106])
# b_2halve = np.asarray(b[104:-7])

# #exponential fit to 1st halve
# log_x_data_1halve = np.log(a_1halve)
# log_y_data_1halve = np.log(b_1halve)

# power1_1H, power2_1H = np.polyfit(a_1halve, log_y_data_1halve, 1)
# y = np.exp(power1_1H) * np.exp(power2_1H*a_1halve)

# plt.plot(a_1halve, y, color ='red', label = '')

# #exponential fit to 2st halve
# log_x_data_2halve = np.log(a_2halve)
# log_y_data_2halve = np.log(b_2halve)

# power1_2H, power2_2H = np.polyfit(a_2halve, log_y_data_2halve, 1)
# y_2H = np.exp(-power1_2H) * np.exp(-(power2_2H*a_2halve))

# plt.plot(a_2halve, y_2H, color = 'red', label = 'fit diff. merged to two exponentials \n outer 3 values removed ')
# #Find FULL WIDHT HALVE MAXIMUM
# peak_top = 0
# for i in b:
#     if i > peak_top:
#         peak_top = i
# half_max = peak_top/2
# left = (np.log(half_max)-power1_1H)/power2_1H
# right = (np.log(half_max)+power1_2H)/(-power2_2H)
# FWHM = right- left
# x_values = [left, right]
# y_value = [half_max, half_max]
# plt.plot(x_values, y_value, label = 'FWHM (=' +str(round(FWHM,3)) + ')', color = 'black')

#del (x_values, y_value, left, right, FWHM, half_max, i)
#plt.legend(loc='center left', bbox_to_anchor=(0.55, 0.85))
# # # Plot the PDF of a gaussian.
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# p = np.multiply(norm.pdf(sorted(difference), muS, stdS), 80)
  
##plot normal distribution::
# plt.plot(sorted(difference), p,  linewidth=2, label = 'normal distr. tech. rep.', color = 'indigo') #, label = 'normal distribution of single, STD=' +  str(round(stdS,3)))

# # Plot the PDFof merged.
# xminM, xmaxM = plt.xlim()
# xM = np.linspace(xminM, xmaxM, 100)
# pM = np.multiply(norm.pdf(sorted(differencemerged_notrounded), muM, stdM), 80)
  
# #plt.plot(sorted(differencemerged_notrounded), pM, color = 'black', linewidth=2, label = 'normal distr. merged') #, label = 'normal distribution of merged, STD=' +  str(round(stdM,3)))

plt.savefig('differences fit to exponential', bbox_inches='tight', dpi=1000)

