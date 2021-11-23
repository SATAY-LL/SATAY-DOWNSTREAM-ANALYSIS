# -*- coding: utf-8 -*-
"""
Created on Tue May 18 15:40:27 2021

@author: floor
"""
import pandas as pd
from statistics import mean
import numpy as np
import csv
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

#STEP III:
#compare our fitness values to cellmap data
#get data from step II and put it in a dataframe
dataframe = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
dataframe['fitness'] = ''
dataframe['fitness_avr'] = ''
dataframe['fitness_variance'] = ''
dataframe['fitness_diff'] = ''
dataframe['fitness_singlemutant'] = ''

dataframe['fitness']= dataframe['fitness'].astype('object')

with open('data_step2_fitnessWTB_excelled.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
dataframe['fitness'] = data

for i in range(0,len(dataframe)):
    if (dataframe.at[i,'fitness']) == []:
        dataframe.at[i,'fitness_avr'] =  'no insertions'
    else:
        dataframe.at[i,'fitness_avr'] = mean(dataframe.at[i,'fitness'])
        dataframe.at[i,'fitness_variance'] = np.var(dataframe.at[i,'fitness'] )

#data from cellmap:
SMFitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\strain_ids_and_single_mutant_fitness.xlsx")
qianfitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\fitnesscian.xls")
Qianlist = []
SMlist = [] 

# Genes from my data and cellmap with same name are put in the same row in dataframe.
for i in range(0, len(SMFitness)):
    print(i)
    for j in range(0, len(qianfitness)):
        if (SMFitness.at[i,'Single mutant fitness (26°)']) == []:
            dataframe.at[i,'fitness_diff'] =  'no insertions'
        elif (SMFitness.at[i,'Systematic gene name']) == qianfitness.at[j,'Systematic gene name'] and 'sn' in SMFitness.at[i,'Strain ID'] and type(SMFitness.at[i,'Single mutant fitness (26°)']) == np.float64 and type(qianfitness.at[j,'SC fitness']) == np.float64 and str(SMFitness.at[i,'Single mutant fitness (26°)'])!= 'nan' and str(qianfitness.at[j,'SC fitness'])!= 'nan':##Only use this for Cellmap, to focus on sn strain?
            # dataframe.at[i,'fitness_diff'] = (qianfitness.at[j,'SC fitness'])-(SMFitness.at[i,'Single mutant fitness (26°)'])
            # dataframe.at[i,'fitness Qian'] = qianfitness.at[j,'SC fitness']
           # dataframe.at[i,'fitness SM'] = SMFitness.at[i,'Single mutant fitness (26°)'
            SMlist.append(SMFitness.at[i,'Single mutant fitness (26°)'])
            Qianlist.append(qianfitness.at[j,'SC fitness'])
            break

#dataSMALL =  dataframe[ 
#      (len(dataframe['fitness'])<6) 
#      ]
#print(dataSMALL)    
#smallVAR = dataframe['fitness_variance'].mean()  


#making plots >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           

# Qianlist = []
# SMlist = [] 
# chromosome  = 'chrI'
# for i in range(0,7149):
#     if  type(dataframe.at[i,'fitness SM']) == np.float64 and type(dataframe.at[i,'fitness Qian']) == np.float64 and str(dataframe.at[i,'fitness SM'])!= 'nan' and str(dataframe.at[i,'fitness Qian'])!= 'nan':# and len(dataframe.at[i,'fitness']) >5 and str(dataframe.at[i,'fitness_singlemutant'])!= 'nan'  and type(dataframe.at[i,'fitness_avr']) == float and  ((dataframe.at[i, 'fitness_avr']) <200) : #and  (sum(dataframe.at[i, 'fitness']) <2000)
#         SMlist.append(dataframe.at[i,'fitness SM'])
#         Qianlist.append((dataframe.at[i,'fitness_avr']))

differencelist = []
for i in range(0,7149):
    if type(dataframe.at[i,'fitness SM']) == np.float64 and str(dataframe.at[i,'fitness SM'])!= 'nan' :
        differencelist.append(dataframe.at[i,'fitness SM']-dataframe.at[i,'fitness Qian']) 
variancelist = []
for i in range(0,7149):
    if type(dataframe.at[i,'fitness_variance']) == np.float64:
        variancelist.append(dataframe.at[i,'fitness_variance'])
        
newlistQ = []
newlistSM = [] 
chromosome  = 'chrI'
for i in range(0,len(SMlist)):
     #if  type(dataframe.at[i,'fitness Qian']) == np.float64 and str(dataframe.at[i,'fitness'])!= 'nan' and str(dataframe.at[i,'fitness Qian'])!= 'nan' and len(dataframe.at[i,'fitness']) >5: #and str(dataframe.at[i,'fitness_singlemutant'])!= 'nan'  and type(dataframe.at[i,'fitness_avr']) == float and  ((dataframe.at[i, 'fitness_avr']) <200) : #and  (sum(dataframe.at[i, 'fitness']) <2000)
        if SMlist[i]<0.98 or SMlist[i] > 1.02:
                 if Qianlist[i]<0.98 or Qianlist[i] > 1.02:
                     newlistQ.append(Qianlist[i])
                     newlistSM.append(SMlist[i])


plt.style.use('Solarize_Light2')

plt.plot(newlistSM, newlistQ,  'x', color='black', alpha=0.3)
#print(mean((differencelist)))
#print(np.var(differencelist))
z = np.polyfit(newlistSM, newlistQ, 1)
p = np.poly1d(z)
y_hat = np.poly1d(z)(newlistSM)

plt.plot(newlistSM, y_hat, "r--", lw=1)
# text = f"$newlistQ={z[0]:0.3f}\;newlistSM{z[1]:+0.3f}$\n$R^2 = {r2_score(newlistQ,y_hat):0.3f}$"
# plt.gca().text(0.05, 0.95, text,transform=plt.gca().transAxes,
#         fontsize=9, verticalalignment='top', horizontalalignment= 'right',color = 'blue')
plt.plot(newlistSM,p(newlistSM),label="trendline", color = 'blue')

plt.xlabel('Fitness CellMap')
plt.ylabel('Fitness Qian et al.')  
# plt.ylim(0.95, 1.06)
# plt.xlim(0.95,1.06)
plt.title('Comparison of fitness CellMap and Qian et al.  ' )#+str(chromosome))
plt.legend()

plt.savefig('afb7.jpg', dpi=1200)
plt.show()

