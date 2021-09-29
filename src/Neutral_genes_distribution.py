# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:05:51 2021

@author: floor
"""
import pandas as pd
import csv
from collections import Counter
import matplotlib.pyplot as plt

neutralgenes = []
SMFitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\strain_ids_and_single_mutant_fitness.xlsx")
for i in range(0,len(SMFitness)):
    if SMFitness.at[i,'Single mutant fitness (26°)'] >0.980 and SMFitness.at[i,'Single mutant fitness (26°)'] <1.020:
        neutralgenes.append(SMFitness.at[i,"Systematic gene name"]) #can also be Stain ID or Allele/gene name
del (i)

tabel = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
tabel['reads']= ''
tabel['reads']= tabel['reads'].astype('object')

with open('data_step2_readsWT1.csv', newline='') as f:
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = list(reader)
tabel['reads'] = data
del(data, reader)

neutral_reads = [ ]
all_reads = []
neutralgene_list = [ ]
for i in range(0,len(tabel)):
    if tabel.at[i,"gene"] in neutralgenes:
        neutral_reads.append(tabel.at[i,"reads"])
    all_reads.append(tabel.at[i,"reads"])
    
neutral_flat_list_ = [item for sublist in neutral_reads for item in sublist]
neutral_totalins = len(neutral_flat_list_)                  
neutral_totalreads = sum(neutral_flat_list_)
avrRperI_neutral = neutral_totalreads/ neutral_totalins

flat_list = [item1 for sublist1 in all_reads for item1 in sublist1]
totalins= len(flat_list)
totalreads= sum(flat_list)
avrRperI_tot = totalreads/totalins

#count abundances
neutral_abundance_reads_counter = Counter(neutral_flat_list_) #abundance neutral reads
tot_abundance_reads_counter = Counter(flat_list) #abundance total 

neutral_abundance_reads = pd.DataFrame.from_dict(neutral_abundance_reads_counter, orient='index').reset_index()
neutral_abundance_reads.columns = ["number of reads","abundance"]
abN =  neutral_abundance_reads.sort_values(by ='number of reads' )
x_data = abN['number of reads'].reset_index(drop=True).truncate(4, 400, copy = False)
y_data = abN['abundance'].reset_index(drop=True).truncate(4, 400, copy = False)

tot_abundance_reads = pd.DataFrame.from_dict(tot_abundance_reads_counter, orient='index').reset_index()
tot_abundance_reads.columns = ["number of reads","abundance"]
abT =  tot_abundance_reads.sort_values(by ='number of reads' )
x_data_TOT = abT['number of reads'].reset_index(drop=True).truncate(4, 400, copy = False)
y_data_TOT = abT['abundance'].reset_index(drop=True).truncate(4, 400, copy = False);


#plots
fig, ax = plt.subplots()

plt.style.use('Solarize_Light2')
ax.bar(x_data_TOT, y_data_TOT,  label='All genes') #bottom=y_data,
ax.bar(x_data, y_data, label='Neutral genes')


ax.set_ylabel('Reads abundance')
ax.set_xlabel("#reads")
ax.set_title('Distribution of reads for neutral genes compared to all genes')

plt.savefig('afb.jpg', dpi=1200)

plt.show()
