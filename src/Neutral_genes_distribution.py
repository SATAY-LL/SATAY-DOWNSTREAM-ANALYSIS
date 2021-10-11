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

all_reads = [ ]
neutralgene_list = [ ]
for i in range(0,len(tabel)):
    if tabel.at[i,"gene"] in neutralgenes:
        all_reads.append(tabel.at[i,"reads"])
        
flat_list = [item for sublist in all_reads for item in sublist]
totalins = len(flat_list)                  
totalreads = sum(flat_list)

#count abundances
abundance_reads_counter = Counter(flat_list)

abundance_reads = pd.DataFrame.from_dict(abundance_reads_counter, orient='index').reset_index()
abundance_reads.columns = ["number of reads","abundance"]
ab =  abundance_reads.sort_values(by ='number of reads' )
x_data = ab['number of reads'].reset_index(drop=True).truncate(4, 400, copy = False)
y_data = ab['abundance'].reset_index(drop=True).truncate(4, 400, copy = False);
plt.plot(x_data,y_data)

