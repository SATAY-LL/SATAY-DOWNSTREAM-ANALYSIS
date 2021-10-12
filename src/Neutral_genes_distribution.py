# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 14:05:51 2021

@author: floor
"""
import pandas as pd
import csv
from collections import Counter
import matplotlib.pyplot as plt
from numpy import trapz
import numpy as np

#readsforeverygenefile = 'data_step2_Dnrp1_1BREADS_exelled.csv' 
def distribution_neutralgenes(readsforeverygenefile):
    neutralgenes = []
    positivegenes = []
    SMFitness = pd.read_excel(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\strain_ids_and_single_mutant_fitness.xlsx")
    for i in range(0,len(SMFitness)):
        if SMFitness.at[i,'Single mutant fitness (26°)'] <1.005 and SMFitness.at[i,'Single mutant fitness (26°)'] >0.995 : #nd SMFitness.at[i,'Single mutant fitness (26°)'] <1.010
            neutralgenes.append(SMFitness.at[i,"Systematic gene name"]) #can also be Stain ID or Allele/gene name
        elif SMFitness.at[i,'Single mutant fitness (26°)'] >1.05 : #nd SMFitness.at[i,'Single mutant fitness (26°)'] <1.010
            positivegenes.append(SMFitness.at[i,"Systematic gene name"]) #can also be Stain ID or Allele/gene name
    del (i)
    
    tabel = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
    tabel['reads']= ''
    tabel['reads']= tabel['reads'].astype('object')
    
    with open(readsforeverygenefile, newline='') as f:
        reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
        data = list(reader)
        tabel['reads'] = data
    del(data, reader)
    
    neutral_reads = [ ]
    positive_reads = []
    mean_per_gene_neutral=[]
    mean_per_gene_pos = []
    all_reads = []
    mean_all_reads = []
    minimal_insertions = 0
    for i in range(0,len(tabel)):
        if tabel.at[i,"gene"] in neutralgenes and len(tabel.at[i,'reads'])>minimal_insertions:
            neutral_reads.append(tabel.at[i,"reads"])
            mean_per_gene_neutral.append(int(np.mean(tabel.at[i,"reads"])))
        elif tabel.at[i,"gene"] in positivegenes and len(tabel.at[i,'reads'])>minimal_insertions:
            positive_reads.append(tabel.at[i,"reads"])
            mean_per_gene_pos.append(np.mean(tabel.at[i,"reads"]))
        if  len(tabel.at[i,'reads'])>0:
            all_reads.append(tabel.at[i,"reads"])
            mean_all_reads.append(int(np.mean(tabel.at[i,"reads"])))
    avr_over_neutral_genes= np.mean(mean_per_gene_neutral)
    avr_over_pos_genes = np.mean(mean_per_gene_pos)
    avr_over_all_genes = np.mean(mean_all_reads)
    
    neutral_flat_list_ = [item for sublist in neutral_reads for item in sublist]
    neutral_totalins = len(neutral_flat_list_)                  
    neutral_totalreads = sum(neutral_flat_list_)
    avrRperI_neutral = neutral_totalreads/ neutral_totalins
    
    flat_list = [item1 for sublist1 in all_reads for item1 in sublist1]
    totalins= len(flat_list)
    totalreads= sum(flat_list)
    avrRperI_tot = totalreads/totalins
    
    #count abundances
    neutral_abundance_reads_counter = Counter((avr_over_neutral_genes)) #abundance neutral reads
    tot_abundance_reads_counter = Counter((avr_over_all_genes)) #abundance total 
    
    neutral_abundance_reads = pd.DataFrame.from_dict(neutral_abundance_reads_counter, orient='index').reset_index()
    neutral_abundance_reads.columns = ["number of reads","abundance"]
    abN =  neutral_abundance_reads.sort_values(by ='number of reads' )
    
    x_data = abN['number of reads'].reset_index(drop=True)#.truncate(0, 94, copy = False)
    y_data = abN['abundance'].reset_index(drop=True)#.truncate(0,94, copy = False)
    areaN= trapz(y_data, dx=1)
    tot_abundance_reads = pd.DataFrame.from_dict(tot_abundance_reads_counter, orient='index').reset_index()
    tot_abundance_reads.columns = ["number of reads","abundance"]
    abT =  tot_abundance_reads.sort_values(by ='number of reads' )
    x_data_TOT = abT['number of reads'].reset_index(drop=True)#.truncate(0, 187, copy = False)
    y_data_TOT = abT['abundance'].reset_index(drop=True)#.truncate(0, 187, copy = False);
    areaT = trapz(y_data_TOT, dx=1)
    y_data = y_data/areaN
    y_data_TOT = y_data_TOT/areaT
    return(x_data, y_data, x_data_TOT, y_data_TOT, avr_over_neutral_genes, avr_over_pos_genes, avrRperI_tot, totalreads)

x_data2, y_data2, x_data_TOT2, y_data_TOT2, avr_over_neutral_genes2, avr_over_pos_genes2, TA, totalA = distribution_neutralgenes('data_step2_Dnrp1_1aREADSexcelled.csv')
x_data1, y_data1, x_data_TOT1, y_data_TOT1, avr_over_neutral_genes1, avr_over_pos_genes1, TB, totalB = distribution_neutralgenes('data_step2_Dnrp1_1BREADS_exelled.csv')
x_dataM, y_dataM, x_data_TOTM, y_data_TOTM, avr_over_neutral_genesM, avr_over_pos_genesM, TM, totalM = distribution_neutralgenes('data_step2_Dnrp1_1abmergedREADSexcelled.csv')

#plots
fig, ax = plt.subplots()

plt.style.use('Solarize_Light2')
ax.bar(x_data_TOT1, y_data_TOT1,  label='All genes') #bottom=y_data,
ax.bar(x_data1, y_data1, label='Neutral genes', color = 'red', alpha = 0.5)


ax.set_ylabel('Normalised reads abundance')
ax.set_xlabel("Mean #reads/tn over genes")
ax.set_title('Distribution of reads, WT B')
ax.yaxis.get_data_interval()
ax.legend()
plt.xlim([0,100])
plt.ylim([0.0,0.17])
plt.savefig('afb3.jpg', dpi=1200)

plt.show()

fig, ax = plt.subplots()

plt.style.use('Solarize_Light2')
ax.bar(x_data_TOT2, y_data_TOT2,  label='All genes') #bottom=y_data,
ax.bar(x_data2, y_data2, label='Neutral genes', color = 'red', alpha = 0.5)
print(totalA, totalB, totalM)

ax.set_ylabel('Normalised reads abundance')
ax.set_xlabel("Mean #reads/tn over genes")
ax.set_title('Distribution of reads, WT A')
ax.yaxis.get_data_interval()
ax.legend()
plt.ylim([0.0,0.17])
plt.xlim([0,100])
plt.savefig('afb1.jpg', dpi=1200)

plt.show()

fig, ax = plt.subplots()

plt.style.use('Solarize_Light2')
ax.bar(x_data_TOTM, y_data_TOTM,  label='All genes') #bottom=y_data,
ax.bar(x_dataM, y_dataM, label='Neutral genes', color = 'red', alpha = 0.5)


ax.set_ylabel('Normalised reads abundance')
ax.set_xlabel(" Mean #reads/tn over genes")
ax.set_title('Distribution of reads, WT merged')
ax.yaxis.get_data_interval()
ax.legend()
plt.xlim([0,100])
plt.ylim([0.0,0.17])
plt.savefig('afb2.jpg', dpi=1200)

plt.show()