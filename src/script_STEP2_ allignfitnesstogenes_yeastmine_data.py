import pandas as pd
from script_STEP1_fitness_calculation_from_reads import fitness1

# STEP II
# The fitness values that we calculated have an insertion location (bp + chr). 
# Using data from cellmap, insertions + fitness values are sorted for each gene.
# For this, we use the function 'fitness', defined in step 1

# Retreive data
genes_start_stop = pd.read_csv(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\ChromosomeRegion_AllGenes.tsv", sep = "\t", names = ["chromosome","gene", "start bp", "stop bp",  "+1 or -1 strand"] )
fitness_insertions = fitness1(r"C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\SATAY-DOWNSTREAM-ANALYSIS\datasets\data nrp1\dnrp1-2_merged-techrep-a_techrep-b_trimmed.sorted.bam.wig")
print('Start STEP2')
#add column to df with all genes, where we are going to put the fitness values
genes_start_stop['fitness'] = ''
genes_start_stop['fitness']= genes_start_stop['fitness'].astype('object')

for i in range (0, len(genes_start_stop)):
    genes_start_stop.at[i,'fitness'] = []

# To make it 17x faster, only look at the right chromosome of yeastmine data.
# if tn startposition is within start and stop bp of gene, fitness value is added to the fitness collumn.
genes_start_stop = genes_start_stop.sort_values('chromosome')
genes_start_stop = genes_start_stop.reset_index()
for i in range(0,len(fitness_insertions)):
    print(i)
    if fitness_insertions.at[i,'chromosome'] == 'chrI':
        for j in range (0, 126):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrII':
        for j in range (126, 607):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrIII':
        for j in range (607, 810):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrIV':
        for j in range (810, 1698):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrIX':
        for j in range (1698, 1959):
           if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrV':
        for j in range (1959, 2317):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrVI':
        for j in range (2317, 2472):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrVII':
        for j in range (2472, 3112):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrVIII':
        for j in range (3112, 3455):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrX':
        for j in range (3455, 3890):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXI':
        for j in range (3890, 4260):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXII':
        for j in range (4260, 4900):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXIII':
        for j in range (4900, 5451):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXIV':
        for j in range (5451, 5910):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXV':
        for j in range (5910, 6549):
            if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrXVI':
        for j in range (6549, 7094):
           if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
    if fitness_insertions.at[i,'chromosome'] =='chrmt':
        for j in range (7094, 7149):
           if fitness_insertions.at[i, 'tn start position'] <= (genes_start_stop.at[j, 'stop bp']) and fitness_insertions.at[i, 'tn start position'] >= (genes_start_stop.at[j, 'start bp']):
                a = fitness_insertions.at[i, 'fitness']
                genes_start_stop.at[j, 'fitness'].append(a)
                break
copy = genes_start_stop          
genes_start_stop = genes_start_stop.sort_values('index')
genes_start_stop.reset_index()
genes_start_stop.to_csv('data_step2_fitnessdNrp1-2-merged_40percentB.csv')
 #sometimes csv file does seperate by comma.
#Therefor: open file in excel and only keep collumn with fitness values, and seperate by comma.
#safe as csv
print('...DONE!!')
import winsound
duration = 1000  # milliseconds
freq = 440  # Hz
winsound.Beep(freq, duration)