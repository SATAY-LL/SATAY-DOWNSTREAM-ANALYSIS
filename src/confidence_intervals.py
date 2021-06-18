# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 17:09:36 2021

@author: floor
"""


# in this function, the method of agresti and coull is used to calculate the binominal proportion interval for each transposon
#I used the two files from leila 'WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene_insertions.txt' and 'WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene_insertions.txt'
#for input, you should fill in the  * r'C:filepath and file name' * of the two files of Leila
#the output is for each gene, all the insertion and their confidence intervals
def confidenceintervals (filedirectorya, filedirectoryb): 

    import pandas as pd
    import numpy as np
    import statsmodels.stats.proportion as sm

#import data from \WT_merged-DpnII-NlaIII-a and b
    dfA = pd.read_csv(filedirectorya, sep = "\t")                           
    dfB = pd.read_csv(filedirectoryb, sep = "\t")                           

    #find total number of reads and number of transposons of a
    ins = []
    sumreadsA = []
    for a in range (0,len(dfA)):
        if dfA.at[a, 'Reads per insertion location'] == '[]':
            nA = 0
        else:
            nA = (list(((dfA.at[a, 'Reads per insertion location'])[1:-1]).split(', ')))
            sumreadsA.extend((np.array(nA).astype(float)))
            ins.append(len(nA))
    nA = sum(sumreadsA) #total reads a
    insertion  = sum(ins)
    averagereadsA = nA/insertion
    del(a, insertion, ins)
    
    #find total number of reads and number of transposons of b
    sumreadsB = []
    ins = []
    for b in range (0,len(dfB)):
        if dfB.at[b, 'Reads per insertion location'] == '[]':
            nB = 0
        else:
            nB = (list(((dfB.at[b, 'Reads per insertion location'])[1:-1]).split(', ')))
            sumreadsB.extend((np.array(nB).astype(float)))
            ins.append(len(nB))
    nB = sum(sumreadsB) #total reads b
    insertion  = sum(ins)
    averagereadsB = nB/insertion
    del(b, insertion, ins)
    
    
    
    dfA['confidence intervals'] = ''
    dfB['confidence intervals'] = ''
    a = []
    Confidence_intervalsA = []
    for i in range (0,len(dfA)):
        if dfA.at[i, 'Reads per insertion location'] == '[]':
            reads = 0
            Confidence_intervalsA.append([0,0])
        else:
            reads = (list(((dfA.at[i, 'Reads per insertion location'])[1:-1]).split(', ')))
            a = []
            for j in range(0,len(reads)):
                confidence_interval = (sm.proportion_confint(int(reads[j]), nA, alpha=0.05, method='agresti_coull')) 
                a.append(list([nA*x for x in confidence_interval]))
            Confidence_intervalsA.append(a)
    dfA['confidence intervals'] = Confidence_intervalsA
    
    del(i, reads, a, j)
    b = []
    Confidence_intervalsB = []
    for i in range (0,len(dfB)):
        if dfB.at[i, 'Reads per insertion location'] == '[]':
            reads = 0
            Confidence_intervalsB.append([0,0])
        else:
            reads = (list(((dfB.at[i, 'Reads per insertion location'])[1:-1]).split(', ')))
            b = []
            for j in range(0,len(reads)):
                confidence_interval = (sm.proportion_confint(int(reads[j]), nB, alpha=0.05, method='agresti_coull')) 
                b.append(list([nB*x for x in confidence_interval]))
            Confidence_intervalsB.append(b)
    dfB['confidence intervals'] = Confidence_intervalsB
    del(b, i, j, reads)
    return(dfA, dfB)

test = confidenceintervals(r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\data\WT_merged-DpnII-NlaIII-a_trimmed.sorted.bam_pergene_insertions.txt', r'C:\Users\floor\OneDrive\Documenten\MASTER\MEP\codes\LaanLab-SATAY-DataAnalysis\Python_scripts\data\WT_merged-DpnII-NlaIII-b_trimmed.sorted.bam_pergene_insertions.txt')
        
        
        
        
        
        
        
        