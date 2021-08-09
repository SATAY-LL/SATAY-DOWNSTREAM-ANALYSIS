# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:17:00 2021

@author: linigodelacruz
"""

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from collections import defaultdict
import os
import seaborn as sns
import scipy 

### Analysis of transposition sites

# def frequency_transposons(data,names_libraries):
#     freq=[]
#     for i in names_libraries.keys():
#         freq.append(data.loc[i]['Nbasepairs'].sum()/data.loc[i]['Ninsertions'].sum())
        
#     return freq

def frequency_transposons(data):
    """Compute how many bp are per transposons in average"""
    
   bp_per_tr=data.loc[:,'Nbasepairs'].sum()/data.loc[:,'Ninsertions'].sum()
        
   return bp_per_tr



def reads_per_transposon(data):
    """Median reads per insertions"""
    
  
    readspertr=data.loc[:,'Nreadsperinsrt'].median()
        
    return readspertr
    
def transposon_density(data):
    tr_per_bp=data.loc[:,'Ninsertions'].sum()/data.loc[:,'Nbasepairs'].sum()
    
    print("There are, in average", 100*tr_per_bp, "transposons", "per each 100 basepairs")
    
    return tr_per_bp

    


    
def median_feature_essentials(data,feature):
    """Compute the median of the desired column for essential genes """
    
   median_feature=data[data.loc[:,'Essentiality']==1][feature].median()  
      
   return median_feature
    
def median_feature_nonessentials(data,feature):
    """Compute the median of the desired column for non essential genes """
    
        
    median_feature=data[data.loc[:,'Essentiality']==0][feature].median()  
    
    return median_feature

# Zoom in per chromosome
def local_variation(chrom,windows,data,column):
    """Computing the mean and std of any feature
    per chromosome with a given rolling window size"""
    
    from scipy.stats import sem 
    
    df=data[data.loc[:,'Chromosome']==chrom]
    
    mean_values = df[column].rolling(windows).mean()[::windows] # series
            
    std_values = df[column].rolling(windows).std()[::windows]/np.sqrt(len(df)) # series 
   
    
    return mean_values,std_values
 

def filter_low_and_biased_reads_genes(target_data,reads_upper_th=25,tr_upper_th=2,tr_density_upper_th=0.01,tr_density_low_th=0.2):
    """Filter dataset to have enough data to work with for further measures
    like fitness per gene"""
    
    #- Discard the genes that has less than 25 reads 
    #- Discard the genes that has more 20% coverage of transposons (centromeres , bias genes), tr-density>0.2 (less than 3% of the data)
    #- Discard the genes that has less than 1% of coverage of transposons , tr-density<0.01 (around 7% of the genes that has less than 20% coverage)
    
    tr_discard_singletr=target_data[target_data['Ninsertions']<tr_upper_th]
    tr_discard_highcoverage=target_data[target_data['tr-density']>tr_density_low_th]
    tr_discard_lowcoverage=target_data[target_data['tr-density']<tr_density_upper_th]
    reads_exceptions=target_data[target_data['Nreads']<reads_upper_th]

    ### adding all the discarded genes 
    d=[tr_discard_singletr.index,tr_discard_highcoverage.index,tr_discard_lowcoverage.index, reads_exceptions.index]
    #exceptions=list(reduce(set.intersection, [set(item) for item in d ])) ## getting the intersection of all the indexes from the discarded genes
    
    d_new=[]
    for i in np.arange(0,len(d)):
        for j in np.arange(0,len(d[i])):
            d_new.append(d[i][j])
        
    bad_unique_indexes=np.unique(d_new)
    bad_df= target_data.index.isin(bad_unique_indexes)
    new_df=target_data[~bad_df] ## dataframe with the right indexes 
    
    
    for i in new_df.index:
        target_data.loc[i,'reads-per-tr']=target_data.loc[i,'Nreads']/(target_data.loc[i,'Ninsertions']-1)#there is one more transposon from the maximum reads

    return target_data