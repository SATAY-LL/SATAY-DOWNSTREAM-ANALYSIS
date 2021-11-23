# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 15:08:17 2021

@author: floor
"""
import wiggelen
import pandas as pd
import random
import collections

## This file is used to randomly delete a portion of the reads
#input: wigfile 
#output: dataframe with new and old number of reads per transposon (location and chr)

#this function is used in script_STEP1 to calculate fitness

def portion_of_allreads(filepath_and_name, percentage): #for 50% percentage= 0.5
    list = []
    for x in wiggelen.walk(open(filepath_and_name)):
        y = [x[0], x[1], x[2]]
        list.append(y)
    data = pd.DataFrame(list,columns=['chromosome','tn start position', '#reads'])
    data['transposon number'] = ''
    
    for i in range(0,len(data)):
        data.at[i,'transposon number']= i
    
    print('Start making read list... (takes a while)')
    transposon = []
    for tn in range (0, len(data)):
            for ftns in range (0, data.at[tn, '#reads']):
                    transposon.append(data.at[tn, 'transposon number'])
    print('...end making read list')
    k = int(len(transposon)*percentage)
    random.seed(2)
    transposon_50percent = random.sample(transposon, k)
    
    RpT50percent = collections.Counter(transposon_50percent)
    df50percent = pd.DataFrame.from_dict(RpT50percent, orient='index')
    merged = data.merge(df50percent, left_index=True, right_index=True, how='inner')
    merged.reset_index(drop=True, inplace = True)
    merged.columns =['chromosome','tn start position', '#reads', 'transposon number', '50% reads']
    print('End of def_portion_of_allreads')
    return(merged)

    
    
    