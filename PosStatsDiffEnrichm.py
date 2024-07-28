#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:50:34 2024

@author: dabrahamsson
"""

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import linregress
import seaborn as sns
import matplotlib as mlp
import matplotlib.pyplot as plt

df = pd.read_csv('BC_Average_AllNeg_PMM_withP+r.csv')
ga = pd.read_csv('Index+Labels.csv')

df.columns.values

loc1 = df.loc[:, 'Alignment ID':'MS/MS spectrum'] #everything between the chem_id and the first samples, -blanks
loc2 = df.loc[:, 'chem_id'] #combination of Average Mz and Average Rt(min) - reinserted after the blanks
loc3 = df.loc[:, 'Cornell_Pos_S-01_LMT_PS_PMM':'Cornell_Pos_S-50_LMH_LCT_PS_PMM'] #samples - does not include blanks

def stat(samples):
    df = pd.concat([loc1, loc2, samples], axis=1)
    
    print(df.shape)
    
    ids = samples.columns.values.tolist()
    mt = pd.melt(df, id_vars=['chem_id', 'Average Mz', 'Average Rt(min)'], value_vars=ids, value_name='Abundance')
    mt.columns = mt.columns.str.replace('variable', 'Sample ID') 
    
    mt = pd.merge(mt, ga, on='Sample ID', how='left') #Index spreadsheet
    mt['RPL'] = np.where(mt['MotMorph'] == 'PMM', 0, 1) #R or C in index
    mt['logA'] = np.log10(mt['Abundance']) #from line 33
    
    def get_element(my_list, position):
        return my_list[position]
    
    mt['RPL'] = mt['RPL'].astype(str) #RPL not in either sheet
    mt['Sample ID'] = mt['Sample ID'].astype(str)  #from index spreadsheet
    mt['identifier'] = mt['RPL'] + '_' + mt['Sample ID']
    mt = mt.drop(['RPL','Sample ID'], axis=1)
    mt = mt.pivot_table('logA','chem_id','identifier')
    
    axisvalues = mt.columns.values
    axisvalues_mt = pd.DataFrame({'identifier': axisvalues}) #indexed from Sample ID
    axisvalues_ = axisvalues_mt['identifier'].str.split('_').apply(lambda x: x[0])
    axisvalues_ = axisvalues_.astype(float)
    
    mt = mt.astype(float)
    def calc_slope(row):
        a = scipy.stats.linregress(axisvalues_ , row)
        return pd.Series(a._asdict())
    
    print (mt.apply(calc_slope,axis=1))
    
    mt = mt.join(mt.apply(calc_slope,axis=1))
    return(mt)


# Run once for serum samples 
mts = stat(loc3)
mts = mts.reset_index()
mts = pd.concat([loc1, mts], axis=1)

#Benjamini-Hochberg filtering p-values
mts = mts.sort_values(by='pvalue') #pvalue calculated in this code?
mts['rank'] = mts.reset_index().index + 1
mts['(I/m)Q'] = (mts['rank']/len(mts))*0.05
mts['(I/m)Q - p'] = mts['(I/m)Q'] - mts['pvalue']
mts['BH_sig'] = np.where(mts['(I/m)Q - p'] < 0, '0', '1')

mts.to_csv('Neg_RPL_serum_Stats05_PMM.csv')

# isolate the features that show a significant relationship with PTB
mts = mts[mts['pvalue'] < 0.05]
mts.to_csv('Neg_RPL_serum_Stats05-2_PMM.csv')


