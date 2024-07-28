# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:57:20 2024

@author: Abrahd05
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

df1 = pd.read_csv('10000Max_NonAvg+NaN_Neg.csv')

df1a = df1.loc[:, 'Alignment ID']
df1b = df1.loc[:, 'Cornell_Pos_S-01_LMT_PS_PMM_R0':'Cornell_Pos_S-50_LMH_LCT_PS_PMM_R1']
#df1c = df1.loc[:, 'Average_Rt(min)':'MS/MS spectrum']

df1bL = np.log10(df1b)

def fillNaN_with_unifrand(df):
    lower, upper = 0, df.min()
    a = df.values
    m = np.isnan(a) # mask of NaNs
    mu, sigma = df.min(), df.std()
    a[m] = stats.truncnorm.rvs(
          (lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=m.sum())
    return df


df1bLmod = df1bL.apply(fillNaN_with_unifrand)

df1bLmod = 10**df1bLmod
#df1 = pd.concat([df1a, df1bLmod, df1c], axis=1)
df1bLmod.to_csv('Imputed_SemUr_Neg_CleanFna_All.csv')

