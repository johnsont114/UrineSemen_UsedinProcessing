# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 17:55:14 2024

@author: JOHNST47
"""

import pandas as pd
import numpy as np

# Load your spreadsheet
df = pd.read_csv('updated_spreadsheet_NaNreplaced.csv')

# Replace all  with NaN
df.replace(0, np.nan, inplace=True)
df[df<1000] = np.NaN

# Save your updated spreadsheet
df.to_csv('updated_spreadsheet_NaNs.csv', index=False)