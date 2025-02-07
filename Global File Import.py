'''
imports files into a pandas dataframe
'''

import os, glob
import pandas as pd

path = r'INSERT YOUR PATH HERE\\'

files = glob.glob(path + '*.csv')

df = pd.concat(
    (
        pd.read_csv(i)          #add appropriate formatting to the read_csv command
                for i in files
    ), 
    axis=1
)

###################################################################################
Here's the same command with boxcar averaging applied to each file.
###################################################################################

df = pd.concat(
    (
        pd.read_csv(i, index_col=0, skiprows=[0,1],header=None).rolling(window=4,win_type='boxcar').mean() 
                for i in files
    ), 
    axis=1
)
