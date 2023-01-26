#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

import pandas as pd
import numpy as np

# construct df with all data from results.tsv
df = pd.read_table('results.tsv')
# print columns of df for dev use in constructing filtered_df later
print('---------')
print(df.columns)
print('---------')

# construct filtered_df with only relevant data points for each protein matche
filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
#print(filtered_df)

for index, row in filtered_df.iterrows():
    print(row['#Query_id'])


