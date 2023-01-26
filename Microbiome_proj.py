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

green_df = df.copy()
green_df = green_df.iloc[0:0]
yellow_df = green_df.copy()
red_df = green_df.copy()

# construct filtered_df with only relevant data points for each protein matche
filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
#print(filtered_df)

# TODO
def isSingleComp(row):
    return True

def qCoverage(row):
    return 0

def hCoverage(row):
    return 0

def eVal(row):
    return 0

def pfamDoms(row):
    common_doms = ""
    return common_doms

'''
for index, row in df.iterrows():
    if(isSingleComp(row)):
        if(eVal(row) <= -10 and pfamDoms().size != 0):
            green_df = green_df.append([row])
        elif(eVal(row) <= -5 and qCoverage() >= 90 and hCoverage >= 90):
            yellow_df = yellow_df.append([row])
        else:
            red_df = red_df.append([row])

'''
