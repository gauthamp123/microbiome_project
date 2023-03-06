#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

import pandas as pd
import numpy as np
import os

# construct df with all data from results.tsv
df = pd.read_table('results.tsv')
# print columns of df for dev use in constructing filtered_df later
print('---------')
print(df.columns)
print('---------')

# create empty dfs for tagging with green, yellow, and red labels
green_df = df.copy()
green_df = green_df.iloc[0:0]
yellow_df = green_df.copy()
red_df = green_df.copy()

# construct filtered_df with only relevant data points for each protein matche
filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
#print(filtered_df)

#Example content:
#  1.A.1.1.1 =>  ["1.A.1.1.1-P0A334"],
#  3.A.1.1.1 =>  ["3.A.1.1.1-P02916", ""3.A.1.1.1-XXXXX", "3.A.1.1.1-YYYYYY", "3.A.1.1.1-ZZZZZZZZZ"],
#  ....
tcdbSystems = {}

# TODO

input = open("tcdb.faa")
tcdbSystems = {}

#This function will fill the dictionary tcdbSystems with data.
def parseTCDBcontent():
    for line in input:
        if(">" in line):
            first = line.split("-")[0][1:]
            if(first in tcdbSystems):
                tcdbSystems.get(first).append(line.strip(">\n"))
            else:
                tcdbSystems[first] = [line.strip(">\n")]

# This function will check the tcdbSystems dictionary for the tcid given 
def isSingleComp(row):
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) == 1):
        return True
    else:
        return False

def FindMembrance(row,df):
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) == 1):
        return [row["Hit_xid"]]
    elif(len(tc_arr) == 0):
        return []
    elif(len(tc_arr) > 1):
        ##print(tc_arr)
        tc_filter_arr= [arr.split("-")[1] for arr in tc_arr]
        tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_filter_arr))
        final_df = pd.concat([ df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])
        max_value = final_df['Hit_n_TMS'].max()
        result = final_df.loc[final_df['Hit_n_TMS'] == max_value]
        ##print(row["Hit_xid"],tc_arr)
        ##pd.set_option('display.max_colwidth', None)
        
        Final_result= result['Hit_xid'].drop_duplicates().tolist()
        print(Final_result)
        return(Final_result)
        ##if(len(Final_result)==3):
            ##print(Final_result)
            ##print(tc_arr)
            ##print(tc_filter_arr)
        
      
    
def isMultiComp(row,df):
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) >1):
        tc_filter_arr= [arr.split("-")[1] for arr in tc_arr]
        tc_filter_arr.remove(row["Hit_xid"])
        for arr in tc_filter_arr:
            if(len(df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}' "))) ==0 :
                ##print(arr,row["Hit_xid"],tc_arr)
                return False

        ##print(row["Hit_xid"],tc_arr)
        return True   
    else:
        return False

def qCoverage(row):
    return row["Query_Coverage"]

def hCoverage(row):
    return row["Hit_Coverage"]

def eVal(row):
    return row["e-value"]

def PfamDoms(row):
    common_doms = ""
    return common_doms

parseTCDBcontent()

##for index, row in df.iterrows():
    ##print(isSingleComp(row))
'''    
GREEN (Best hits)
   a) High coverage (e.g., >=75% in both proteins; set as a threshold variable)
   b) low E-value (e.g., <= 1e-10)
   c) shared Pfam domains
   d) If the value is larger (e.g., > 1e-10) but there is reasonable coverage (e.g. > 50%) in both proteins AND they have shared domains.

   YELLOW (not so sure)
   a)lower coverage (e.g. <75% in both proteins)
   b)Ok value (< 1e-3)
   c)there are no shared domains
   d) if one protein has very low coverage (e.g., ~10%) and there is other protein in the genome matching other region of the TCDB protein, select this proteins as potential fusions. 
   e) 2 or more proteins in TCDB can fuse to form a single component system in your reference genome.
      ** d and e required the coordinates that Clark is going to add to the file.


   RED (no good)
   a) One protein has very low coverage (e.g. ~10%), there are no common domains AND there are no other proteins in the genome matching the same protein in TCDB.
'''

for index, row in df.iterrows():
    ##isMultiComp(row,df)
    FindMembrance(row,df)
    if(isSingleComp(row)):
        if(eVal(row) <= float("1e-10") and qCoverage(row) >= 75 and hCoverage(row) >= 75):
            ##green_df = green_df.append([row])
            green_df = pd.concat([green_df,row.to_frame().T])
        elif(eVal(row) <= float("1e-3")and qCoverage(row) <75 and hCoverage(row) <75):
            ##yellow_df = yellow_df.append([row])
            yellow_df = pd.concat([yellow_df,row.to_frame().T])
        elif(qCoverage(row) <=10 and hCoverage(row) <=10):
            ##red_df = red_df.append([row])
            red_df = pd.concat([red_df,row.to_frame().T])
        
    '''
        if(eVal(row) <= -10 and pfamDoms().size != 0):
            green_df = green_df.append([row])
        elif(eVal(row) < -3 and qCoverage() >= 90 and hCoverage >= 90):
            yellow_df = yellow_df.append([row])
        else:
            red_df = red_df.append([row])
    '''
###print(green_df)
