#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

import pandas as pd
import numpy as np
import os
import re
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

with open('hmmtop.out') as f:
    lines = f.readlines()
    hmmtop_df = pd.DataFrame(columns=['Hit_tcid', 'Hit_xid', 'Hit_n_TMS'])
    for line in lines:
        fields = re.split(r'[ -]+', line.strip())
        new_col=[fields[2],fields[3],fields[5]]
        new_row = pd.Series([fields[2],fields[3],fields[5]], index=hmmtop_df.columns)
        hmmtop_df.loc[len(hmmtop_df)] = new_row
    hmmtop_df['Hit_n_TMS'] = hmmtop_df['Hit_n_TMS'].astype(int)

def FindMembraneProtein(row,df):
    if isSingleComp(row):
        return ([])
    ##get rid of single comp system first
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    if tcid.split(".")[1]=="B":
        return (tc_all_arr)
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_Notfind_arr=list(set(tc_all_arr)-set(tc_filter_arr))
    find_df = pd.concat([df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])
    unfind_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])

    if find_df['Hit_n_TMS'].nunique() == 1 and unfind_df['Hit_n_TMS'].nunique() == 1 :
        if str(find_df['Hit_n_TMS'].unique()[0]) == str(unfind_df['Hit_n_TMS'].unique()[0]):
            
            return (tc_all_arr)
    
    Found_arr_gblast = find_df[(find_df['Hit_n_TMS'] >= 3) & (abs(find_df['Hit_n_TMS'] - find_df['Query_n_TMS']) <= 2)]["Hit_xid"].tolist()
    tcdb_arr= Found_arr_gblast +tc_Notfind_arr
    if len(tcdb_arr)==0:
        return([])
    tcdb_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tcdb_arr])
    Final_return = tcdb_df[(tcdb_df['Hit_n_TMS']>= 3)]  
    ##print(tc_all_arr)
    ##print(Final_return)
    return(Final_return["Hit_xid"].tolist())  
    ##final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])
'''
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    ##This is the actual proteins we find in the gblast result
    ##tc_Notfind_arr=list(set(tc_all_arr)-set(tc_filter_arr))
    final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])
    ##final df is the whole line that tc_filter_arr locate
    if final_df.empty:
        return ([])
    
    if final_df['Hit_n_TMS'].nunique() == 1: ##and str(final_df['Hit_n_TMS'].unique()[0]) != "0":
        ##This means if all the proteins have same hit TMS in their system
        if abs(final_df['Hit_n_TMS'].iloc[0] - final_df['Query_n_TMS'].iloc[0]) <= 2: ##If the Query TMS has less than 3 difference
            
            return(final_df)
    ##Otherwise we check the protein that has TMS greater than 2, and if the systems contain such protein and has
    ##query TMS difference less than 3, then we say they are membrane protein
    filtered_df = final_df[(final_df['Hit_n_TMS'] >= 3) & (abs(final_df['Hit_n_TMS'] - final_df['Query_n_TMS']) <= 2)]
    ##print(tc_arr)
    ##print(final_df)
    ##print(filtered_df[['Query_n_TMS','Hit_n_TMS']])
    return(filtered_df)
'''        
      
'''    
def isMultiComp(row,df,int input):
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
'''
def isMultiComp3(row,df,input):
    if isSingleComp(row):
        return (False)
    ##get rid of single comp system first
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    if(set(tc_all_arr)==set(tc_filter_arr)):
        return("green",tc_filter_arr)
    MembraneProteins= FindMembraneProtein(row, df)
    if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
        print(tc_all_arr)
        print("MP",MembraneProteins)
        print("Yellow",tc_filter_arr)
        return("Yellow", tc_filter_arr)
    else:
        return ("Red", tc_filter_arr)


def isMultiComp(row,df):
    ##result=FindMembraneProtein(row,df)
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    tc_arr = list(set(tc_arr))
    if(len(tc_arr) >1):
        tc_filter_arr = [arr.split("-")[1] for arr in tc_arr]
        mutiComp_arr=[]
        mutiComp_arr_xid=[]
        for arr in tc_filter_arr:
            query_result_df=df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}' ")
            query_result_list=query_result_df.values.tolist()
            query_result_xid=query_result_df["Hit_xid"].tolist()
            if(len(query_result_list)) !=0 :
                mutiComp_arr += query_result_list
                mutiComp_arr_xid += query_result_xid
        ##print(tc_filter_arr,query_result_xid)
        if set(tc_filter_arr) == set(mutiComp_arr_xid):
            ##print("green")
            ##print(tc_filter_arr,mutiComp_arr_xid)
            return("green", mutiComp_arr)
        elif(len(mutiComp_arr)==0):
            return("red",mutiComp_arr)
        else:
            tc_filter_arr_list= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_filter_arr))
            mutiComp_arr_df = pd.concat([ df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr_list]).drop_duplicates()
            if mutiComp_arr_df.loc[mutiComp_arr_df['Hit_n_TMS'] >=3].values.tolist() != []:
                return("yellow",mutiComp_arr_df.values.tolist())
            else:
                return("red",mutiComp_arr_df.values.tolist())

def isMultiComp2(row,df):
    ##result=FindMembraneProtein(row,df)
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    tc_arr = list(set(tc_arr))
    if(len(tc_arr) >1):
        tc_filter_arr = [arr.split("-")[1] for arr in tc_arr]
        mutiComp_arr_df=pd.DataFrame(columns=df.columns)
        mutiComp_arr_xid=[]
        for arr in tc_filter_arr:
            query_result_df=df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}' ")
            query_result_xid=query_result_df["Hit_xid"].tolist()
            if(len(query_result_df)) !=0 :
                ##print(mutiComp_arr_df,query_result_df)
                mutiComp_arr_df = pd.concat([mutiComp_arr_df,query_result_df])
                mutiComp_arr_xid += query_result_xid
        if set(tc_filter_arr) == set(mutiComp_arr_xid):
            return("green", mutiComp_arr_df)
        elif(len(mutiComp_arr_df)==0):
            return("red",mutiComp_arr_df)
        else:
            tc_filter_arr_list= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_filter_arr))
            mutiComp_arr_df = pd.concat([ df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr_list]).drop_duplicates()
            if len(mutiComp_arr_df.loc[mutiComp_arr_df['Hit_n_TMS'] >=3]) != 0:
                return("yellow",mutiComp_arr_df)
            else:
                print("red",mutiComp_arr_df)

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
    isMultiComp3(row,df,0.5)
   
   
    ##FindMembraneProtein(row,df)
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
            
   
   
###print(green_df)
