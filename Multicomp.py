#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

import pandas as pd
import numpy as np
import os
import re
from fusion_dist_dir import isFusion,genDict
from pprint import pprint
# construct df with all data from results.tsv
df = pd.read_table('data/results.tsv')
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

Output_df= df[['Hit_tcid','Hit_xid','#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]
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
    hmmtop_df = pd.DataFrame(columns=['Hit_tcid', 'Hit_xid', 'Hit_n_TMS','Match_length'])
    ##Create a dataframe that include above columns 
    for line in lines:
        fields = re.split(r'[ -]+', line.strip())
        ##split the system and protein names
        new_col=[fields[2],fields[3],fields[5],fields[1]]
        new_row = pd.Series([fields[2],fields[3],fields[5],fields[1]], index=hmmtop_df.columns)
        hmmtop_df.loc[len(hmmtop_df)] = new_row
    hmmtop_df['Hit_n_TMS'] = hmmtop_df['Hit_n_TMS'].astype(int)
#print(hmmtop_df)
def FindMembraneProtein(row,df):
    if isSingleComp(row): ##get rid of single comp system first
        return ([])
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)##This contains all the proteins with their system name

    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    
    ##if tcid.split(".")[1]=="B":##and 
    if tcid[0:3]=="1.B":##If there is a beta barrel, we assume they are membrane proteins
        return (tc_all_arr)
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    ##tc_filter_arr includes the proteins that showed up in the gblast result
    tc_Notfind_arr=list(set(tc_all_arr)-set(tc_filter_arr))##This is the missing proteins in actual result
    find_df = pd.concat([df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])
    ##This is the whole line for tc_filter_arr in gblast result
    unfind_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tc_all_arr])
    ##This is the whole line for tc_all_arr in hmmtop file
    if find_df['Hit_n_TMS'].nunique() == 1 and unfind_df['Hit_n_TMS'].nunique() == 1 :
        ##If all the proteins have same Hit TMS in both files, we say they are membrane proteins
        if str(find_df['Hit_n_TMS'].unique()[0]) == str(unfind_df['Hit_n_TMS'].unique()[0]):
           
            return (tc_all_arr)
    Found_arr_gblast = find_df[(find_df['Hit_n_TMS'] >= 3) & (abs(find_df['Hit_n_TMS'] - find_df['Query_n_TMS']) <= 2)]["Hit_xid"].tolist()
    ##If there are proteins that has Hit TMS >=3 and difference with its query TMS <=2, we say it's membrane proteins
    tcdb_arr= Found_arr_gblast +tc_Notfind_arr##This arr contains the possible output proteins 
    if len(tcdb_arr)==0:
        return([])
    tcdb_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tcdb_arr])
    Final_return = tcdb_df[(tcdb_df['Hit_n_TMS']>= 3)]  
    ##print(tc_all_arr)
    ##print(find_df)
    ##print(unfind_df)
    ##print(Final_return["Hit_xid"].tolist())
    return(Final_return["Hit_xid"].tolist())  
    ##final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])


def isMultiComp(row,df,input):
    if isSingleComp(row):
        return (False)
    ##get rid of single comp system first
    tcid = row["Hit_tcid"]
    Fusion_Add=[]
    Fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    if(len(Fusion_results) !=1):
        sortedGeneArr = sorted(Fusion_results, key=lambda x: x['sstart'])
        Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]             
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))

    if(set(tc_all_arr)==set(tc_filter_arr)):##If all the proteins in that system can be found, then green
        #print(tc_arr)
        #print("green",tc_filter_arr,"Fusion Results:",Fusion_Add)
       
        if(eVal(row) <= float("1e-10") and (len(Fusion_Add)!=0 or (qCoverage(row) >= 75 and hCoverage(row) >= 75))):
            return({"color":"Green",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0})

    MembraneProteins= FindMembraneProtein(row, df)
    # if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
    if(input*len(tc_all_arr)<=len(tc_filter_arr)) or len(set(MembraneProteins) & set(tc_filter_arr))>0:
        ##given some proteins can be found while containing the membrane proteins 
       if(eVal(row) <= float("1e-3") and ((qCoverage(row) >=50 and hCoverage(row) >=50) or len(Fusion_Add)!=0 )):
            #print("Yellow")
            return({"color":"Yellow",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0})

    
    #print(tc_arr)
    #print("MP",MembraneProteins)    
    #print("Red",tc_filter_arr)
    #print(("Red", tc_filter_arr,"Fusion Results:",Fusion_Add,eVal(row),hCoverage(row)))
    return({"color":"Red",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0})
   


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
'''
def Write_multicomp(Output_dict,Output_df_row):
    Intermediate=Output_df_row.copy()
    Intermediate['isFusion']=Output_dict['isFusion']
    missing_proteins= "NA" if(len(Output_dict['Missing_proteins'])==0) else ",".join(Output_dict['Missing_proteins'])
    Intermediate['Missing_components']= missing_proteins
    filename=f"{Output_dict['color']}.tsv"
    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)
'''
def Write_multicomp(Output_dict,Output_df_row):
    Intermediate=Output_df_row.copy()
    Intermediate['isFusion']=Output_dict['isFusion']
    filename=f"{Output_dict['color']}.tsv"
    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)

    for hit_xid in Output_dict['Missing_proteins']:
        _Intermediate=Output_df_row.copy()
        _Intermediate=_Intermediate.applymap(lambda x: 'NA')
        #print(_Intermediate)

        _Intermediate["Hit_tcid"]=Output_df_row["Hit_tcid"]
        _Intermediate["Hit_xid"]=hit_xid
        Missing_infor = hmmtop_df.loc[(hmmtop_df['Hit_xid'] ==  hit_xid)]

        _Intermediate["Match_length"]=Missing_infor["Match_length"].iloc[0] if not Missing_infor.empty else "NA"
        _Intermediate["Hit_n_TMS"]=Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
        with open(filename, mode="a", encoding='utf-8') as f:
            _Intermediate.to_csv(f, sep='\t', header=False, index=False)
            
geneFusions={}
genDict(geneFusions, df)
for filename in ["Green.tsv","Red.tsv","Yellow.tsv"]:
    if os.path.exists(filename):
        os.remove(filename)
for index, row in df.iterrows():

    Output_dict= isMultiComp(row, df, 0.5)
    if(Output_dict):
        Write_multicomp(Output_dict,Output_df.loc[[index],Output_df.columns])
    #print(Output_dict)

for filename in ["Green.tsv","Red.tsv","Yellow.tsv"]:
    if os.path.exists(filename):
        df = pd.read_csv(filename, sep='\t')
        df = df.fillna('NA')
        NA_count = df.eq('NA').sum(axis=1).rename('NA_count')
        df['na_sort'] = NA_count
        df = df.sort_values(by=['Hit_tcid','na_sort']).drop(columns=['na_sort'])
        df.to_csv(filename, sep='\t', index=False)