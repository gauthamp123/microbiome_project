# Take the results.tsv file and parse it
import pandas as pd
import csv
import numpy as np
import os
import re
# from fusion_distribution import isFusion,geneFusions
# from fusion_dist_dir import isFusion,genDict
from fusion import isFusion, geneFusions, genDict
from parseXML import parse
from pprint import pprint

GENOME = 'GCF_009648975.1'  # we can make this a command line argument
# construct df with all data from results.tsv
df = pd.read_table(GENOME + '/results.tsv')
# print columns of df for dev use in constructing filtered_df later


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

def adjustOverlapScore():
    overlap_dict = parse(GENOME + '/hmmtop.db', GENOME + '/xml/' ,GENOME + '/results.tsv', 8)
    score_dict = {}
    score_dict['TM_Overlap_Score'] = []
    for k in overlap_dict:
        score_dict['TM_Overlap_Score'].append(overlap_dict[k]['alignedTMS'])
    
    score_df = pd.DataFrame(score_dict)
    df['TM_Overlap_Score'] = score_df['TM_Overlap_Score']




adjustOverlapScore()


green = {}
yellow = {}
red = {}
# To distinguish between single and multiple component systems. 
# Example content:
#  1.A.1.1.1 =>  ["1.A.1.1.1-P0A334"],
#  3.A.1.1.1 =>  ["3.A.1.1.1-P02916", ""3.A.1.1.1-XXXXX", "3.A.1.1.1-YYYYYY", "3.A.1.1.1-ZZZZZZZZZ"],
#  ....
tcdbSystems = {}

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

#This is a helper function for finding common pFam domains and can be used to check if a value is a float
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
    
def qCoverage(row):
    return float(row.get(key='Query_Coverage'))

def hCoverage(row):
    return float(row.get(key='Hit_Coverage'))

def eVal(row):
    return float(row.get(key='e-value'))

def pfamDoms(row):
    doms = []
    if isfloat(row.get(key='Query_Pfam')):
        return doms
    elif isfloat(row.get(key='Subject_Pfam')):
        return doms

    q_pfam = row.get(key='Query_Pfam').split(',')
    s_pfam = row.get(key='Subject_Pfam').split(',')


    for q_domain in q_pfam:
        for s_domain in s_pfam:
            if q_domain == s_domain:
                doms.append(q_domain)
    common_doms = [*set(doms)]
    return common_doms


parseTCDBcontent()

# This function will check the tcdbSystems dictionary for the tcid given 
def isSingleComp(row):
    tcid = row["Hit_tcid"]
    hit_id = row['Hit_xid']
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) == 1):
        return True
    else:
        return False
    
def categorizeSingleComp(row):
    # check if it is a fusion
    fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    if(len(fusion_results) !=1):
        sortedGeneArr = sorted(fusion_results, key=lambda x: x['sstart'])
        if(len(isFusion(sortedGeneArr))!=0):
            row_to_write = row.tolist()
            row_to_write.append(True)
            return ('yellow', row_to_write)
    
    row_to_write = row.tolist()
    row_to_write.append(False)
    # check evalue
    if eVal(row) <= 1e-10 and len(pfamDoms(row)) != 0:
        green[row.get(key='#Query_id')] = row_to_write
        return ('Green', row_to_write)
    elif eVal(row) < 1e-3 and qCoverage(row) >= 90 and hCoverage(row) >= 90:
        yellow[row.get(key='#Query_id')] = row_to_write
        return ('Yellow', row_to_write)
    else:
        red[row.get(key='#Query_id')] = row_to_write
        return ('Red', row_to_write)
    

    
with open('hmmtop.out') as f:
    lines = f.readlines()
    hmmtop_df = pd.DataFrame(columns=['Hit_tcid', 'Hit_xid', 'Hit_n_TMS','Match_length'])
    ##Create a dataframe that include above columns 
    for line in lines:
        fields = re.split(r'[ -]+', line.strip())
        ##split the system and protein names
        new_col=[fields[2],fields[3],fields[5]]
        new_row = pd.Series([fields[2],fields[3],fields[5],fields[1]], index=hmmtop_df.columns)
        hmmtop_df.loc[len(hmmtop_df)] = new_row
    hmmtop_df['Hit_n_TMS'] = hmmtop_df['Hit_n_TMS'].astype(int)

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
       if(eVal(row) <= float("1e-3") and ((qCoverage(row) <75 and hCoverage(row) <75) or len(Fusion_Add)!=0 )):
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
    
        
        

def write_singlecomp(Output_dict,Output_df_row):
    color = output_dict[0]
    dictionary = output_dict[1]
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = dictionary[len(dictionary) - 1]
    missing_proteins = "NA"
    Intermediate['Missing_components']= missing_proteins
    filename=f"{color}.tsv"

    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)
    




geneFusions={}
genDict(geneFusions, df)

for filename in ["Green.tsv","Red.tsv","Yellow.tsv"]:
    if os.path.exists(filename):
        os.remove(filename)
for index, row in df.iterrows():


    Output_dict= isMultiComp(row, df, 0.5)

    if(Output_dict):
        Write_multicomp(Output_dict,Output_df.loc[[index],Output_df.columns])
    else:
        output_dict = categorizeSingleComp(row)
        write_singlecomp(output_dict, Output_df.loc[[index],Output_df.columns])


    

    #print(Output_dict)
