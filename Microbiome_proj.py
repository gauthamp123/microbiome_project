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
    ##Create a dataframe that include above columns 
    for line in lines:
        fields = re.split(r'[ -]+', line.strip())
        ##split the system and protein names
        new_col=[fields[2],fields[3],fields[5]]
        new_row = pd.Series([fields[2],fields[3],fields[5]], index=hmmtop_df.columns)
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
    if(set(tc_all_arr)==set(tc_filter_arr)):##If all the proteins in that system can be found, then green
        
        return("green",tc_filter_arr)
    MembraneProteins= FindMembraneProtein(row, df)
    if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
        ##given some proteins can be found while containing the membrane proteins 
        
        return("Yellow", tc_filter_arr)
    else:
        print(tc_arr)
        print("MP",MembraneProteins)
        print("Red",tc_filter_arr)
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

#This is a helper function for finding common pFam domains and can be used to check if a value is a float
def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def qCoverage(row):
    return float(row[7])

def hCoverage(row):
    return float(row[8])

def eVal(row):
    return float(row[2])

def PfamDoms(row):
    doms = []
    if isfloat(row[10]):
        return doms
    elif isfloat(row[11]):
        return doms
        
    q_pfam = row[10].split(',')
    s_pfam = row[11].split(',')


    for q_domain in q_pfam:
        for s_domain in s_pfam:
            if q_domain == s_domain:
                doms.append(q_domain)
    common_doms = [*set(doms)]
    return common_doms


parseTCDBcontent()

for index, row in df.iterrows():
    print(isSingleComp(row))

'''
for index, row in df.iterrows():
    if(isSingleComp(row)):
        if(eVal(row) <= -10 and pfamDoms().size != 0):
            green_df = green_df.append([row])
        elif(eVal(row) < -3 and qCoverage() >= 90 and hCoverage >= 90):
            yellow_df = yellow_df.append([row])
        else:
            red_df = red_df.append([row])

'''
