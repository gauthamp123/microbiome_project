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

GENOME = 'MicrobiomeResults/GCF_000013425.1'  # we can make this a command line argument
TMS_THRESH_GREEN = 0.8
TMS_THRESH_YELLOW = 0.5
EXCELLENT_EVAL = 1e-30
E_VAL_GREEN = 1e-10
E_VAL_YELLOW = 1e-3 
Q_COV_THRESH = 80
S_COV_THRESH = 80
LOW_COV = 50
AUTO_RED = 20


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



def adjustOverlapScore():
    overlap_dict = parse(GENOME + '/hmmtop.db', GENOME + '/xml/' ,GENOME + '/results.tsv', 8)
    score_dict = {}
    score_dict = {}
    for k in overlap_dict:
        score_dict[k] = overlap_dict[k]['alignedTMS']
    
    for index, row in df.iterrows():
        if row['#Query_id'] in score_dict:
            df.at[index, 'TM_Overlap_Score'] = score_dict[row['#Query_id']]
        else:
            df.at[index, 'TM_Overlap_Score'] = 0
            
    


adjustOverlapScore()

Output_df= df[['Hit_tcid','Hit_xid','#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]

multi_df = df[['Hit_xid', '#Query_id', 'Hit_tcid', 'Match_length', 'e-value', '%_identity', 'Query_Length', 'Hit_Length', 'Q_start',
               'Q_end', 'S_start', 'S_end', 'Query_Coverage', 'Hit_Coverage', 'Query_n_TMS', 'Hit_n_TMS', 'TM_Overlap_Score', 'Family_Abrv',
               'Predicted_Substrate', 'Query_Pfam', 'Subject_Pfam']]

for index, row in multi_df.iterrows():
    multi_df.at[index, 'Hit_xid'] = row['Hit_tcid'] + '-' + row['Hit_xid']

multi_df = multi_df.groupby('Hit_xid').first()

multi_dict = multi_df.to_dict('index')

df.to_csv('adj.csv', index=False)

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

    if len(q_pfam) == 0 and len(s_pfam) == 0:
        return ['N/A']
    for q_domain in q_pfam:
        for s_domain in s_pfam:
            if q_domain == s_domain:
                doms.append(q_domain)
    common_doms = [*set(doms)]
    return common_doms

def tm_overlap_decision(row):
    q_tmss = int(row.get(key='Query_n_TMS'))
    t_tmss = int(row.get(key='Hit_n_TMS'))
    # if there are no tms regions in the tcdb protein tms regions are irrelevant in decision making
    if t_tmss == 0:
        return False
    # calculates the percent overlap of tms regions relative to tcdb proteins
    def tmoverlap_percent(row):
        overlap_percent = int(row.get(key='TM_Overlap_Score')) / t_tmss
        return overlap_percent

    # if there are under 3 tmss then theres a possibility of there being mishits
    if t_tmss > 3:
        # if there is great overlap return true
        if tmoverlap_percent(row) >= TMS_THRESH_GREEN:
            return True
    else:
        # could potentially have a fusion candid so should be placed into yellow
        if q_tmss >= 1 and t_tmss >= 1:
            return True
    
    return False
        
    
def overlap_percent(row):
    q_tmss = int(row.get(key='Query_n_TMS'))
    t_tmss = int(row.get(key='Hit_n_TMS'))
    if q_tmss == 0 or t_tmss == 0:
        return 0

    overlap_percent = int(row.get(key='TM_Overlap_Score')) / max(q_tmss, t_tmss)
    return overlap_percent


def assign_if_fusion(fusions, tcdb_tms):
    # need to check that the combined tms count from all fusions meets threshold if ritvik can add to output
    tot_cov = 0
    tot_tmoverlap = 0
    for i in range(0, len(fusions)):
        tot_cov = fusions[0]['send'] - fusions[0]['sstart']
        for i in range(1, len(fusions)):
            start = 0
            end = fusions[i]['send'] 
            if fusions[i]['sstart'] < fusions[i-1]['send']:
                start = fusions[i-1]['send']
            else:
                start = fusions[i]['sstart']
            
            if fusions[i]['sstart'] <= fusions[i-1]['sstart'] and fusions[i]['send'] <= fusions[i-1]['send']:
                return 'yellow'

            tot_cov += end - start
    
    if tcdb_tms == 0:
        tot_tmoverlap = 1
    else:
        tot_tmoverlap = float(tot_tmoverlap / tcdb_tms)

    tot_cov = float(tot_cov / fusions[0]['hit_length']) * 100
        

    if tot_cov >= S_COV_THRESH and tot_tmoverlap >= TMS_THRESH_GREEN:
        return 'green'
    elif tot_cov >= LOW_COV:
        return 'yellow'
    else:
        return 'red'
    


def make_decision(eval, qcov, scov, pfam_doms, fusions, tmoverlap, tcdb_tms):
    # automatic green if e-val is 0 or extremely good with okay coverage and matching pfam
    if (eval == 0.0 or eval <= EXCELLENT_EVAL) and (qcov > LOW_COV and scov > LOW_COV) and len(pfam_doms) > 0:
        return 'green'
    # if coverage is less than a very low number it is impossible for it to be a good hit
    if qcov <= AUTO_RED or scov <= AUTO_RED:
        return 'red'
    # if there are less than 3 tms, there is a possiblity of mischaracterizing tms regions
    if tcdb_tms > 3:
        # if cov is too low it is automatically red
        if qcov < LOW_COV or scov < LOW_COV:
            return 'red'
        # great e val, there are common pfam domains, good cov and high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and qcov >= Q_COV_THRESH and scov >= S_COV_THRESH and tmoverlap >= TMS_THRESH_GREEN:
            return 'green'
        # great e val, common doms, has good coverage or there is a high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # great e value, has good coverage or there is a high tms overlap means good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # okay e value, common doms, has good coverage or there is a high tms overlap means good hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_YELLOW):
            return 'green'
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_GREEN and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        # ok e-val, no common doms, has good coverage and there is a ok tms overlap means yellow hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and tmoverlap >= TMS_THRESH_YELLOW):
            return 'yellow'
        # okay e val, and com dom tms overlap is good and the coverage is not too terrible > 50% its a yellow hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and (qcov > LOW_COV or scov > LOW_COV) and tmoverlap >= TMS_THRESH_YELLOW:
            return 'yellow'
        # if low e-val, has low coverage and there are low tms overlap
        if eval <= E_VAL_YELLOW and (qcov < LOW_COV or scov < LOW_COV) and tmoverlap < TMS_THRESH_YELLOW:
            return 'red'
    else: # considering possibility of mischaracterization of tms
        # is a fusion candidate, has great e value and has good cov is good hit
        if eval <= E_VAL_GREEN and (qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and len(pfam_doms) > 0:
            return 'green'
        # good e value, common doms match, good cov is good hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0 and (qcov >= Q_COV_THRESH and scov >= S_COV_THRESH):
            return 'green'
        # has ok e value, common doms match, good cov is good hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and (qcov >= LOW_COV and scov >= LOW_COV):
            return 'green'
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_GREEN and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        # is fusion, ok e value and good coverage but not good tms overlap
        if eval <= E_VAL_YELLOW and qcov >= Q_COV_THRESH and scov >= S_COV_THRESH:
            return 'yellow'
        # great e val, and com dom but not good coverage is yellow hit
        if eval <= E_VAL_GREEN and len(pfam_doms) > 0:
            return 'yellow'
        # okay e val, and com dom but not good coverage is yellow hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0:
            return 'yellow'
        # okay e val, no common doms, good tms overlap
        if eval <= E_VAL_YELLOW and tmoverlap >= TMS_THRESH_GREEN:
            return 'yellow'
        
    return 'red'

def final_decision(eval, qcov, scov, pfam_doms, fusions, tmoverlap, tcdb_tms):
    decision = make_decision(eval, qcov, scov, pfam_doms, fusions, tmoverlap, tcdb_tms)
    if decision == 'fusion found':
        return assign_if_fusion(fusions, tcdb_tms)
    else:
        return decision

parseTCDBcontent()

# This function will check the tcdbSystems dictionary for the tcid given 
def isSingleComp(row):
    tcid = row["Hit_tcid"]
    tc_arr = tcdbSystems.get(tcid)
    if(len(tc_arr) == 1):
        return True
    else:
        return False
    
def categorizeSingleComp(row):
    # check if it is a fusion
    row_to_write = row.tolist()
    sortedGeneArr = []
    fusions = ''
    fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    if(len(fusion_results) !=1):
        sortedGeneArr = sorted(fusion_results, key=lambda x: x['sstart'])
        if(len(isFusion(sortedGeneArr))!=0):
            sortedGeneArr = isFusion(sortedGeneArr)
            for k in range(0, len(sortedGeneArr)):
                fusions += sortedGeneArr[k]['query'] + ' '
                
            row_to_write.append(fusions)
        else:
            row_to_write.append(fusions)
    else:
        row_to_write.append(fusions)
    

    tcdb_tms = row['Hit_n_TMS']

    if final_decision(eVal(row), qCoverage(row), hCoverage(row), pfamDoms(row), sortedGeneArr, overlap_percent(row), tcdb_tms) == 'green':
        green[row.get(key='#Query_id')] = row_to_write
        return ('Green', row_to_write)
    elif final_decision(eVal(row), qCoverage(row), hCoverage(row), pfamDoms(row), sortedGeneArr, overlap_percent(row), tcdb_tms) == 'yellow':
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

# returns empty list when there are no membrane proteins found. 
# otherwise returns list of membrane protein accessions.
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
    return((Final_return['Hit_tcid'] + '-' + Final_return["Hit_xid"]).tolist())  

    # return Final_return.to_dict()
    ##final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])

def multicomp_decision(tcdb_proteins, genome_missing_proteins, membrane_proteins):
    tcid = tcdb_proteins[0].split('-')[0]
    # max_protein = {'Protein': 0}
    # temp_key = 'Protein'
    # for protein in tcdb_proteins:
    #     curr_row = pd.DataFrame([multi_dict[protein]])
    #     ttms = int(row.get(key='Hit_n_TMS'))
    #     if ttms > max_protein[temp_key]:
    #         max_protein = {protein: ttms}
    #         temp_key = protein
        

    # check to see how many missing comps and what is missing
    num_membrane_proteins = len(membrane_proteins)
    if num_membrane_proteins == 0:
        return 'No membrane proteins'
    
    membrane_accessions = []
    for i in range(0, num_membrane_proteins):
        if '-' in membrane_proteins[i]:
            accession = membrane_proteins[i].split('-')[1]
            membrane_accessions.append(accession)
        else:
            accession = membrane_proteins[i]
            membrane_accessions.append(accession)

    # for i in range(0, membrane_accessions):
    #     if membrane_accessions[i] in genome_missing_proteins:
    #         num_membrane_proteins -= 1
    genome_membrane_accessions = []
    for protein in membrane_accessions:
        if protein not in genome_missing_proteins:
            genome_membrane_accessions.append(protein)

    num_membrane_proteins = len(genome_membrane_accessions)            

    if num_membrane_proteins == 0:
        return 'red'
    
    for protein in genome_membrane_accessions:
        key = tcid + '-' + protein
        curr_row = pd.DataFrame([multi_dict[key]])
        if tm_overlap_decision(curr_row) == True:
            return 'yellow'
    
    return 'red'

    # protein_tmoverlap_dict = {}
    # # collect the tms overlap score for all the proteins 
    # # PRECONDITION: all proteins here are in the genome
    # for protein in membrane_proteins:
    #     curr_row = pd.DataFrame([multi_dict[protein]])
    #     protein_tmoverlap_dict[protein] = {'overlap_percent': overlap_percent(curr_row)}


    # return 'red'

def isMultiComp(row,df,input):
    if isSingleComp(row):
        return (False)
    ##get rid of single comp system first
    tcid = row["Hit_tcid"]
    Fusion_Add=[]
    Fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]
    sortedGeneArr = []
    if(len(Fusion_results) !=1):
        sortedGeneArr = sorted(Fusion_results, key=lambda x: x['sstart'])
        Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]             
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))

    fusions = ''
    id = row['#Query_id']
    if id == 'YP_500900.1':
        print('here')
    tcdb_tms = row['Hit_n_TMS']
    if(set(tc_all_arr)==set(tc_filter_arr)):##If all the proteins in that system can be found, then green
        #print(tc_arr)
        #print("green",tc_filter_arr,"Fusion Results:",Fusion_Add)
       
        if(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True):
            return({"color":"Green",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0})

    MembraneProteins= FindMembraneProtein(row, df)
    decision = multicomp_decision(tc_arr, tc_missing_arr, MembraneProteins)
    # if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
    if(input*len(tc_all_arr)<=len(tc_filter_arr)) or len(set(MembraneProteins) & set(tc_filter_arr))>0:
        ##given some proteins can be found while containing the membrane proteins 
       if final_decision(eVal(row), qCoverage(row), hCoverage(row), pfamDoms(row), sortedGeneArr, overlap_percent(row), tcdb_tms) == 'green' or final_decision(eVal(row), qCoverage(row), hCoverage(row), pfamDoms(row), sortedGeneArr, overlap_percent(row), tcdb_tms) == 'yellow':
            if decision == 'yellow' or decision == 'No membrane proteins':
                return({"color":"Yellow",
                        "Found_proteins":tc_filter_arr,
                        "All_proteins":tc_arr,
                        'Missing_proteins':tc_missing_arr,
                        "Fusion_results":Fusion_Add,
                        "isFusion":len(Fusion_Add)>0})


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
    fusions = ''
    for fusion in Output_dict['Fusion_results']:
        fusions += fusion + ' '
    Intermediate=Output_df_row.copy()
    Intermediate['isFusion'] = fusions
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
        Missing_infor = hmmtop_df.loc[(hmmtop_df['Hit_xid'] ==  hit_xid)]   # issue here

        _Intermediate["Match_length"]=Missing_infor["Match_length"].iloc[0] if not Missing_infor.empty else "NA"
        _Intermediate["Hit_n_TMS"]=Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
        with open(filename, mode="a", encoding='utf-8') as f:
            _Intermediate.to_csv(f, sep='\t', header=False, index=False)
    
        
        

def write_singlecomp(output_dict,Output_df_row):
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
    single = isSingleComp(row)

    # if row['Hit_tcid'] == '1.C.3.4.2':
    #     print('here')
    if(not single):
        Output_dict = isMultiComp(row, df, 0.5)
        Write_multicomp(Output_dict,Output_df.loc[[index],Output_df.columns])
    else:
        output_dict = categorizeSingleComp(row)
        write_singlecomp(output_dict, Output_df.loc[[index],Output_df.columns])


    

    #print(Output_dict)