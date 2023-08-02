#!/usr/bin/env python

import argparse
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import csv
import numpy as np
import os
import re
import pickle as pic
import subprocess
# from fusion_distribution import isFusion,geneFusions
# from fusion_dist_dir import isFusion,genDict
from fusion import isFusion, geneFusions, genDict, setGenome, check_overlap, fusion_TMS_count
from parseXML import parse
from pprint import pprint

GENOME = ''
TMS_THRESH_GREEN = 0.8
TMS_THRESH_YELLOW = 0.5
EXCELLENT_EVAL = 1e-30
E_VAL_GREEN = 1e-10
E_VAL_YELLOW = 1e-3 
Q_COV_THRESH = 80
S_COV_THRESH = 80
LOW_COV = 50
AUTO_RED = 20
Membraneprotein_threshold=3
Hit_TMS_Diff= 2
MIN_PROTEINS = 0.5

tcdbSystems = {}
green = {}
yellow = {}
red = {}
hmmtop_df = pd.DataFrame(columns=['Hit_tcid', 'Hit_xid', 'Hit_n_TMS','Match_length'])
geneFusions={}
df = None
id_dict = {}
tcid_assignments = {}

def create_id_dict(df):
    if df == None:
        return {}
    id_df = df[['Hit_tcid', 'Hit_xid', '#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
     'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
     ,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]
    
    id_dict = {}
    for index, row in id_df.iterrows():
        if row['Hit_tcid'] in id_dict:
            id_dict[row['Hit_tcid']].append(row)
        else:
            id_dict[row['Hit_tcid']] = [row]
    return id_dict


def adjustOverlapScore(df):
    print(GENOME)
    overlap_dict = parse(GENOME + '/hmmtop.db', GENOME + '/xml/' ,GENOME + '/results.tsv', 8)
    score_dict = {}
    for k in overlap_dict:
        score_dict[k] = overlap_dict[k]['alignedTMS']
    for index, row in df.iterrows():
        if row['#Query_id'] in score_dict:
            df.at[index, 'TM_Overlap_Score'] = score_dict[row['#Query_id']]
        else:
            df.at[index, 'TM_Overlap_Score'] = 0
    print('finished')

#This function will fill the dictionary tcdbSystems with data.
def parseTCDBcontent(input):
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
    # print(fusions)
    if len(fusions) == 0:
        return 'red'
    tot_cov = int(fusions[0]['send']) - int(fusions[0]['sstart'])
    for i in range(1, len(fusions)):
        # caluculate overlap
        overlap = 0
        if fusions[i]['sstart'] < fusions[i-1]['send']:
            overlap = fusions[i-1]['send'] - fusions[i]['sstart']
        tot_cov += fusions[i]['send'] - fusions[i]['sstart'] - overlap
        if fusions[i]['send'] < fusions[i-1]['send']:
            return 'red'
        if fusions[i]['sstart'] == fusions[i-1]['sstart']:
            return 'red'
        '''
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
        '''
    if tcdb_tms == 0:
        overlap_percent_fus = 1
    else:
        net_overlap = fusion_TMS_count(fusions)
        overlap_percent_fus = net_overlap / tcdb_tms
        

    tot_cov = float(tot_cov / int(fusions[0]['hit_length'])) * 100
    
    if tot_cov >= S_COV_THRESH and overlap_percent_fus >= TMS_THRESH_GREEN:
        return 'green'
    elif tot_cov >= LOW_COV and overlap_percent_fus >= TMS_THRESH_YELLOW:
        return 'yellow'
    else:
        return 'red'
    

def make_decision(row, fusions, tcdb_tms):
    eval = eVal(row)
    qcov = qCoverage(row)
    scov = hCoverage(row)
    pfam_doms = pfamDoms(row)
    tmoverlap = overlap_percent(row)
    # automatic green if e-val is 0 or extremely good with okay coverage and has some tms coverage
    if (eval == 0.0 or eval <= EXCELLENT_EVAL) and (qcov > LOW_COV and scov > LOW_COV) and tcdb_tms > 0:
        if tm_overlap_decision(row) == True or len(pfam_doms) > 0:
            return 'green'
        else:
            return 'yellow'
    elif (eval == 0.0 or eval <= EXCELLENT_EVAL) and (qcov > LOW_COV and scov > LOW_COV) and len(pfam_doms) > 0:
        return 'green'
    
    # if coverage is less than a very low number it is impossible for it to be a good hit
    if qcov <= AUTO_RED or scov <= AUTO_RED and len(fusions) == 0:
        return 'red'
    # if there are less than 3 tms, there is a possiblity of mischaracterizing tms regions
    if tcdb_tms > 3:
        # if cov is too low it is automatically red
        if qcov < LOW_COV or scov < LOW_COV and len(fusions) == 0:
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
        # okay e value, has good coverage and there is a high tms overlap means good hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and tm_overlap_decision(row) == True):
            return 'green'
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_GREEN and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_GREEN):
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        # ok e-val, no common doms, has good coverage or there is a good tms overlap means yellow hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) or tmoverlap >= TMS_THRESH_YELLOW):
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
        # great e-val no matching pfam domains but has goood coverage or good tmoverlap its a good hit
        if eval <= E_VAL_YELLOW and ((qcov >= Q_COV_THRESH and scov >= S_COV_THRESH) and tm_overlap_decision(row) == True):
            return 'green'
        # if there is a fusion found we need to do further analysis
        if len(fusions) > 0:
            return 'fusion found'
        # has ok e value, common doms match, good cov is good hit
        if eval <= E_VAL_YELLOW and len(pfam_doms) > 0 and (qcov >= LOW_COV and scov >= LOW_COV):
            return 'yellow'
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

def final_decision(row, fusions, tcdb_tms):
    eval = eVal(row)
    qcov = qCoverage(row)
    scov = hCoverage(row)
    pfam_doms = pfamDoms(row)
    tmoverlap = overlap_percent(row)
    decision = make_decision(row, fusions, tcdb_tms)
    if decision == 'fusion found':
        # if assign_if_fusion returns red its a false fusion or unlikely
        if assign_if_fusion(fusions, tcdb_tms) == 'red':
            return make_decision(eval, qcov, scov, pfam_doms, [], tmoverlap, tcdb_tms)
        else:
            return assign_if_fusion(fusions, tcdb_tms)
    else:
        return decision


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
    tcid = row['Hit_tcid']
    tcdb_tms = row['Hit_n_TMS']    
    fus_list = []
    fusion_results= geneFusions[row["Hit_tcid"] + "-" + row["Hit_xid"]]

    if(len(fusion_results) !=1):
        sortedGeneArr = sorted(fusion_results, key=lambda x: x['sstart'])
        if(len(isFusion(sortedGeneArr))!=0):
            sortedGeneArr = isFusion(sortedGeneArr)
            for k in range(0, len(sortedGeneArr)):
                fusions += sortedGeneArr[k]['query'] + ' '
            row_to_write.append(fusions)
        else:
            sortedGeneArr = []
            row_to_write.append(fusions)
    else:
        row_to_write.append(fusions)
    if assign_if_fusion(sortedGeneArr, tcdb_tms) == 'red':
        sortedGeneArr = []
    
    if tcid == '2.A.6.6.11':
        print(tm_overlap_decision(row))
        #print(tcdb_tms)
    if len(fusions) > 0:
        fus_list = fusions.split(' ')
        if row['#Query_id'] not in fus_list:
            row_to_write.pop()
            row_to_write.append('')
    fus_color = ''
    if len(fusions) != 0:
        fus_color = assign_if_fusion(sortedGeneArr, tcdb_tms)
    if tcid in tcid_assignments:
        return (tcid_assignments[tcid], row_to_write)
    
    decision = final_decision(row, sortedGeneArr, tcdb_tms)
    if tcid == '1.B.12.3.3':
        print(decision)
    if decision == 'green':
        green[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Green'
        return ('Green', row_to_write)
    elif decision == 'yellow':
        yellow[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Yellow'
        return ('Yellow', row_to_write)
    elif decision == 'red':
        red[row.get(key='#Query_id')] = row_to_write
        tcid_assignments[row.get(key='Hit_tcid')] = 'Red'
        return ('Red', row_to_write)

    

    
def execute_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    process.wait()
    if process.returncode == 0:
        print(f"command'{command}'success")
    else:
        print(f"command'{command}'fail")

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
    Found_arr_gblast = find_df[(find_df['Hit_n_TMS'] >= Membraneprotein_threshold) & (abs(find_df['Hit_n_TMS'] - find_df['Query_n_TMS']) <= Hit_TMS_Diff)]["Hit_xid"].tolist()
    ##If there are proteins that has Hit TMS >=3 and difference with its query TMS <=2, we say it's membrane proteins
    tcdb_arr= Found_arr_gblast +tc_Notfind_arr##This arr contains the possible output proteins 
    if len(tcdb_arr)==0:
        return([])
    tcdb_df = pd.concat([hmmtop_df.query(f"Hit_xid=='{arr}' and Hit_tcid=='{tcid}'") for arr in tcdb_arr])
    Final_return = tcdb_df[(tcdb_df['Hit_n_TMS']>= Membraneprotein_threshold)]  
    ##print(tc_all_arr)
    ##print(find_df)
    ##print(unfind_df)
    ##print(Final_return["Hit_xid"].tolist())
    # return((Final_return['Hit_tcid'] + '-' + Final_return["Hit_xid"]).tolist())  

    return Final_return.to_dict()
    ##final_df = pd.concat([df.query(f"Hit_xid=='{arr}'") for arr in tc_filter_arr])

def isEmpty(dict):
    if type(dict) is list:
        return len(dict) == 0
    count = 0
    for key in dict:
        if not bool(dict[key]):
            count += 1
    
    return count == 4

def multicomp_decision(tcdb_proteins, protein_type, mem_dict, tc_filter_arr):
    max_tms = 0
    tcid = tcdb_proteins[0].split('-')[0]
    if (isEmpty(mem_dict) or len(mem_dict) == 0) and (len(tc_filter_arr)/len(tcdb_proteins)) < MIN_PROTEINS:
       return 'red'
    elif isEmpty(mem_dict) or len(mem_dict) == 0:
       return 'yellow'
    
    
    membrane_protein_tms = {}

    membrane_proteins = []
    if type(mem_dict) is list:
        membrane_proteins = mem_dict
    else:
        membrane_proteins_accessions = list(mem_dict["Hit_xid"].values())
        membrane_proteins_tcids = list(mem_dict["Hit_tcid"].values())
        membrane_tmss = list(mem_dict["Hit_n_TMS"].values())
        
        membrane_proteins = membrane_proteins_accessions
        for i in range(0, len(membrane_proteins_accessions)):
            protein = membrane_proteins_tcids[i] + '-' + membrane_proteins_accessions[i]
            membrane_protein_tms[protein] = membrane_tmss[i]
            membrane_tmss = list(mem_dict['Hit_n_TMS'].values()) 
        # sorts the dictionary based on ascending tms values
        membrane_protein_tms = {k: membrane_protein_tms[k] for k in sorted(membrane_protein_tms, key=membrane_protein_tms.get, reverse=True)}
        for key in membrane_protein_tms:
            max_tms = membrane_protein_tms[key]
            break
    
     
    # if beta barrell and most proteins are in the system should be yellow->green
    if protein_type == '1.B' and abs(len(tcdb_proteins) - len(tc_filter_arr)) <= 2:
        return 'yellow->green'
    all_mems_found = True
    if type(mem_dict) is list:
        for protein in membrane_proteins:
            if protein not in tc_filter_arr:
                all_mems_found = False
                break
    else:
        for protein in membrane_proteins:
            accession = protein
            if accession not in tc_filter_arr:
                all_mems_found = False
                break
    if all_mems_found == True:
        return 'yellow->green'
    
    
    if len(tcdb_proteins) > 3 and (len(tc_filter_arr)/len(tcdb_proteins)) < MIN_PROTEINS:
        return 'red'
    # membrane_accessions = []
    # num_membrane_proteins = len(membrane_proteins)
    # for protein in genome_missing_proteins:
    #     key = tcid + '-' + protein
    #     if key in membrane_protein_tms:[GENOME + "Green.tsv",GENOME + "Red.tsv",GENOME + "Yellow.tsv"]   
    #         membrane_protein_tms.pop(key)
    #         num_membrane_proteins -= 1

    # for i in range(0, num_membrane_proteins):
    #     if '-' in membrane_proteins[i]:
    #         accession = membrane_proteins[i].split('-')[1]
    #         membrane_accessions.append(accession)
    #         membrane_protein_tms.pop(membrane_proteins[i])
    #     else:
    #         accession = membrane_proteins[i]
    #         membrane_accessions.append(accession)

    # for i in range(0, membrane_accessions):
    #     if membrane_accessions[i] in genome_missing_proteins:
    #         num_membrane_proteins -= 1
    genome_membrane_accessions = tc_filter_arr

    num_membrane_proteins = len(genome_membrane_accessions)            

    if num_membrane_proteins == 0:
        return 'red'
    
    for protein in genome_membrane_accessions:
        key = tcid + '-' + protein
        if key in membrane_protein_tms:
            if max_tms - membrane_protein_tms[key] < 2:
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
        if len(isFusion(sortedGeneArr)) != 0:
            Fusion_Add=[x["query"] for x in isFusion(sortedGeneArr)]
        else:
            sortedGeneArr = []
    tc_arr = tcdbSystems.get(tcid)
    tc_all_arr= [arr.split("-")[1] for arr in tc_arr] ##All the proteins that TCDB provide
    tc_filter_arr= list(filter(lambda x: (len(df.query(f"Hit_xid=='{x}' and Hit_tcid=='{tcid}'"))) !=0, tc_all_arr))
    tc_missing_arr= list(set(tc_all_arr) - set(tc_filter_arr))

    fusions = '' 
    qid = row['#Query_id']
    
    if len(Fusion_Add) > 0:
        if qid not in Fusion_Add:
            Fusion_Add = []       

    protein_type = row['Hit_tcid'].split('.')[0] + '.' + row['Hit_tcid'].split('.')[1]
    tcdb_tms = row['Hit_n_TMS']
    if tcid == '2.A.2.3.3':
        if tcid in tcid_assignments:
            print(tcid_assignments[row['Hit_tcid']])
        else:
            print('not assigned yet')


    if row['Hit_tcid'] in tcid_assignments:
        if tcid_assignments[tcid] == 'yellow->green':
            return({"color": 'Green',
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0,
                    "Initial_decision": 'Yellow'})
        return({"color": tcid_assignments[row['Hit_tcid']],
                "Found_proteins":tc_filter_arr,
                "All_proteins":tc_arr,
                'Missing_proteins':tc_missing_arr,
                "Fusion_results":Fusion_Add,
                "isFusion":len(Fusion_Add)>0,
                "Initial_decision": tcid_assignments[row['Hit_tcid']]})
    
    if(set(tc_all_arr)==set(tc_filter_arr)):##If all the proteins in that system can be found, then green
        #print(tc_arr)
        #print("green",tc_filter_arr,"Fusion Results:",Fusion_Add)
        #if tcid == '3.A.1.12.4':
        #    print(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True)
        #if(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True):
        if(eVal(row) <= E_VAL_GREEN and qCoverage(row) >= Q_COV_THRESH and hCoverage(row) >= S_COV_THRESH and len(pfamDoms(row)) != 0 and tm_overlap_decision(row) == True):    
            tcid_assignments[row['Hit_tcid']] = 'Green'
            return({"color":"Green",
                    "Found_proteins":tc_filter_arr,
                    "All_proteins":tc_arr,
                    'Missing_proteins':tc_missing_arr,
                    "Fusion_results":Fusion_Add,
                    "isFusion":len(Fusion_Add)>0,
                    "Initial_decision": 'green'})

    MembraneProteins= FindMembraneProtein(row, df)
    count_mem = 0
    total_tcmem = 0
    if len(MembraneProteins) > 0 or not isEmpty(MembraneProteins):
        if type(MembraneProteins) is not list:
            mem_proteins = list(MembraneProteins['Hit_xid'].values())
            total_tcmem = len(mem_proteins)
        else:
            mem_proteins = MembraneProteins
            total_tcmem = len(MembraneProteins)

        for p in mem_proteins:
            if p in tc_filter_arr:
                count_mem += 1
    # temp = isEmpty(MembraneProteins)
    if row['Hit_tcid'] == '2.A.63.1.1':
        print(MembraneProteins)
    decision = multicomp_decision(tc_arr, protein_type, MembraneProteins, tc_filter_arr)
    if row['Hit_tcid'] == '2.A.63.1.1':
        print(decision)
        print(len(set(tc_filter_arr)))
<<<<<<< Updated upstream
    # if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
    if (total_tcmem > 0 and input <= float(count_mem / total_tcmem)) or (len(set(tc_filter_arr))>0 and isEmpty(MembraneProteins)):
        ##given some proteins can be found while containing the membrane proteins 
=======
    if len(tc_missing_arr) == 0 and final_decision(row, sortedGeneArr, tcdb_tms) == 'yellow':
        tcid_assignments[row['Hit_tcid']] = 'Yellow'
        return({"color":"Yellow",
                "Found_proteins":tc_filter_arr,
                "All_proteins":tc_arr,
                'Missing_proteins':tc_missing_arr,
                "Fusion_results":Fusion_Add,
                "isFusion":len(Fusion_Add)>0,
                'Initial_decision': 'Yellow'})
    # if(input*len(tc_all_arr)<=len(tc_filter_arr)) and len(set(MembraneProteins) & set(tc_filter_arr))>0:
    if (total_tcmem > 0 and input <= float(count_mem / total_tcmem)) or (len(set(tc_filter_arr))>0 and isEmpty(MembraneProteins)):
        ##given some proteins can be found while containing the membrane proteins
        if final_decision(row, sortedGeneArr, tcdb_tms) == 'yellow' or final_decision(row, sortedGeneArr, tcdb_tms) == 'green':
>>>>>>> Stashed changes
            if decision == 'yellow':
                tcid_assignments[row['Hit_tcid']] = 'Yellow'
                return({"color":"Yellow",
                        "Found_proteins":tc_filter_arr,
                        "All_proteins":tc_arr,
                        'Missing_proteins':tc_missing_arr,
                        "Fusion_results":Fusion_Add,
                        "isFusion":len(Fusion_Add)>0, 
                        'Initial_decision': 'Yellow'})
            elif decision == 'yellow->green':
                tcid_assignments[tcid] = 'yellow->green'
                return({"color":"Green",
                        "Found_proteins":tc_filter_arr,
                        "All_proteins":tc_arr,
                        'Missing_proteins':tc_missing_arr,
                        "Fusion_results":Fusion_Add,
                        "isFusion":len(Fusion_Add)>0,
                        'Initial_decision': 'Yellow'})

    tcid_assignments[row['Hit_tcid']] = 'Red'
    return({"color":"Red",
            "Found_proteins":tc_filter_arr,
            "All_proteins":tc_arr,
            'Missing_proteins':tc_missing_arr,
            "Fusion_results":Fusion_Add,
            "isFusion":len(Fusion_Add)>0,
            'Initial_decision': 'Red'})
   


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
    # if len(Output_dict['Missing_proteins']) != 0:
    #     print('here')
    fusions = ''
    for fusion in Output_dict['Fusion_results']:
        fusions += fusion + ' '
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = fusions
    if len(Output_dict['Missing_proteins']) == 0:
        Intermediate['Missing_components'] = 'NA'
    else:
        Intermediate['Missing_components'] = str(Output_dict['Missing_proteins'])
    Intermediate['Initial_decision'] = Output_dict['Initial_decision']
    filename=f"{Output_dict['color']}.tsv"
    filename = GENOME + '/analysis/' + filename
    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)
    
    for hit_xid in Output_dict['Missing_proteins']:
        _Intermediate = Output_df_row.copy()
        _Intermediate = _Intermediate.applymap(lambda x: 'NA')

        _Intermediate["Hit_tcid"] = Output_df_row["Hit_tcid"]
        _Intermediate["Hit_xid"] = hit_xid
        Missing_infor = hmmtop_df.loc[hmmtop_df['Hit_xid'] == hit_xid]

        _Intermediate["Match_length"] = Missing_infor["Match_length"].iloc[0] if not Missing_infor.empty else "NA"
        _Intermediate["Hit_n_TMS"] = Missing_infor["Hit_n_TMS"].iloc[0] if not Missing_infor.empty else "NA"
        with open(filename, mode="a", encoding='utf-8') as f:
            _Intermediate.to_csv(f, sep='\t', header=False, index=False)
    
        
        

def write_singlecomp(output_dict,Output_df_row):
    color = output_dict[0]
    dictionary = output_dict[1]
    Intermediate = Output_df_row.copy()
    Intermediate['isFusion'] = dictionary[len(dictionary) - 1]
    missing_proteins = "NA"
    Intermediate['Missing_components']= missing_proteins
    Intermediate['Initial_decision'] = 'NA'
    filename=f"{color}.tsv"
    filename = GENOME + '/analysis/' + filename

    filemode='a' if os.path.exists(filename) else 'w' 
    #print(Intermediate)
    with open(filename, mode=filemode, encoding='utf-8') as f:
        Intermediate.to_csv(f, sep='\t', header=filemode=='w', index=False)
    




def remove_duplicates(input_file, output_file):
    seen = set()  # Track seen rows
    with open(input_file, 'r', newline='') as input_file, \
         open(output_file, 'w', newline='') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')

        for index, row in enumerate(reader):
            if index == 0:  # Write the first row
                writer.writerow(row)
                continue

            row_tuple = tuple(row)  # Convert row to a hashable tuple

            if row_tuple not in seen:
                seen.add(row_tuple)
                writer.writerow(row)

def custom_sort(row):
    parts = row.split('.')
    num1 = int(parts[0])
    letter = parts[1]
    num2 = int(parts[2])
    num3 = int(parts[3])
    num4 = int(parts[4])
    return num1, letter, num2, num3, num4

        

    #print(Output_dict)

def main():
#def run():
    # GENOME = 'MicrobiomeResults/GCF_000013425.1'  # we can make this a command line argument
    
    parser = argparse.ArgumentParser(description='Microbiome analysis script')
    parser.add_argument('-g', '--genome', type=str, help='Genome ID')
    parser.add_argument('-q', '--qcov', type=float, help='Query coverage threshold')
    parser.add_argument('-s', '--scov', type=float, help='Subject coverage threshold')
    parser.add_argument('-r', '--autored', type=float, help='Lowest possible coverage allowed')
    parser.add_argument('-m', '--membrane', type=int, help='Membrane proteins threshold')
    parser.add_argument('-t', '--tmsdiff', type=int, help='Lowest amount of TMSs to be a Membrane Protein')
    args = parser.parse_args()
    '''
    if not os.path.exists(args.genome):
        print(f"Genome folder not found: {args.genome}")
        exit()
    '''
    global GENOME
    global Q_COV_THRESH
    global S_COV_THRESH
    global AUTO_RED
    global Membraneprotein_threshold
    global Hit_TMS_Diff
    global df

    if args.genome is not None:
        GENOME = args.genome
    elif len(GENOME) == 0:
        return 'Genome not provided'
    if args.qcov is not None:
        Q_COV_THRESH = args.qcov
    if args.scov is not None:
        S_COV_THRESH = args.scov
    if args.autored is not None:
        AUTO_RED = args.autored
    if args.membrane is not None:
        Membraneprotein_threshold = args.membrane
    if args.tmsdiff is not None:
        Hit_TMS_Diff = args.tmsdiff
    print('starting')
    # construct df with all data from results.tsv
    df = pd.read_table(GENOME + '/results.tsv')
    setGenome(GENOME)
    # print columns of df for dev use in constructing filtered_df later
    print(GENOME)
    if not os.path.exists(GENOME + '/analysis/'):
        os.mkdir(GENOME + '/analysis')

    # create empty dfs for tagging with green, yellow, and red labels
    green_df = df.copy()
    green_df = green_df.iloc[0:0]
    yellow_df = green_df.copy()
    red_df = green_df.copy()

    # construct filtered_df with only relevant data points for each protein matche
    filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
    'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
    print('adjusting overlap score')
    adjustOverlapScore(df)

    Output_df= df[['Hit_tcid','Hit_xid','#Query_id','Match_length','e-value','%_identity','Query_Length','Hit_Length','Q_start',
    'Q_end','S_start','S_end','Query_Coverage','Hit_Coverage','Query_n_TMS','Hit_n_TMS','TM_Overlap_Score','Family_Abrv'
    ,'Predicted_Substrate','Query_Pfam','Subject_Pfam']]

    df.to_csv('adj.csv', index=False)

    # To distinguish between single and multiple component systems. 
    # Example content:
    #  1.A.1.1.1 =>  ["1.A.1.1.1-P0A334"],
    #  3.A.1.1.1 =>  ["3.A.1.1.1-P02916", ""3.A.1.1.1-XXXXX", "3.A.1.1.1-YYYYYY", "3.A.1.1.1-ZZZZZZZZZ"],
    #  ....

    input = open("tcdb.faa")


    parseTCDBcontent(input)

    ##with open('hmmtop.out') as f:
    hmmtop_path=GENOME+ "/analysis/hmmtop.out"
    ##Test db file
    '''
    with open(hmmtop_path_db, "rb") as file:
            
            db_data = pic.load(file)
            print(type(db_data))
            csv_file_path = "hmmtop_test_db.csv"  
            pd.DataFrame.from_dict(db_data).to_csv(csv_file_path, index=False)
    '''
    print(hmmtop_path)
    print(os.path.exists(hmmtop_path))        
    if not os.path.exists(hmmtop_path):
        TCDB_seqs = GENOME + "/analysis/mcs_tcids"
        if os.path.exists(TCDB_seqs):
            os.system(f"rm -r {TCDB_seqs}")
        os.mkdir(TCDB_seqs)
        data = pd.read_csv(GENOME + '/results.tsv', sep='\t')
        hit_tcid_array = data['Hit_tcid'].unique()
        try:
            with open(GENOME + '/analysis/mcs_tcids.txt', 'w') as f:
                for tcid in hit_tcid_array:
                    f.write(tcid + '\n')
                f.close()
        except Exception as e:
            print("An error occurred:", str(e))
        command_1=[f"extractTCDB.pl -i {GENOME}/analysis/mcs_tcids.txt -o {GENOME}/analysis/mcs_tcids -f fasta"]
        print(command_1)
        joined_commands = ';'.join(command_1)
        execute_command(joined_commands)
        command2=f"cd {GENOME};cat analysis/mcs_tcids/*faa > analysis/all.faa &&hmmtop -if=analysis/all.faa -of=analysis/hmmtop.out -sf=FAS -pi=spred -is=pseudo"
        print(command2)
        execute_command(command2)
        
        
    with open(hmmtop_path) as f:
        lines = f.readlines()
        ##Create a dataframe that include above columns 
        for line in lines:
            fields = re.split(r'[ -]+', line.strip())
            ##split the system and protein names
            new_col=[fields[2],fields[3],fields[5]]
            new_row = pd.Series([fields[2],fields[3],fields[5],fields[1]], index=hmmtop_df.columns)
            hmmtop_df.loc[len(hmmtop_df)] = new_row
        hmmtop_df['Hit_n_TMS'] = hmmtop_df['Hit_n_TMS'].astype(int)
    
    parseTCDBcontent(input)
    Missing_protein_list=[]
    genDict(geneFusions, GENOME)
    for filename in [GENOME + "/analysis/Green.tsv",GENOME + "/analysis/Red.tsv",GENOME + "/analysis/Yellow.tsv"]:
        if os.path.exists(filename):
            os.remove(filename)
            

    for index, row in df.iterrows():
        single = isSingleComp(row)

        # if row['Hit_tcid'] == '1.C.3.4.2':
        #     print('here')
        if(not single):
            Output_dict = isMultiComp(row, df, 0.5)
            # tcid_assignments[row['Hit_tcid']] = Output_dict['color']
            Write_multicomp(Output_dict,Output_df.loc[[index],Output_df.columns])
        else:
            output_dict = categorizeSingleComp(row)
            if row['Hit_tcid'] == '2.A.2.3.3':
                print(output_dict)
            write_singlecomp(output_dict, Output_df.loc[[index],Output_df.columns])

    input_files = ['Green.tsv', 'Red.tsv', 'Yellow.tsv']
    output_files = ['Green_adj.tsv', 'Red_adj.tsv', 'Yellow_adj.tsv']

    for i in range(len(input_files)):
        in_file = GENOME + '/analysis/' + input_files[i]
        out = GENOME + '/analysis/' + output_files[i]
        remove_duplicates(in_file, out)


    print([GENOME + "Green.tsv",GENOME + "Red.tsv",GENOME + "Yellow.tsv"])
    for filename in [GENOME + "/analysis/Green.tsv",GENOME + "/analysis/Red.tsv",GENOME + "/analysis/Yellow.tsv"]: 
        if os.path.exists(filename):
            df = pd.read_csv(filename, sep='\t')
            df = df.fillna('NA')
            NA_count = df.eq('NA').sum(axis=1).rename('NA_count')
            df['na_sort'] = NA_count
            df = df.sort_values(by=['Hit_tcid','na_sort']).drop(columns=['na_sort'])
            mask = (df['#Query_id'] == 'NA') & (df['e-value'] == 'NA')
            df_temp = df.loc[mask] 
            df_temp.drop_duplicates(subset=['Hit_tcid','Hit_xid'], keep='first', inplace=True) 
            df.loc[mask] = df_temp
            df.dropna(inplace=True)
            df.to_csv(filename, sep='\t', index=False)
            # df.sort_values(by=df.columns[0], key=lambda x: x.map(custom_sort))
            command = f"sort -o {filename} -t '.' -k1,1n -k2,2 -k3,3n -k4,4n -k5,5n {filename}"
            subprocess.run(command, shell=True)

main()
