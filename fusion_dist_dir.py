import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import os

import seaborn as sns
import seaborn.objects as so

'''
WP_024544048.1(c).png
'''

def fusion_dist(fusion_candidate_list):

    plt.cla()
    y_count = len(fusion_candidate_list) + 0.5

    plt.plot(0,0)
    plt.plot([0,fusion_candidate_list[0]['hit_length']],[y_count+1,y_count+1],label = "TCDB Reference" + "(0-" + str(fusion_candidate_list[0]['hit_length']) + ")")

    coverage_list = []
    for candidate in fusion_candidate_list:
        coverage_list = coverage_list + (list(range(candidate['sstart'], candidate['send'])))
        plt.plot([candidate['sstart'],candidate['send']], [y_count, y_count], label = candidate['query'] + "(" + str(candidate['sstart']) + "-" + str(candidate['send']) + ")")
        y_count -= 1

    residue_cov = 0
    frequency_list = {x:coverage_list.count(x) for x in coverage_list}
    for key in frequency_list:
        if(frequency_list[key] > 0):
            residue_cov += 1
    
    if(residue_cov/fusion_candidate_list[0]['hit_length'] > 0.8):
        # print(fusion_candidate_list[0]['query'])
        # print("COV LIST")
        # print(coverage_list)
        # print("-------")
        # print("FREQ LIST")
        # print(frequency_list)
        # print("-------")
        #x= input("next")
        plt.plot(frequency_list.keys(), frequency_list.values(), '.--',label = "Frequency Distribution")
        #for point in frequency_list:
            #plt.plot(point, frequency_list[point])
        plt.savefig("output/" + fusion_candidate_list[0]['query'] + "(c).png", dpi=300)
        plt.legend()
        plt.savefig("output/" + fusion_candidate_list[0]['query'] + "(l).png", dpi=300)


# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_FUSION_COV_OVERLAP = 20
MIN_QUERY_COV_THRESHOLD = 80
MIN_SUBJECT_COV_THRESHOLD = 5
MAX_SUBJECT_COV_THRESHOLD = 60 # keep at 60
MIN_FUSION_COVERAGE = 80

# generates a dictionary from the dataframe
def genDict(dict, input_df):
    for index, row in input_df.iterrows():
        new_entry = {"query": row["#Query_id"], "qcov": row["Query_Coverage"], "sstart": row["S_start"], "send": row["S_end"], "scov": row["Hit_Coverage"], "hit_length": row["Hit_Length"]}
        tcid = row["Hit_tcid"] + "-" + row["Hit_xid"]
        if tcid in dict:
            dict.get(tcid).append(new_entry)
            #print("Appending : " + tcid)
        else:
            dict[tcid] = [new_entry]
            #print("New entry")

# takes in sorted array (subset of dictionary) and 
# outputs list of dictionary (each dictionary is a fusion candidate)
def isFusion(sortedArr):

    #print("-------")
    #print("Input is " + str(sortedArr))

    # output list
    fus_list = []
    overlap_length = 0
    
    # iteratate through each fusion candidate 
    for i in range(len(sortedArr)):
        if sortedArr[i]['scov'] > MAX_SUBJECT_COV_THRESHOLD:
            #print("i rejected, candidate scov too long")
            continue
        elif sortedArr[i]['scov'] < MIN_SUBJECT_COV_THRESHOLD:
            #print("i rejected, candidate scov was too small")
            continue
        elif sortedArr[i]['qcov'] < MIN_QUERY_COV_THRESHOLD:
            #print("i rejected, candidate qcov was too small")
            continue
        else:
            #print("Accepting i protein: " + sortedArr[i]["query"])
            fus_list.append(sortedArr[i])
    tot_length = 0
    overlap_length = 0

    #checks for extra fusions if proteins are lined up end to end
    for i in range(len(fus_list)):
        tot_length += (fus_list[i]['send'] - fus_list[i]['sstart'])
        if i < len(fus_list)-1:
            overlap_length += min(fus_list[i+1]['send'], fus_list[i]['send']) - max(fus_list[i+1]['sstart'], fus_list[i]['sstart'])
    tot_length -= overlap_length
    
    if len(fus_list) == 0:
        #print("No candidates identified")
        return []
    elif len(fus_list) < 2:
        #print("Insufficient Candidates identified")
        #print(fus_list)
        return []
    elif (tot_length/fus_list[0]['hit_length'])*100 < MIN_FUSION_COVERAGE:
        #print("Didn't meet MIN_THRESHOLD")
        #print("total length: " + str(tot_length))
        #print("hit length: " + str(fus_list[0]['hit_length']))
        #print("-------")
        return []
    else:
        #print("output is: " + str(fus_list))
        #print("total length: " + str(tot_length))
        #print("hit length: " + str(fus_list[0]['hit_length']))
        return fus_list


def fusion_analysis(input_tsv):
    df = pd.read_table(input_tsv)

    geneFusions = {}

    genDict(geneFusions, df)

    #sorts the genomeFusions dictionary by the genome start index
    num_out = 0
    num_in = 0
    for id in geneFusions:
        if(len(geneFusions[id]) == 1):
            continue
        sortedGeneArr = sorted(geneFusions[id], key=lambda x: x['sstart'])
        # print(sortedGeneArr)
        # x = input()
        num_in+=1
        if len(isFusion(sortedGeneArr)) != 0:
            #print(isFusion(sortedGeneArr))
            fusion_dist(isFusion(sortedGeneArr))
            num_out+=1
        #else:
            #print("FAIL")
            #isGraph = input()
            #if(isGraph == "Y"):
                #fusion_dist(isFusion(sortedGeneArr))
        #y = input()
    #print("-------")
    #print(str(num_out) + " fusions found out of " + str(num_in))

directory = 'data'
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        fusion_analysis(f)