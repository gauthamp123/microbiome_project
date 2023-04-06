import pandas as pd
import numpy as np
import os

df = pd.read_table('results.tsv')

geneFusions = {}
test_fus =  [{'query': 'A', 'qcov': 69.6, 'sstart': 295, 'send': 713, 'scov': 58.2}, {'query': 'B', 'qcov': 65.4, 'sstart': 405, 'send': 709, 'scov': 42.3}, {'query': 'C', 'qcov': 41.4, 'sstart': 481, 'send': 717, 'scov': 32.9}]

# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_FUSION_COV_OVERLAP = 20
MIN_QUERY_COV_THRESHOLD = 80
MIN_SUBJECT_COV_THRESHOLD = 5
MAX_SUBJECT_COV_THRESHOLD = 90
MIN_FUSION_COVERAGE = 20

# generates a dictionary from the dataframe
def genDict(dict, input_df):
    for index, row in input_df.iterrows():
        new_entry = {"query": row["#Query_id"], "qcov": row["Query_Coverage"], "sstart": row["S_start"], "send": row["S_end"], "scov": row["Hit_Coverage"]}
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

    print("-------")
    print("Input is " + str(sortedArr))

    # output list
    fus_list = []
    overlap_length = 0
    
    # iteratate through each fusion candidate 
    for i in range(len(sortedArr)):
        if sortedArr[i]['scov'] > MAX_SUBJECT_COV_THRESHOLD:
            print("i rejected, candidate scov too long")
            continue
        elif sortedArr[i]['scov'] < MIN_SUBJECT_COV_THRESHOLD:
            print("i rejected, candidate scov was too small")
            continue
        elif sortedArr[i]['qcov'] < MIN_QUERY_COV_THRESHOLD:
            print("i rejected, candidate qcov was too small")
            continue
        else:
            print("Analyzing i protein: " + sortedArr[i]["query"])
            # invalid fusion counter
            inv_fus_count = 0
            for j in range(i+1, len(sortedArr)):
                print("Comparing with j protein " + sortedArr[j]["query"])
                # This looks for proteins that have too much overlap with the current protein
                overlap_size = sortedArr[i]['send'] - sortedArr[j]['sstart']
                # vars represent the percent overlap
                perc_protein1 = (overlap_size / (sortedArr[i]['send'] - sortedArr[i]['sstart'])) * 100
                perc_protein2 = (overlap_size / (sortedArr[j]['send'] - sortedArr[j]['sstart'])) * 100
                
                # making sure theres not total overlap between 2 proteins (actually, this might be fine? )
                # if total overlap, reassign the overlap 
                if sortedArr[i]['send'] >= sortedArr[j]['send']:
                    overlap_size = sortedArr[j]['send'] - sortedArr[j]['sstart']

                # checking to see if the minimum overlap among comparable proteins is less than the threshold
                if (min(perc_protein1, perc_protein2)) > MAX_FUSION_COV_OVERLAP:
                    print("Exceeds MAX_OVERLAP between i and j")
                    inv_fus_count += 1
                # if theres no overlap move on
                elif ((sortedArr[i]['send'] - sortedArr[j]['sstart']) / (sortedArr[i]['sstart'] - sortedArr[j]['send'])) * 100 == 0:
                    print("No overlap between i and j; ending future j comparisons")
                    break
                else:
                    print("overlap between i and j added to net overlap_size var")
                    overlap_length += overlap_size
                
            # if there is a potential fusion add to fus_list
            print("comparing inv_fus of " + str(inv_fus_count) + " with " + str((len(sortedArr[i+1:]))))
            if i == len(sortedArr)-1:
                if inv_fus_count < 1:
                    print("adding i protein (last i protein in list)")
                    fus_list.append(sortedArr[i])
            elif inv_fus_count < (len(sortedArr[i+1:])):
                print("adding i protein (general)")
                fus_list.append(sortedArr[i])
            else:
                print("not adding i protein (rejected)")

    tot_length = 0

    #checks for extra fusions if proteins are lined up end to end
    for i in range(len(fus_list)):
        tot_length += (fus_list[i]['send'] - fus_list[i]['sstart'])
    tot_length -= overlap_length
    
    
    if len(fus_list) < 2:
        print("Only ONE CANDIDATE identified")
        return []
    elif tot_length > MIN_FUSION_COVERAGE:
        print("output is: " + str(fus_list))
        return fus_list
    else:
        print("Didn't meet MIN_THRESHOLD")
        print("-------")
        return []
#print(isFusion(test_fus))



genDict(geneFusions, df)


#sorts the genomeFusions dictionary by the genome start index
num_out = 0
for id in geneFusions:
    if(len(geneFusions[id]) == 1):
        continue
    sortedGeneArr = sorted(geneFusions[id], key=lambda x: x['sstart'])
    # print(sortedGeneArr)
    # x = input()
    if len(isFusion(sortedGeneArr)) != 0:
        print("SUCCESS")
        num_out+=1
    #y = input()
print(str(num_out) + " fusions found")
