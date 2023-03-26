import pandas as pd
import numpy as np
import os

df = pd.read_table('results.tsv')

geneFusions = {}
test_fus = [{'query': 'YP_501037.1', 'qcov': 63.9, 'sstart': 242, 'send': 356, 'scov': 31.8}, {'query': 'YP_501311.1', 'qcov': 72.0, 'sstart': 247, 'send': 356, 'scov': 30.4}]

# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_OVERLAP = 20
MAX_THRESHOLD = 80
MIN_THRESHOLD = 10

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

    # output list
    fus_list = []
    overlap_length = 0
    
    # iteratate through each fusion candidate 
    for i in range(len(sortedArr)):
        if sortedArr[i]['qcov'] > MAX_THRESHOLD:
            return []
        else:
            print("Analyzing i protein: " + sortedArr[i]["query"])
            # invalid fusion counter
            inv_fus_count = 0
            for j in range(i+1, len(sortedArr)):
                # This looks for proteins that have too much overlap with the current protein
                overlap_size = sortedArr[i]['send'] - sortedArr[j]['sstart']
                # vars represent the percent overlap
                perc_protein1 = (overlap_size / (sortedArr[i]['send'] - sortedArr[i]['sstart'])) * 100
                perc_protein2 = (overlap_size / (sortedArr[j]['send'] - sortedArr[j]['sstart'])) * 100
                
                # making sure theres not total overlap between 2 proteins
                if sortedArr[i]['send'] >= sortedArr[j]['send']:
                    print("invalid identified")
                    inv_fus_count += 1
                # checking to see if the minimum overlap among comparable proteins is less than the threshold
                elif (min(perc_protein1, perc_protein2)) < MAX_OVERLAP:
                    inv_fus_count += 1
                # if theres no overlap move on
                elif ((sortedArr[i]['send'] - sortedArr[j]['sstart']) / (sortedArr[i]['sstart'] - sortedArr[j]['send'])) * 100 == 0:
                    break
                else:
                    overlap_length += overlap_size
            # if there is a potential fusion add to fus_list
            print("comparing inv_fus of " + str(inv_fus_count) + " with " + str((len(sortedArr[i+1:]))))
            if inv_fus_count < (len(sortedArr[i+1:])):
                print("adding i protein ")
                fus_list.append(sortedArr[i])
            else:
                print("not adding i protein")

    tot_length = 0

    #checks for extra fusions if proteins are lined up end to end
    for i in range(len(fus_list)):
        tot_length += (fus_list[i]['send'] - fus_list[i]['sstart'])
    tot_length -= overlap_length
    
    if tot_length > MIN_THRESHOLD:
        return fus_list
    else:
        return []
    
print(isFusion(test_fus))


'''
genDict(geneFusions, df)


#sorts the genomeFusions dictionary by the genome start index
for id in geneFusions:
    if(len(geneFusions[id]) == 1):
        continue
    sortedGeneArr = sorted(geneFusions[id], key=lambda x: x['sstart'])
    # print(sortedGeneArr)
    # x = input()
    if len(isFusion(sortedGeneArr)) != 0:
        print(isFusion(sortedGeneArr))
    #y = input()

'''