import pandas as pd
import numpy as np
import os

df = pd.read_table('results.tsv')

geneFusions = {}
test_fus = {'query': 'YP_501170.1', 'qcov': 50, 'sstart': 1, 'send': 10, 'scov': 99.7}, {'query': 'YP_500726.1', 'qcov': 77.5, 'sstart': 6, 'send': 20, 'scov': 87.1}


# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_OVERLAP = 20
MAX_THRESHOLD = 80
MIN_THRESHOLD = 10

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

def isFusion(sortedArr):
    fus_list = []
    overlap_length = 0

    
    for i in range(len(sortedArr)):
        if sortedArr[i]['qcov'] > MAX_THRESHOLD:
            return []
        else:
            inv_fus_count = 0
            for j in range(i+1, len(sortedArr)):
                # This looks for proteins that have too much overlap with the current protein
                if sortedArr[i]['send'] >= sortedArr[j]['send']:
                    inv_fus_count += 1
                elif ((sortedArr[i]['send'] - sortedArr[j]['sstart']) / (sortedArr[j]['send'] - sortedArr[i]['sstart'])) * 100 < MAX_OVERLAP:
                    inv_fus_count += 1
                elif ((sortedArr[i]['send'] - sortedArr[j]['sstart']) / (sortedArr[i]['sstart'] - sortedArr[j]['send'])) * 100 == 0:
                    i += 1
                else:
                    overlap_length += abs(sortedArr[i]['send'] - sortedArr[j]['sstart'])
            if inv_fus_count < len(sortedArr) - i:
                fus_list.append(sortedArr[i])

    tot_length = 0

    for i in range(len(fus_list)):
        tot_length += (fus_list[i]['send'] - fus_list[i]['sstart'])
    tot_length -= overlap_length
    
    if tot_length > MIN_THRESHOLD:
        return fus_list
    else:
        return []
    

    
print(isFusion(test_fus))
x = input()        

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

#print(genomeFusions)
