import pandas as pd
import numpy as np
import os

df = pd.read_table('results.tsv')

geneFusions = {}

# compare overlap between individual fusion candidates to ensure there is not over a 20% overlap 
# (if an individual candidate overlaps with all others with greater than thresh, reject)
MAX_OVERLAP = 20

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
    for i in range(len(sortedArr) - 1):
        

#genDict(geneFusions, df)


#sorts the genomeFusions dictionary by the genome start index
for id in geneFusions:
    if(len(geneFusions[id]) == 1):
        continue
    sortedGeneArr = sorted(geneFusions[id], key=lambda x: x['sstart'])
    print(sortedGeneArr)
    x = input()

#print(genomeFusions)
