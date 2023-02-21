import pandas as pd
import numpy as np
import os

df = pd.read_table('results.tsv')

genomeFusions = {}

def genDict(dict, input):
    for index, row in input.iterrows():
        new_entry = {"query": row["#Query_id"], "qcov": row["Query_Coverage"], "sstart": row["S_start"], "send": row["S_end"], "scov": row["Hit_Coverage"]}
        tcid = row["Hit_tcid"] + "-" + row["Hit_xid"]
        if tcid in dict:
            dict.get(tcid).append(new_entry)
            print("Appending : " + tcid)
        else:
            dict[tcid] = [new_entry]
            print("New entry")

genDict(genomeFusions, df)