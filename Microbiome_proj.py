
#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

import pandas as pd
import numpy as np

# construct df with all data from results.tsv
df = pd.read_table('results.tsv')
# print columns of df for dev use in constructing filtered_df later
print('---------')
print(df.columns)
print('---------')

# construct filtered_df with only relevant data points for each protein matche
filtered_df = df[['#Query_id', '%_identity', 'e-value', 'Q_start', 'Q_end', 'S_start', 'S_end', 
'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
#print(filtered_df)


# for index, row in filtered_df.iterrows():
#     #print(row['#Query_id'])

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

#print(filtered_df[2])
def find_fusion(row):
    q_start = int(row[9])
    q_end = int(row[10])

    s_start = int(row[11])
    s_end = int(row[12])

    if q_start < s_start and q_end > s_end:
        return True
    elif s_start < q_start and s_end > q_end:
        return True

    return False

# For each query this function will return common pFam domains between the query and the subject
def PfamDoms(row):
    doms = []
    if isfloat(row[20]):
        return doms
    elif isfloat(row[21]):
        return doms
        
    q_pfam = row[20].split(',')
    s_pfam = row[21].split(',')

    print('q: ')
    print(q_pfam)

    print('s: ')
    print(s_pfam)

    for q_domain in q_pfam:
        for s_domain in s_pfam:
            if q_domain == s_domain:
                doms.append(q_domain)
    #print(common_doms)
    common_doms = [*set(doms)]
    print('res')
    print(common_doms)
    return common_doms


for index, row in df.iterrows():
    #print(PfamDoms(row))
    PfamDoms(row)


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