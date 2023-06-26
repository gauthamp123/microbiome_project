import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns
import seaborn.objects as so

'''

TEST CASES:

Reject, Unsure (Cake shape) -
test_input = [{'query': 'YP_499975.1', 'qcov': 92.9, 'sstart': 258, 'send': 1229, 'scov': 55.9, 'hit_length': 1738}, 
                {'query': 'YP_499935.1', 'qcov': 87.6, 'sstart': 505, 'send': 839, 'scov': 19.2, 'hit_length': 1738}, 
                {'query': 'YP_499741.1', 'qcov': 80.5, 'sstart': 642, 'send': 1582, 'scov': 54.1, 'hit_length': 1738}, 
                {'query': 'YP_500508.1', 'qcov': 93.0, 'sstart': 757, 'send': 1020, 'scov': 15.1, 'hit_length': 1738}]

Reject (low coverage) -
test_input = [{'query': 'YP_499455.1', 'qcov': 90.2, 'sstart': 152, 'send': 310, 'scov': 50.5, 'hit_length': 313}, {'query': 'YP_499456.1', 'qcov': 86.5, 'sstart': 187, 'send': 310, 'scov': 39.3, 'hit_length': 313}]

Reject (should be single component) how to reject this case without the rare 90-10 case?
test_input = [{'query': 'YP_500738.1', 'qcov': 85.4, 'sstart': 14, 'send': 221, 'scov': 89.6, 'hit_length': 231}, {'query': 'YP_500401.1', 'qcov': 85.3, 'sstart': 34, 'send': 224, 'scov': 82.3, 'hit_length': 231}]

Reject (same as previous)
test_input = [{'query': 'YP_498717.1', 'qcov': 94.7, 'sstart': 22, 'send': 595, 'scov': 89.5, 'hit_length': 640}, {'query': 'YP_501198.1', 'qcov': 86.6, 'sstart': 282, 'send': 530, 'scov': 38.8, 'hit_length': 640}]

Reject (same as previous)
test_input = [{'query': 'YP_501493.1', 'qcov': 95.4, 'sstart': 10, 'send': 133, 'scov': 87.9, 'hit_length': 140}, {'query': 'YP_501494.1', 'qcov': 89.7, 'sstart': 10, 'send': 113, 'scov': 73.6, 'hit_length': 140}]



'''

#test_input = [{'query': 'YP_501406.1', 'qcov': 96.9, 'sstart': 84, 'send': 420, 'scov': 59.2, 'hit_length': 568}, {'query': 'YP_499273.1', 'qcov': 90.9, 'sstart': 87, 'send': 419, 'scov': 58.5, 'hit_length': 568}, {'query': 'YP_499225.1', 'qcov': 90.8, 'sstart': 108, 'send': 422, 'scov': 55.3, 'hit_length': 568}]

def fusion_dist(fusion_candidate_list):

    plt.cla()
    y_count = len(fusion_candidate_list) + 0.5

    plt.plot([0,fusion_candidate_list[0]['hit_length']],[y_count+1,y_count+1],label = "TCDB Reference" + "(0-" + str(fusion_candidate_list[0]['hit_length']) + ")")

    coverage_list = []
    for candidate in fusion_candidate_list:
        coverage_list = coverage_list + (list(range(candidate['sstart'], candidate['send'])))
        plt.plot([candidate['sstart'],candidate['send']], [y_count, y_count], label = candidate['query'] + "(" + str(candidate['sstart']) + "-" + str(candidate['send']) + ")")
        y_count -= 1

    # print("COV LIST")
    # print(coverage_list)
    # print("-------")

    # print("FREQ LIST")
    frequency_list = {x:coverage_list.count(x) for x in coverage_list}
    # print(frequency_list)
    # print("-------")
    next = input("Next output")

    plt.plot(frequency_list.keys(), frequency_list.values(), label = "Frequency Distribution")
    plt.legend()
    plt.show()

df = pd.read_table('GCF_009648975.1/results.tsv')

geneFusions = {}
test_fus =  [{'query': 'A', 'qcov': 69.6, 'sstart': 295, 'send': 713, 'scov': 58.2}, {'query': 'B', 'qcov': 65.4, 'sstart': 405, 'send': 709, 'scov': 42.3}, {'query': 'C', 'qcov': 41.4, 'sstart': 481, 'send': 717, 'scov': 32.9}]

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

    # print("-------")
    # print("Input is " + str(sortedArr))

    # output list
    fus_list = []
    overlap_length = 0
    
    # iteratate through each fusion candidate 
    for i in range(len(sortedArr)):
        if sortedArr[i]['scov'] > MAX_SUBJECT_COV_THRESHOLD:
            # print("i rejected, candidate scov too long")
            continue
        elif sortedArr[i]['scov'] < MIN_SUBJECT_COV_THRESHOLD:
            # print("i rejected, candidate scov was too small")
            continue
        elif sortedArr[i]['qcov'] < MIN_QUERY_COV_THRESHOLD:
            # print("i rejected, candidate qcov was too small")
            continue
        else:
            # print("Accepting i protein: " + sortedArr[i]["query"])
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
        # print("No candidates identified")
        return []
    elif len(fus_list) < 2:
        # print("Insufficient Candidates identified")
        # print(fus_list)
        return []
    elif (tot_length/fus_list[0]['hit_length'])*100 < MIN_FUSION_COVERAGE:
        # print("Didn't meet MIN_THRESHOLD")
        # print("total length: " + str(tot_length))
        # print("hit length: " + str(fus_list[0]['hit_length']))
        # print("-------")
        return []
    else:
        # print("output is: " + str(fus_list))
        # print("total length: " + str(tot_length))
        # print("hit length: " + str(fus_list[0]['hit_length']))
        return fus_list


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
    # num_in+=1
    # if len(isFusion(sortedGeneArr)) != 0:
    #     # print(isFusion(sortedGeneArr))
    #     # print("SUCCESS - Would you like to see the graph? (Y/N)")
    #     isGraph = input()
    #     if(isGraph == "Y"):
    #         fusion_dist(isFusion(sortedGeneArr))
    #     num_out+=1
        # print("FAIL")
        #isGraph = input()
        #if(isGraph == "Y"):
            #fusion_dist(isFusion(sortedGeneArr))
    #y = input()
# print("-------")
# print(str(num_out) + " fusions found out of " + str(num_in))
