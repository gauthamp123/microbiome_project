import xml.etree.ElementTree as ET
import xml
import pandas as pd
from mmseqs_hmmtop import overlapDict
import pickle

# with open("GCF_009648975.1/hmmtop.db", "rb") as file:
#     # Deserialize the data using pickle.load()
#     data = pickle.load(file)

<<<<<<< Updated upstream
dictionary = {}
for row in df.itertuples(index=False):
    query_id = row._0
    values = row._asdict()
    del values['_0']
    
    # Add the values to the dictionary
    if query_id not in dictionary:
        dictionary[query_id] = {}
    
    dictionary[query_id].update(values)

sequences = {}
for key in dictionary:
    xmlfile = xml_format + key + '.xml'
    tree = ET.parse(xmlfile)
    root = tree.getroot()

    # ['gnl', 'TC-DB', 'P60778', '2.A.1.7.14 Protein tsgA OS=Escherichia coli (strain K12) GN=tsgA PE=1 SV=1']
    for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
        hit_info = item.find('Hit_def').text.split('|')
        #print(hit_info[3].split(' ')[0])

        if hit_info[2] == dictionary[key]['Hit_xid'] and hit_info[3].split(' ')[0] == dictionary[key]['Hit_tcid']:
            j = item.findall('Hit_hsps/Hsp')
            for h_item in j:
                query_seq = h_item.find('Hsp_qseq').text
                subject_seq = h_item.find('Hsp_hseq').text
                sequences[key] = {'query_sequence': query_seq, 'subject_sequence': subject_seq}

print(sequences)
        
=======
# # Extracts all query data and target data from hmmtop file
# query_data = data['queries']
# target_data = {}
# for key in data['tcdb']:
#     key_elems = key.split('|')
#     new_key = key_elems[3] + '-' + key_elems[2]

#     target_data[new_key] = data['tcdb'][key]
>>>>>>> Stashed changes



# xml_format = "GCF_009648975.1/xml/"
# tsv_file = 'GCF_009648975.1/results.tsv'
# df = pd.read_table(tsv_file)
# file_to_find = df['#Query_id']

def parse(hmmtop_file, xml_dir, results_file, minRes):
    with open(hmmtop_file, "rb") as file:
        # Deserialize the data using pickle.load()
        data = pickle.load(file)

    # Extracts all query data and target data from hmmtop file
    query_data = data['queries']
    target_data = {}
    for key in data['tcdb']:
        key_elems = key.split('|')
        new_key = key_elems[3] + '-' + key_elems[2]

        target_data[new_key] = data['tcdb'][key]



    xml_format = xml_dir
    tsv_file = results_file
    df = pd.read_table(tsv_file)
    file_to_find = df['#Query_id']

    dictionary = {}


    sequences = {}
    mmseqsDict = {}
    hmmTopDict = {}
    # for each id in the dataframe results.tsv
    for row in df.itertuples(index=False):
        key = row._0
        hmm_key = row._0
        xmlfile = xml_format + key + '.xml'
        tree = ET.parse(xmlfile)
        root = tree.getroot()

        # ['gnl', 'TC-DB', 'P60778', '2.A.1.7.14 Protein tsgA OS=Escherichia coli (strain K12) GN=tsgA PE=1 SV=1']
        # goes through each branch of the xml file that contains 'hit'
        for item in root.findall('./BlastOutput_iterations/Iteration/Iteration_hits/Hit'):
            hit_info = item.find('Hit_def').text.split('|')
            # looks for the specific target accession that correlates to the results.tsv data
            if hit_info[2] == row.Hit_xid and hit_info[3].split(' ')[0] == row.Hit_tcid:
                j = item.findall('Hit_hsps/Hsp')
                # if it exists, then look to see if there is alignment sequences in xml
                for h_item in j:
                    query_seq = h_item.find('Hsp_qseq').text
                    subject_seq = h_item.find('Hsp_hseq').text
                    target_id = row.Hit_tcid + '-' + row.Hit_xid
                    # only include in mmseqs if the target accession and query exists in results.tsv
                    if key in query_data and target_id in target_data:
                        mmseqsDict[key] = {'qaln': query_seq, 'taln': subject_seq, 'target': target_id,
                                        'qstart': row.Q_start, 'qend': row.Q_end, 'tstart': row.S_start, 'tend': row.S_end}
                    
                    
                    # puts the overlaps in the format needed by overlapDict
                    qtms = {}
                    if key in query_data:
                        qtms['tms'] = list(query_data[key].values())
                        if key not in hmmTopDict:
                            hmmTopDict[key] = qtms


                    ttms = {}
                    if target_id in target_data:
                        ttms['tms'] = list(target_data[target_id].values())
                        if target_id not in hmmTopDict:
                            hmmTopDict[target_id] = ttms

        # there are some keys in df that are not in query or target data what to do about those?

    overlap_dict = overlapDict(mmseqsDict, hmmTopDict, minRes)
    return overlap_dict

hmmtop_file, xml_dir, results_file, minRes = "GCF_009648975.1/hmmtop.db", "GCF_009648975.1/xml/", 'GCF_009648975.1/results.tsv', 8
parse(hmmtop_file, xml_dir, results_file, minRes)

'''
#For Testing
tempMmseqs = {}
tempData = {}
tempData['qaln'] = 'MGFDIGGDIGKPLKDAFDKFGADIKMTFLTVLNWMK--WISIG------ILIVISVI-------LICKIIKVLFQCGKCLLSCFGFCKK'
tempData['taln'] = 'MGFSINFD---PIINKFREFQTNINHNINEQLDKLKMVWINLGSHIKYWFIIIISILTILFILFLLIKITKLILNCKKIFSCCCNVCCK'
tempData['target'] = '1.A.100.1.1-B2X7D9'
tempData['qstart'] = 1
tempData['qend'] = 74
tempData['tstart'] = 22
tempData['tend'] = 107
tempMmseqs['1.A.95.2.5-YP_009361958'] = tempData

tempHmmtop = {}
tempQhmmData = {}
tempQhmmData['tms'] = [[26, 49], [54,72]]
tempThmmData = {}
tempThmmData['tms'] = [[68, 92]]
tempHmmtop['1.A.95.2.5-YP_009361958'] = tempQhmmData
tempHmmtop['1.A.100.1.1-B2X7D9'] = tempThmmData
overlap_dict = overlapDict(tempMmseqs, tempHmmtop, minRes)
'''