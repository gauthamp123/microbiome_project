import xml.etree.ElementTree as ET
import xml
import pandas as pd

xml_format = "GCF_009648975.1/xml/"
tsv_file = 'GCF_009648975.1/results.tsv'
df = pd.read_table(tsv_file)
file_to_find = df['#Query_id']

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
        


