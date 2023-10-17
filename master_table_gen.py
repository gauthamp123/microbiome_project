from genome_comparison import *
import pandas as pd
import os

genomes = []

# generates the columns needed for our Master Dataframe
def master_df_column_gen(directory):
    columns = ['#TCID', 'Acc', 'CE', 'Role', 'hit_tms_no']
    per_genome_columns = ['query', 'q_tms', 'evalue', 'pident', 'qcov', 'scov']
    num_genomes = 0
    items = os.listdir(directory)
    
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]
    columns.extend(genomes)
    for genome in genomes:
        columns.extend(per_genome_columns)

    return columns

master_dict = {}

# TODO: basically the main script LOL
def master_dict_generation(green_dir):
    for filename in os.listdir(green_dir):
        if os.path.isfile(os.path.join(green_dir, filename)):

            file_parts = filename.split('_')[:2]
            genome = '_'.join(file_parts)

            green_df = pd.read_csv(os.path.join(green_dir, filename), delimiter='\t')
            
            for index, row in green_df.iterrows():
                tcid_acc = row['Hit_tcid'] + '-' + row['Hit_xid']
                if tcid_acc in master_dict and master_dict[tcid_acc]:
                    if genome in master_dict[tcid_acc] and master_dict[tcid_acc][genome]:
                        master_dict[tcid_acc][genome][row['#Query_id']] = {}
                    else:
                        master_dict[tcid_acc][genome] = {row['#Query_id']:{}}
                else:
                    
                    master_dict[tcid_acc] = {genome: {row['#Query_id']:{}}}

    for key, value in master_dict.items():
        print(key)
        print(value)

df_columns = master_df_column_gen('test_comparisons')

# necessary calls as pre-reqs for gtoms script (prolly gonna have to change the dir)
#getSmithWaterman('/Users/gautham/microbiome_project/test_comparisons')
#parse_sw('/Users/gautham/microbiome_project/test_comparisons')

master_dict_generation('test_comparisons/greens')
