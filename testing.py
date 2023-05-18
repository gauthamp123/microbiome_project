import pandas as pd
df = pd.read_csv('Red.tsv', sep='\t')
df = df.fillna('NA')
NA_count = df.eq('NA').sum(axis=1).rename('NA_count')
df['na_sort'] = NA_count
df = df.sort_values(by=['Hit_tcid','na_sort']).drop(columns=['na_sort'])
df.to_csv('Red.1.tsv', sep='\t', index=False)