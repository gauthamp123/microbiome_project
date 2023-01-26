#!/usr/bin/env python
# coding: utf-8

# Take the results.tsv file and parse it

# In[10]:


import pandas as pd
import numpy as np

hash_map = {}


# In[1]:


df = pd.read_table('results.tsv')
#print(df)

print(df.columns)
#filtered_df = df[['#Query_id, 'e-value', 'Query_Coverage', 'Hit_Coverage', 'Query_Pfam']]
filtered_df = df[['#Query_id', 'Query_Coverage', 'e-value', 'Hit_Coverage', 'Query_Pfam']]
print(filtered_df)


# In[ ]:




