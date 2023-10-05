#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import random


# In[2]:


region_to_check = 'nc'


# In[3]:


data_path = '/work/project/ext_014/EMTAB10972/CSIDE/vol_plot_data/'
region_list = ['r','nc','dm','t','cp','m','e']
cell_type = ['mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','mg','retina']


# In[4]:


list_in_scTEref =  pd.read_csv(data_path+region_to_check+'/'+region_to_check+'_combined_scTEref_iw.csv')
list_in_syn =  pd.read_csv(data_path+region_to_check+'/'+region_to_check+'_combined_im.csv')
list_in_scTEref_sig = list_in_scTEref[list_in_scTEref['sig_cat'] == 'sig']
list_in_syn_sig = list_in_syn[list_in_syn['sig_cat'] == 'sig']
list_in_scTEref_nsig = list_in_scTEref[list_in_scTEref['sig_cat'] == 'non_sig']
list_in_syn_nsig = list_in_syn[list_in_syn['sig_cat'] == 'non_sig']
for ct in cell_type:
    globals()[ct+'_scTEref_sig'] = list_in_scTEref_sig[list_in_scTEref_sig['cell_name'] == ct]
    globals()[ct+'_syn_sig'] = list_in_syn_sig[list_in_syn_sig['cell_name'] == ct]
    globals()[ct+'_scTEref_sig'] = globals()[ct+'_scTEref_sig']['genes'].tolist()
    globals()[ct+'_syn_sig'] = globals()[ct+'_syn_sig']['genes'].tolist()
for ct in cell_type:
    globals()[ct+'_scTEref_nsig'] = list_in_scTEref_nsig[list_in_scTEref_nsig['cell_name'] == ct]
    globals()[ct+'_syn_nsig'] = list_in_syn_nsig[list_in_syn_nsig['cell_name'] == ct]
    globals()[ct+'_scTEref_nsig'] = globals()[ct+'_scTEref_nsig']['genes'].tolist()
    globals()[ct+'_syn_nsig'] = globals()[ct+'_syn_nsig']['genes'].tolist()


# In[5]:


len(list_in_scTEref_sig)


# In[6]:


len(list_in_syn_sig)


# In[7]:


#False -ve
in_scTEref_sig_not_in_syn_sig = []
for ct in cell_type:
    for element in globals()[ct+'_scTEref_sig']:
        if element not in globals()[ct+'_syn_sig']:
            in_scTEref_sig_not_in_syn_sig.append(element)
len(in_scTEref_sig_not_in_syn_sig)


# In[8]:


#False +ve
in_syn_sig_not_in_scTEref_sig = []
for ct in cell_type:
    for element in globals()[ct+'_syn_sig']:
        if element not in globals()[ct+'_scTEref_sig']:
            in_syn_sig_not_in_scTEref_sig.append(element)
len(in_syn_sig_not_in_scTEref_sig)


# In[75]:


#True -ve
in_scTEref_nsig_in_syn_nsig = []
for ct in cell_type:
    for element in globals()[ct+'_scTEref_nsig']:
        if element in globals()[ct+'_syn_nsig']:
            in_scTEref_nsig_in_syn_nsig.append(element)
len(in_scTEref_nsig_in_syn_nsig)


# In[76]:


#True +ve
in_scTEref_sig_in_syn_sig = []
for ct in cell_type:
    for element in globals()[ct+'_scTEref_sig']:
        if element in globals()[ct+'_syn_sig']:
            in_scTEref_sig_in_syn_sig.append(element)
len(in_scTEref_sig_in_syn_sig)


# In[ ]:




