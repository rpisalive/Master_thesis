#!/usr/bin/env python
# coding: utf-8

# In[37]:


import scanpy as sc
from scipy import io
import pandas as pd
import numpy as np
import gzip


# In[9]:


get_ipython().system('mkdir CT1')
get_ipython().system('mkdir CT2')
get_ipython().system('mkdir CT3')
get_ipython().system('mkdir AD1')
get_ipython().system('mkdir AD2')
get_ipython().system('mkdir AD3')
get_ipython().system('mkdir Reference')


# In[38]:


c2lpath = '/work/project/ext_014/GSE220442/cell2location/'
cspath = '/work/project/ext_014/GSE220442/CSIDE/'
rdpath = '/work/project/ext_014/GSE220442/rawdata/'


# In[39]:


#spot must be named as barcodes and barcodes must not contain '.' to fit the function read.SpatialRNA in spacexr
adata = sc.read_h5ad(c2lpath+'saved_h5ad/spatial_te_with_cluster_anno.h5ad')
adata.obs.index.name = 'barcodes'
adata.obs.index = adata.obs_names.str[:19]


# In[40]:


#remove noise
adata = adata[adata.obs['Layer']!='Noise']


# In[41]:


#for the import_weights function
sample = ['CT1','CT2','CT3','AD1','AD2','AD3']

#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names] = adata[adata.obs['sample'] == names]
    globals()[names] = globals()[names].obs.drop(columns = ['in_tissue','array_row',
                                                            'array_col','sample','in_ori','_indices','_scvi_batch','_scvi_labels','total RNA counts',
                                                            'RNA detection sensitivity (y_s)','leiden','region_cluster','in_tdata','Layer'])
    for i in range(len(globals()[names])):
        globals()[names]['others'] = globals()[names].iloc[i]['Total cell abundance (sum_f w_sf)']-globals()[names].iloc[i]['Astro']-globals()[names].iloc[i]['Endo']-globals()[names].iloc[i]['Exc']-globals()[names].iloc[i]['Inh']-globals()[names].iloc[i]['Micro']-globals()[names].iloc[i]['OPC']-globals()[names].iloc[i]['Oligo']
    globals()[names] = globals()[names].drop(columns = ['Total cell abundance (sum_f w_sf)'])
    globals()[names] = globals()[names].div(globals()[names].sum(axis=1), axis=0)
for names in sample:
    globals()[names].to_csv(cspath+names+'/ct_proportion.csv.gz', compression='gzip')


# In[42]:


#Explanatory variables preparation for C-SIDE
for names in sample:
    globals()[names+'_spte'] = adata[adata.obs['sample'] == names]
for names in sample:
    bar_reg = globals()[names+'_spte'].obs['Layer'].to_frame(name = 'regions')
    bar_reg.to_csv(cspath+names+'/build_designmatrix_regions.csv.gz',compression='gzip')


# In[43]:


#Counts matrix preparation
for names in sample:
    barcodes = globals()[names+'_spte'].obs_names.tolist()
    features =  globals()[names+'_spte'].var_names.tolist()
    globals()[names] = pd.DataFrame(globals()[names+'_spte'].X.T.todense())
    globals()[names] = globals()[names].set_axis(barcodes, axis=1)
    globals()[names] = globals()[names].set_axis(features, axis=0)
    globals()[names].to_csv(cspath+names+'/counts.csv.gz',compression='gzip')


# In[44]:


#Location preparation
for names in sample:
    position = globals()[names+'_spte'].obs[['array_row','array_col']]
    position.to_csv(cspath+names+'/location.csv.gz',compression='gzip')


# In[ ]:




