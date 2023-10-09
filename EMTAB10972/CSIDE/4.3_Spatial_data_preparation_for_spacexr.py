#!/usr/bin/env python
# coding: utf-8

# In[10]:


import scanpy as sc
from scipy import io
import pandas as pd
import gzip


# In[8]:


get_ipython().system('mkdir S1')
get_ipython().system('mkdir S2')
get_ipython().system('mkdir S3')
get_ipython().system('mkdir Reference')


# In[11]:


c2lpath = '/work/project/ext_014/EMTAB10972/cell2location/'
cspath = '/work/project/ext_014/EMTAB10972/CSIDE/'


# In[12]:


#spot must be named as barcodes and barcodes must not contain '.' to fit the function read.SpatialRNA in spacexr
adata = sc.read_h5ad(c2lpath+'saved_h5ad/spatial_te_with_cluster_anno.h5ad')
adata.obs.index.name = 'barcodes'
adata.obs.index = adata.obs_names.str[:19]


# In[13]:


#for the import_weights function
sample = ['S1','S2','S3']

#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names] = adata[adata.obs['sample'] == names]
    globals()[names] = globals()[names].obs.drop(columns = ['in_tissue','array_row',
                                                            'array_col','sample','in_meta','_indices','_scvi_batch','_scvi_labels','total RNA counts',
                                                            'RNA detection sensitivity (y_s)','leiden','region_cluster','in_tdata','max_celltype','final_cluster'])
    for i in range(len(globals()[names])):
        globals()[names]['others'] = globals()[names].iloc[i]['Total cell abundance (sum_f w_sf)']-globals()[names].iloc[i]['choroid_plexus']-globals()[names].iloc[i]['dien_NPC']-globals()[names].iloc[i]['dien_neuron']-globals()[names].iloc[i]['epi']-globals()[names].iloc[i]['hind_NPC']-globals()[names].iloc[i]['hind_neuron']-globals()[names].iloc[i]['mesen']-globals()[names].iloc[i]['mg']-globals()[names].iloc[i]['nc']-globals()[names].iloc[i]['pns']-globals()[names].iloc[i]['retina']-globals()[names].iloc[i]['tel_NPC']-globals()[names].iloc[i]['tel_neuron']
    globals()[names] = globals()[names].drop(columns = ['Total cell abundance (sum_f w_sf)'])
    globals()[names] = globals()[names].div(globals()[names].sum(axis=1), axis=0)
for names in sample:
    globals()[names].to_csv(cspath+names+'/ct_proportion.csv.gz', compression='gzip')


# In[14]:


#Explanatory variables preparation for C-SIDE
adata.var_names_make_unique()
for names in sample:
    globals()[names+'_spte'] = adata[adata.obs['sample'] == names]
for names in sample:
    bar_reg = globals()[names+'_spte'].obs['final_cluster'].to_frame(name = 'regions')
    bar_reg.to_csv(cspath+names+'/build_designmatrix_regions.csv.gz',compression='gzip')


# In[15]:


#Counts matrix preparation of spatial data
for names in sample:
    barcodes = globals()[names+'_spte'].obs_names.tolist()
    features =  globals()[names+'_spte'].var_names.tolist()
    globals()[names] = pd.DataFrame(globals()[names+'_spte'].X.T.todense())
    globals()[names] = globals()[names].set_axis(barcodes, axis=1)
    globals()[names] = globals()[names].set_axis(features, axis=0)
    globals()[names].to_csv(cspath+names+'/counts.csv.gz',compression='gzip')


# In[16]:


#Location preparation of spatial data
for names in sample:
    position = globals()[names+'_spte'].obs[['array_row','array_col']]
    position.to_csv(cspath+names+'/location.csv.gz',compression='gzip')


# In[ ]:




