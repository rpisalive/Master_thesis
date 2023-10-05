#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
from anndata import AnnData
import numpy as np


# In[2]:


#path definition
te_path = '/work/project/ext_014/EMTAB10972/snakerangerout/scte/'
raw_data_path = '/work/project/ext_014/EMTAB10972/rawdata/'
reference_path = '/work/project/ext_014/EMTAB10972/CSIDE/Reference/'


# In[3]:


#Concatenate scTE objects
sc_meta = pd.read_csv(raw_data_path+'EMTAB10974_annotation.csv.gz')
sample = pd.read_csv(raw_data_path+'EMTAB10974_sample_name.csv.gz').iloc[:,0].unique()
def concatenate_teobjects(sample_name, path=te_path):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    teobj = sc.read_h5ad(path+sample_name+'.h5ad')
    teobj.obs['sample'] = sample_name
    # add sample name to obs names
    teobj.obs_names = teobj.obs["sample"]+ '_' +teobj.obs_names
    teobj.obs.index.name = 'barcodes'
    return teobj
#append scte objedcts into a list
te_objects = []
for i in sample:
    te_objects.append(concatenate_teobjects(i, path=te_path))
#combined all_te_objects
sample_data = pd.read_csv(raw_data_path+'EMTAB10974_sample_name.csv.gz') #to read it a list of sample name, creaetd manually, make sure the name is matching with the names in sc_meta
sc_te = sc.concat(te_objects)
sc_te.obs.index = sc_te.obs.index.str[:-2]
#Matching the size
sc_te.obs['in_meta'] =  sc_te.obs.index.isin(sc_meta['barcode'])
sc_te = sc_te[sc_te.obs['in_meta']==True]


# In[4]:


#Preparation for spacexr RCTD run
barcodes = sc_te.obs_names.tolist()
features = sc_te.var_names.tolist()
matrix = pd.DataFrame(sc_te.X.T.todense())
matrix = matrix.set_axis(barcodes, axis=1)
matrix = matrix.set_axis(features, axis=0)
matrix.to_csv(reference_path+'sc_counts.csv.gz',compression='gzip')


# In[ ]:




