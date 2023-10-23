#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import numpy as np
import gzip


# In[7]:


raw_data_path = '/work/project/ext_014/GSE220442/rawdata/'
c2lpath = '/work/project/ext_014/GSE220442/cell2location/'
cspath = '/work/project/ext_014/GSE220442/CSIDE/'


# In[3]:


#Generation of TE list
repeat_list = []
unique_list = []

with gzip.open(raw_data_path+'rmsk.txt.gz', 'rt') as file:
    for line in file:
        columns = line.rstrip().split('\t')
        name = columns[10]
        classification = columns[11]
        entry = (name, classification)
        repeat_list.append(entry)  # Add each entry to the repeat_list

        if entry not in unique_list:
            unique_list.append(entry)  # Add unique entries to unique_list
            
scTE_list = []
for entry in unique_list:
    name, classification = entry
    scTE_list.append(name)


# In[4]:


#Get a normalized spot-gene matrix from visium anndata 
def get_N_vis_matrix(adata):
    vis_matrix = pd.DataFrame(adata.X.T.todense())
    vis_matrix.columns = adata.obs_names
    vis_matrix.index = adata.var_names
    N_vis_matrix = vis_matrix.div(vis_matrix.sum(axis=0), axis=1)
    return N_vis_matrix


# In[5]:


#Get the overall TE fraction
def get_gene_fraction(matrix,TE_list):
    TE_check = []
    for i in range(len(matrix)):
        if matrix.index[i] in TE_list:
            TE_check.append('TE')
        else:
            TE_check.append('non-TE')
    matrix['TE_check'] = TE_check
    #Fraction check
    all_spots_sum = matrix.groupby('TE_check').sum()
    all_spots_sum = all_spots_sum.sum(axis=1).to_frame(name="sub_total")
    Total = all_spots_sum.sum()
    for i in range(len(all_spots_sum)):
        all_spots_sum.iloc[i] = all_spots_sum.iloc[i]/Total
    for i in range(len(all_spots_sum)):
        if all_spots_sum.index[i] == 'non-TE':
            weight = all_spots_sum.iloc[i]
    weight=weight[0]
    return weight


# In[22]:


#Reference dataframe generation
adata = sc.read_10x_h5(raw_data_path+'mtg_nature_pure.h5')
celltype_anno = pd.read_csv(raw_data_path+'MTG_pure_annotation.csv.gz', index_col=0)
overlap_barcode = np.intersect1d(adata.obs.index.tolist(), celltype_anno.index.tolist())
celltype_anno = celltype_anno.loc[overlap_barcode, :]
adata = adata[overlap_barcode, :]
adata.obs['celltype'] = celltype_anno['celltype'].copy()
#Cell type profiling
cell_type = adata.obs['celltype'].unique()
matrix = pd.DataFrame(adata.X.T.todense())
ct_col = adata.obs['celltype'].values
matrix.columns = ct_col
for cell in cell_type:
    matrix[cell+'_sum'] = matrix.filter(like=cell).sum(1)
    matrix = matrix.drop(columns=cell)
matrix.index = adata.var_names


# In[45]:


#Check if genes in reference contains repeats/TEs in spatial TE dataframe, if yes, remove from reference
adata_vis = sc.read_h5ad(c2lpath+'saved_h5ad/spatial_te_with_cluster_anno.h5ad')
adata_vis.obs.index.name = None
adata_vis.obs.index = adata_vis.obs_names.str[:19]
vis_mat = get_N_vis_matrix(adata_vis)
get_gene_fraction(vis_mat,scTE_list)
TE_vis_mat = vis_mat[vis_mat['TE_check']=='TE']
TE_vis_mat = TE_vis_mat.drop(columns=['TE_check'])
TE_vis_mat = TE_vis_mat.iloc[:,0:7]
for genes in matrix.index:
    if genes in TE_vis_mat.index:
        matrix.drop(labels = genes, axis=0, inplace=True)
#Normallize reference matrix and assign weight
N_matrix = matrix.div(matrix.sum(axis=0), axis=1)
N_matrix = N_matrix*get_gene_fraction(vis_mat,scTE_list)


# In[46]:


#Combine the TE rows from spatial data with the reference matrix and assigned the same value for cell-type specific TE expressions (i.e. regular genes + TEs per cell type = 1)
n_slot = TE_vis_mat.shape[0] * TE_vis_mat.shape[1]
remain_fraction = (1-N_matrix.iloc[:, [0]].sum()[0])*N_matrix.shape[1]
TE_vis_mat[TE_vis_mat<1] = remain_fraction/n_slot
TE_vis_mat.columns = N_matrix.columns
combined = pd.concat([N_matrix, TE_vis_mat], axis=0)
combined.columns = cell_type


# In[48]:


combined.to_csv(cspath+'Reference/cell_type_profiles.csv.gz',compression='gzip')


# In[ ]:




