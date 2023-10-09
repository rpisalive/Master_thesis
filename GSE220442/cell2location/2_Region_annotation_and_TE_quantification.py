#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import logging, matplotlib, os, sys
import scanpy as sc
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
import scipy as sp
from anndata import AnnData
import numpy as np
import seaborn as sns
import anndata as ad
import math


# In[2]:


#path definition
h5ad_path = '/work/project/ext_014/GSE220442/cell2location/saved_h5ad/'
te_path = '/work/project/ext_014/GSE220442/snakerangerout/scte/'
raw_data_path = '/work/project/ext_014/GSE220442/rawdata/'


# In[3]:


#Concatenate scTE objects
adata_vis = sc.read_h5ad(h5ad_path+'combined_regional_cluster.h5ad')
sample = adata_vis.obs['sample'].unique()
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
    teobj.obs.index.name = 'spot_id'
    return teobj
#append scte objedcts into a list
te_objects = []
for i in sample:
    te_objects.append(concatenate_teobjects(i, path=te_path))
#combined all_te_objects
sample_data = pd.read_csv(raw_data_path+'GSE220442visium.csv.gz')
sp_te = te_objects[0].concatenate(
    te_objects[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)
#Matching the size
sp_te.obs['in_adata'] =  sp_te.obs.index.isin(adata_vis.obs.index)
sp_te = sp_te[sp_te.obs['in_adata']==True]
if len(sp_te.obs)<len(adata_vis.obs):
        adata_vis.obs['in_tdata'] = adata_vis.obs.index.isin(sp_te.obs.index)
        adata_vis = adata_vis[adata_vis.obs['in_tdata'] == True]
#copy all other layers
sp_te.obs=adata_vis.obs.copy()
sp_te.uns=adata_vis.uns.copy()
sp_te.obsm=adata_vis.obsm.copy()
sp_te.obsp=adata_vis.obsp.copy()
res = []
for sub in sp_te.obs_names:
    modified_sub = sub[:-2]
    res.append(modified_sub)
sp_te.obs_names = res


# In[4]:


#Copy annotation from original study
meta = pd.read_csv(raw_data_path+'GSE220442_metadata.csv.gz',index_col =0 )
patient_ID_sample = {'2_5': 'CT3','1_1':'CT1','18_64':'CT2','2_3':'AD1','2_8':'AD2','T4857':'AD3'}
sam_name = []
for i in range(len(meta)):
    sam_name.append(patient_ID_sample[meta.iloc[i]['patientID']])
meta['sam_name'] = sam_name
meta.index = meta["sam_name"]+ '_' +meta.index
res = []
for sub in meta.index:
    modified_sub = sub[0:20]
    res.append(modified_sub)
meta.index = res
meta['in_adata'] = meta.index.isin(sp_te.obs.index)
meta = meta[meta['in_adata']==True]
sp_te.obs['Layer'] = meta['Layer']


# In[5]:


for names in sample:
    ts = sp_te[sp_te.obs['sample'] == names].copy()
    sc.pl.spatial(ts, color=['Layer'], size=1.3, img_key='hires',library_id=names,alpha=0.5,title=names)


# In[6]:


sp_te.write(h5ad_path+'spatial_te_with_cluster_anno.h5ad')


# In[7]:


sp_te.var["mt"] = sp_te.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(sp_te, qc_vars=["mt"], inplace=True)
sc.pl.violin(sp_te, ['n_genes_by_counts', 'total_counts','pct_counts_mt'], groupby='sample',rotation= 45)


# In[8]:


#Splitting the whole object into individual tissue object
sp_te = sp_te[sp_te.obs['Layer'] != 'Noise']
sp_te.var_names_make_unique()
for names in sample:
    globals()[names+'te'] = sp_te[sp_te.obs['sample'] == names]


# In[9]:


#CT1 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: CT1")
sns.distplot(CT1te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(CT1te.obs["total_counts"][CT1te.obs["total_counts"] < 3000], kde=False, bins=40, ax=axs[1])
sns.distplot(CT1te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(CT1te.obs["n_genes_by_counts"][CT1te.obs["n_genes_by_counts"] < 2000], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(CT1te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {CT1te.n_obs}")
sc.pp.filter_cells(CT1te, min_counts=1900)
sc.pp.filter_cells(CT1te, max_counts=20000)
#sc.pp.filter_genes(CT1te, min_cells = 8)
CT1te = CT1te[CT1te.obs['pct_counts_mt'] <27, :]
print(f"#spots after filter: {CT1te.n_obs}")


# In[10]:


#CT2 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: CT2")
sns.distplot(CT2te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(CT2te.obs["total_counts"][CT2te.obs["total_counts"] < 3500], kde=False, bins=40, ax=axs[1])
sns.distplot(CT2te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(CT2te.obs["n_genes_by_counts"][CT2te.obs["n_genes_by_counts"] < 2000], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(CT2te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {CT2te.n_obs}")
sc.pp.filter_cells(CT2te, min_counts=2700)
sc.pp.filter_cells(CT2te, max_counts=18000)
#sc.pp.filter_genes(CT2te, min_cells = 10)
CT2te = CT2te[CT2te.obs['pct_counts_mt'] <20, :]
print(f"#spots after filter: {CT2te.n_obs}")


# In[11]:


#CT3 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: CT3")
sns.distplot(CT3te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(CT3te.obs["total_counts"][CT3te.obs["total_counts"] < 1500], kde=False, bins=40, ax=axs[1])
sns.distplot(CT3te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(CT3te.obs["n_genes_by_counts"][CT3te.obs["n_genes_by_counts"] < 2000], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(CT3te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {CT3te.n_obs}")
sc.pp.filter_cells(CT3te, min_counts=500)
sc.pp.filter_cells(CT3te, max_counts=30000)
#sc.pp.filter_genes(CT3te, min_cells = 8)
CT3te = CT3te[CT3te.obs['pct_counts_mt'] <20, :]
print(f"#spots after filter: {CT3te.n_obs}")


# In[12]:


#AD1 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: AD1")
sns.distplot(AD1te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(AD1te.obs["total_counts"][AD1te.obs["total_counts"] < 3000], kde=False, bins=40, ax=axs[1])
sns.distplot(AD1te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(AD1te.obs["n_genes_by_counts"][AD1te.obs["n_genes_by_counts"] < 2500], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(AD1te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {AD1te.n_obs}")
sc.pp.filter_cells(AD1te, min_counts=2000)
sc.pp.filter_cells(AD1te, max_counts=23000)
#sc.pp.filter_genes(AD1te, min_cells = 10)
AD1te = AD1te[AD1te.obs['pct_counts_mt'] <27, :]
print(f"#spots after filter: {AD1te.n_obs}")


# In[13]:


#AD2 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: AD2")
sns.distplot(AD2te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(AD2te.obs["total_counts"][AD2te.obs["total_counts"] < 2500], kde=False, bins=40, ax=axs[1])
sns.distplot(AD2te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(AD2te.obs["n_genes_by_counts"][AD2te.obs["n_genes_by_counts"] < 1800], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(AD2te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {AD2te.n_obs}")
sc.pp.filter_cells(AD2te, min_counts=1900)
sc.pp.filter_cells(AD2te, max_counts=23000)
#sc.pp.filter_genes(AD2te, min_cells = 8)
AD2te = AD2te[AD2te.obs['pct_counts_mt'] <20, :]
print(f"#spots after filter: {AD2te.n_obs}")


# In[14]:


#AD3 QC
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
fig.suptitle(f"Covariates for filtering: AD3")
sns.distplot(AD3te.obs["total_counts"], bins=100,kde=False, ax=axs[0])
sns.distplot(AD3te.obs["total_counts"][AD3te.obs["total_counts"] < 2000], kde=False, bins=40, ax=axs[1])
sns.distplot(AD3te.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(AD3te.obs["n_genes_by_counts"][AD3te.obs["n_genes_by_counts"] < 2000], kde=False, bins=60, ax=axs[3])
sc.pl.scatter(AD3te, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
print(f"#spots before MT filter: {AD3te.n_obs}")
sc.pp.filter_cells(AD3te, min_counts=1300)
sc.pp.filter_cells(AD3te, max_counts=20000)
#sc.pp.filter_genes(AD3te, min_cells = 10)
AD3te = AD3te[AD3te.obs['pct_counts_mt'] <20, :]
print(f"#spots after filter: {AD3te.n_obs}")


# In[15]:


#Combined filtered objects
combined = sc.concat([CT1te,CT2te,CT3te,AD1te,AD2te,AD3te])
sp_te.obs['in_filtered'] = sp_te.obs.index.isin(combined.obs.index)
sp_te = sp_te[sp_te.obs['in_filtered'] == True]


# In[16]:


sc.pp.normalize_total(sp_te, target_sum=1e4,inplace=True)


# In[17]:


#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names+'_spte'] = sp_te[sp_te.obs['sample'] == names]


# In[18]:


#Pulls all genes from adata into a list
all_genes = []
for names in sample:# in range(len(sample)):
    all_genes.append(globals()[names+'_spte'].var.index)

import gzip

repeat_list = []
unique_list = []

#Create list for repeats and TEs
with gzip.open(raw_data_path+'rmsk.txt.gz', 'rt') as file:
    for line in file:
        columns = line.rstrip().split('\t')
        name = columns[10]
        classification = columns[11]
        entry = (name, classification)
        repeat_list.append(entry)  # Add each entry to the repeat_list

        if entry not in unique_list:
            unique_list.append(entry)  # Add unique entries to unique_list


# In[19]:


#Define list of different TE families
DNAte = []
LINEte = []
LTRte = []
SINEte = []
Retroposonte = []
classification_list_te = ['DNA', 'LINE', 'LTR', 'SINE', 'Retroposon']
for te_type in classification_list_te:
    for entry in unique_list:
        name, classification = entry
        if classification == te_type:
            globals()[te_type+'te'].append(name)


# In[20]:


#Check TEs in sample anndata for each family 
concernte = ["AluYh9","L1MB4_5","AluSp","HERV-FC1","AluYc5","THER2","PRIMA4_LTR","LTR77","PB1D11",'LTR5_Hs','LTR5A','LTR5B', 'HERVK-int', 'HERVK11-int', 'HERVK22-int','HERVK3-int','HERVK9-int']
DNAtel = []
LTRtel = []
concerntel = []
LINEtel = []
SINEtel = []
Retroposontel = []
for i in range(len(sample)):
    DNAtemp = []
    for te in DNAte:
        if te in all_genes[i]:
            DNAtemp.append(te)
    DNAtel.append(DNAtemp)
    LTRtemp = []
    for te in LTRte:
        if te in all_genes[i]:
            LTRtemp.append(te)
    LTRtel.append(LTRtemp)
    LINEtemp = []
    for te in LINEte:
        if te in all_genes[i]:
            LINEtemp.append(te)
    LINEtel.append(LINEtemp)
    SINEtemp = []
    for te in SINEte:
        if te in all_genes[i]:
            SINEtemp.append(te)
    SINEtel.append(SINEtemp)
    Retrotemp = []
    for te in Retroposonte:
        if te in all_genes[i]:
            Retrotemp.append(te)
    Retroposontel.append(Retrotemp)
    concerntemp = []
    for te in concernte:
        if te in all_genes[i]:
            concerntemp.append(te)
    concerntel.append(concerntemp)


# In[21]:


#For dotplots pdf generation
#for i in range(len(DNAtel)):
#    n = 30  # Number of elements per sublist
#    DNAtel[i] = [DNAtel[i][j:j + n] for j in range(0, len(DNAtel[i]), n)]
#    
#for i in range(len(LTRtel)):
#    n = 30  # Number of elements per sublist
#    LTRtel[i] = [LTRtel[i][j:j + n] for j in range(0, len(LTRtel[i]), n)]
#
#for i in range(len(LINEtel)):
#    n = 30  # Number of elements per sublist
#    LINEtel[i] = [LINEtel[i][j:j + n] for j in range(0, len(LINEtel[i]), n)]
#
#for i in range(len(SINEtel)):
#    n = 30  # Number of elements per sublist
#    SINEtel[i] = [SINEtel[i][j:j + n] for j in range(0, len(SINEtel[i]), n)]
#
#i=0
#for names in sample:
#    for j in range(len(DNAtel[i])):
#        sc.pl.dotplot(globals()[names+'te'], DNAtel[i][j], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False, swap_axes = True, save='DNA'+str(j)+names+'.pdf')
#    for k in range(len(LTRtel[i])):
#        sc.pl.dotplot(globals()[names+'te'], LTRtel[i][k], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False, swap_axes = True, save='LTR'+str(k)+names+'.pdf')
#    for l in range(len(LINEtel[i])):
#        sc.pl.dotplot(globals()[names+'te'], LINEtel[i][l], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False, swap_axes = True, save='LINE'+str(l)+names+'.pdf')
#    for m in range(len(SINEtel[i])):
#        sc.pl.dotplot(globals()[names+'te'], SINEtel[i][m], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False, swap_axes = True,save='SINE'+str(m)+names+'.pdf')
#    sc.pl.dotplot(globals()[names+'te'], Retroposontel[i], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False,swap_axes = True,save='Retro'+names+'.pdf')
#    sc.pl.dotplot(globals()[names+'te'], concerntel[i], groupby='final_cluster', dot_max=0.7,
#                      dendrogram=True, standard_scale='var', show=False, swap_axes = True, save='Concerned'+names+'.pdf')
#    i+=1


# In[22]:


#Create a list for existing TEs in each sample anndata
i=0
list = [DNAtel, LTRtel, LINEtel, SINEtel, Retroposontel]
for names in sample:
    globals()['all_te_in_'+names] = []
    for telist in list:
        for te in telist[i]:
            globals()['all_te_in_'+names].append(te)
    i+=1


# In[23]:


#Subset each sample anndata for only TEs
i=0
for names in sample:
    globals()['allte'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(globals()['all_te_in_'+names])]
    globals()['DNA'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(DNAtel[i])]
    globals()['LTR'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(LTRtel[i])]
    globals()['LINE'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(LINEtel[i])]
    globals()['SINE'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(SINEtel[i])]
    globals()['Retro'+names] = globals()[names+'te'][:, globals()[names+'te'].var.index.isin(Retroposontel[i])]
    i+=1


# In[24]:


#Prepare empty list for TE counts calculation
allte_nc = []
DNA_nc = []
LTR_nc = []
LINE_nc = []
SINE_nc = []
Retro_nc = []


# In[71]:


for names in sample:
    for i in range(len(globals()['allte'+names])):
        allte_nc.append(globals()['allte'+names].X[i].sum())


# In[72]:


for names in sample:
    for i in range(len(globals()['DNA'+names])):
        DNA_nc.append(globals()['DNA'+names].X[i].sum())


# In[73]:


for names in sample:
    for i in range(len(globals()['LTR'+names])):
        LTR_nc.append(globals()['LTR'+names].X[i].sum())


# In[74]:


for names in sample:
    for i in range(len(globals()['LINE'+names])):
        LINE_nc.append(globals()['LINE'+names].X[i].sum())


# In[75]:


for names in sample:
    for i in range(len(globals()['SINE'+names])):
        SINE_nc.append(globals()['SINE'+names].X[i].sum())


# In[76]:


for names in sample:
    for i in range(len(globals()['Retro'+names])):
        Retro_nc.append(globals()['Retro'+names].X[i].sum())


# In[77]:


#Create normalized count columns in sp_te
te_subset = ['allte','DNA','LTR','LINE','SINE','Retro']
for type in te_subset:
    sp_te.obs[type+'_nc'] = globals()[type+'_nc']


# In[78]:


alltelog10_plus_1= []
DNAlog10_plus_1 = []
LTRlog10_plus_1 = []
LINElog10_plus_1 = []
SINElog10_plus_1 = []
Retrolog10_plus_1 = []
for type in te_subset:
    for i in range(len(globals()[type+'_nc'])):
        globals()[type+'log10_plus_1'].append(math.log10(globals()[type+'_nc'][i]+1))
    sp_te.obs[type] = globals()[type+'log10_plus_1']


# In[79]:


#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names+'_spte'] = sp_te[sp_te.obs['sample'] == names]


# In[80]:


for names in sample:
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["allte"], title=names+'_'+'allte'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_allte_log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["DNA"], title=names+'_'+'DNA'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_DNA_log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["LTR"], title=names+'_'+'LTR'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_LTR_log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["LINE"], title=names+'_'+'LINE'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_LINE_log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["SINE"], title=names+'_'+'SINE'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_SINE_log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["Retro"], title=names+'_'+'Retro'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_Retro_log10_normalized_total_counts.pdf')


# In[81]:


sp_te.write(h5ad_path+'normalized_teobj.h5ad')


# In[82]:


sp_te = sc.read_h5ad('normalized_teobj.h5ad')
sample = sp_te.obs['sample'].unique()
te_subset = ['allte','DNA','LTR','LINE','SINE','Retro']


# In[83]:


for type in te_subset:
    globals()[type] = sp_te.obs[type]
    sp_te.obs = sp_te.obs.drop(columns =[type])
    sp_te.obs = sp_te.obs.drop(columns =[type+'_nc'])


# In[84]:


for type in te_subset:
    globals()['spte_'+type] = sp_te.obs.copy()
    globals()['spte_'+type]['values'] = globals()[type]
    globals()['spte_'+type]['TE_family'] = type


# In[85]:


violin_slice_com = pd.concat([spte_allte,spte_DNA,spte_LTR,spte_LINE,spte_SINE,spte_Retro])


# In[86]:


fig, ax = plt.subplots(figsize=(30,5))
sns.violinplot(data=violin_slice_com, y = violin_slice_com['values'], x= violin_slice_com['sample'], hue='TE_family',ax=ax, inner='box', width=0.5, cut=3, rotation = 45)
ax.set_title('Log10_(Normalized_counts+1)_for_TE_famliies_in_each_sample')
ax.set_xlabel("samples")
ax.set_ylabel("Log10_(Normalized_counts+1)")
#plt.xticks(rotation=45)
plt.show()


# In[87]:


for names in sample:
    globals()[names] = violin_slice_com[violin_slice_com['sample'] == names]


# In[ ]:


for names in sample:
    fig, ax = plt.subplots(figsize=(30,5))
    sns.violinplot(data=globals()[names], y = globals()[names]['values'], x= globals()[names]['Layer'], hue='TE_family',ax=ax, inner='box', width=0.5, cut=3, rotation = 45)
    ax.set_title('Log10_(Normalized_counts+1)_for_TE_famliies_in_each_region_in_'+names)
    ax.set_xlabel("regions")
    ax.set_ylabel("Log10_(Normalized_counts+1)")
    plt.show()


# In[ ]:




