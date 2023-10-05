#!/usr/bin/env python
# coding: utf-8

# In[112]:


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
import matplotlib as mpl
import math


# In[113]:


#path definition
h5ad_path = '/work/project/ext_014/EMTAB10972/cell2location/saved_h5ad/'
te_path = '/work/project/ext_014/EMTAB10972/snakerangerout/scte/'
raw_data_path = '/work/project/ext_014/EMTAB10972/rawdata/'


# In[114]:


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
sample_data = pd.read_csv(raw_data_path+'EMTAB10972visium.csv.gz') #to read it a list of sample name, creaetd manually
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


# In[115]:


#Max cell type, sum of dien, hind & tel
cols = ['choroid_plexus', 'dien_NPC', 'dien_neuron', 'epi', 'hind_NPC', 'hind_neuron', 'mesen', 'mg', 'nc', 'pns', 'retina', 'tel_NPC', 'tel_neuron']
subset = sp_te.obs[cols]
subset['max_cell_type'] = subset.idxmax(axis=1)
sp_te.obs['max_celltype'] = subset['max_cell_type']
#sp_te.obs.to_csv('sp_te_meta.csv')


# In[116]:


cell_cluster = {'mesen':'Mesenchymal',
                'retina':'Retina',
                'tel_neuron':'Telencephalon',
                'tel_NPC':'Telencephalon',
                'dien_neuron':'Diencephalon–mesencephalon',
                'hind_NPC':'Rhombencephalon',
                'pns':'Neural crest',
                'dien_NPC':'Diencephalon–mesencephalon',
                'epi':'Epithelial',
                'hind_neuron':'Rhombencephalon',
                'choroid_plexus':'Choroid Plexus',
                'mg':'microglia',
                'nc':'Neural crest'}
final_cluster = []
for i in range(len(sp_te.obs)):
    for ct in cell_cluster:
        if sp_te.obs.iloc[i]['max_celltype'] == ct:
                final_cluster.append(cell_cluster[ct])
sp_te.obs['final_cluster'] = final_cluster


# In[117]:


for names in sample:
    ad = sp_te[sp_te.obs['sample'] == names].copy()
    sc.pl.spatial(ad, color=['final_cluster'], size=1.3, img_key='hires',library_id=names,alpha=0.5,title=names)


# In[118]:


sp_te.write(h5ad_path + 'spatial_te_with_cluster_anno.h5ad')


# In[119]:


#Splitting the whole object into individual tissue object
sp_te.var_names_make_unique()
for names in sample:
    globals()[names+'_spte'] = sp_te[sp_te.obs['sample'] == names]


# In[120]:


#fig, axs = plt.subplots(1, 4, figsize=(15, 4))
for names in sample:
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    globals()[names+'_spte'].var["mt"] = globals()[names+'_spte'].var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(globals()[names+'_spte'],qc_vars=["mt"], inplace=True)
    fig.suptitle(f"Covariates for filtering: {names}")
    sns.distplot(globals()[names+'_spte'].obs["total_counts"], bins=100,kde=False, ax=axs[0])
    sns.distplot(globals()[names+'_spte'].obs["total_counts"][globals()[names+'_spte'].obs["total_counts"] < 6000], kde=False, bins=40, ax=axs[1])
    sns.distplot(globals()[names+'_spte'].obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    sns.distplot(globals()[names+'_spte'].obs["n_genes_by_counts"][globals()[names+'_spte'].obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])
sp_te.var["mt"] = sp_te.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(sp_te, qc_vars=["mt"],inplace=True)
sc.pl.violin(sp_te, ['n_genes_by_counts', 'total_counts','pct_counts_mt'], groupby='sample',rotation= 45)


# In[96]:


#QC
sc.pp.filter_cells(sp_te, min_counts=4000)
sc.pp.filter_cells(sp_te, max_counts=40000)
sp_te = sp_te[sp_te.obs["pct_counts_mt"] < 20]
#sc.pp.filter_genes(sp_te, min_cells = 3)
print(f"#cells after MT filter: {sp_te.n_obs}")


# In[97]:


sc.pp.normalize_total(sp_te, target_sum=1e4,inplace=True)


# In[98]:


#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names+'_spte'] = sp_te[sp_te.obs['sample'] == names]


# In[99]:


#Pull all genes from adata into a list
all_genes = []
for names in sample:# in range(len(sample)):
    all_genes.append(globals()[names+'_spte'].var.index)

#Grab repeats and TEs into a list
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


# In[100]:


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


# In[101]:


#Check and categorize TEs in the sample anndata
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


# In[102]:


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


# In[103]:


#Generate a list of all TEs in sample anndata
i=0
list = [DNAtel, LTRtel, LINEtel, SINEtel, Retroposontel]
for names in sample:
    globals()['all_te_in_'+names] = []
    for telist in list:
        for te in telist[i]:
            globals()['all_te_in_'+names].append(te)
    i+=1


# In[104]:


#Establish dataframes for each TE family in each sample
i=0
for names in sample:
    globals()['allte'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(globals()['all_te_in_'+names])]
    globals()['DNA'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(DNAtel[i])]
    globals()['LTR'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(LTRtel[i])]
    globals()['LINE'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(LINEtel[i])]
    globals()['SINE'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(SINEtel[i])]
    globals()['Retro'+names] = globals()[names+'_spte'][:, globals()[names+'_spte'].var.index.isin(Retroposontel[i])]
    i+=1


# In[105]:


#Calculate number of counts grouped by TE families
te_subset = ['allte','DNA','LTR','LINE','SINE','Retro']
allte_nc = []
DNA_nc = []
LTR_nc = []
LINE_nc = []
SINE_nc = []
Retro_nc = []
for names in sample:
    for te_type in te_subset:
        for i in range(len(globals()[te_type+names])):
                           globals()[te_type+'_nc'].append(globals()[te_type+names].X[i].sum())
for type in te_subset:
    sp_te.obs[type+'_nc'] = globals()[type+'_nc']


# In[106]:


#Adding 1 to the normalized total counts
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


# In[73]:


#Generate spatial plots for TE families
for names in sample:
    globals()[names+'_spte'] = sp_te[sp_te.obs['sample'] == names]
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["allte"], title=names+'_'+'allte'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'allte'+'_'+'log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["DNA"], title=names+'_'+'DNA'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'DNA'+'_'+'log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["LTR"], title=names+'_'+'LTR'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'LTR'+'_'+'log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["LINE"], title=names+'_'+'LINE'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'LINE'+'_'+'log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["SINE"], title=names+'_'+'SINE'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'SINE'+'_'+'log10_normalized_total_counts.pdf')
    sc.pl.spatial(globals()[names+'_spte'], img_key="hires",library_id=names, color=["Retro"], title=names+'_'+'Retro'+'_'+'log10_normalized_total_counts',cmap='magma',show=False,save=names+'_'+'Retro'+'_'+'log10_normalized_total_counts.pdf')


# In[107]:


#Save h5ad object
sp_te.write(h5ad_path + 'normalized_teobj.h5ad')


# In[108]:


sp_te = sc.read_h5ad(h5ad_path + 'normalized_teobj.h5ad')
sample = sp_te.obs['sample'].unique()
te_subset = ['allte','DNA','LTR','LINE','SINE','Retro']
#Create individual dataframe for each TE family
for type in te_subset:
    globals()[type] = sp_te.obs[type]
    sp_te.obs = sp_te.obs.drop(columns =[type])
    sp_te.obs = sp_te.obs.drop(columns =[type+'_nc'])
for type in te_subset:
    globals()['spte_'+type] = sp_te.obs.copy()
    globals()['spte_'+type]['values'] = globals()[type]
    globals()['spte_'+type]['TE_family'] = type


# In[109]:


#Combine all TE anndata objects
violin_slice_com = pd.concat([spte_allte,spte_DNA,spte_LTR,spte_LINE,spte_SINE,spte_Retro])


# In[110]:


#Violin plots in all samples
fig, ax = plt.subplots(figsize=(15,5))
sns.violinplot(data=violin_slice_com, y = violin_slice_com['values'], x= violin_slice_com['sample'], hue='TE_family',ax=ax, inner='box', width=0.5, cut=3, rotation = 45)
ax.set_title('Log10_(Normalized_counts+1)_for_TE_famliies_in_each_sample')
ax.set_xlabel("samples")
ax.set_ylabel("Log10_(Normalized_counts+1)")
#plt.xticks(rotation=45)
plt.show()


# In[111]:


#Violin plots in all regions per sample
for names in sample:
    globals()[names] = violin_slice_com[violin_slice_com['sample'] == names]
    fig, ax = plt.subplots(figsize=(30,5))
    sns.violinplot(data=globals()[names], y = globals()[names]['values'], x= globals()[names]['final_cluster'], hue='TE_family',ax=ax, inner='box', width=0.5, cut=3, rotation = 45)
    ax.set_title('Log10_(Normalized_counts+1)_for_TE_famliies_in_each_region_in_'+names)
    ax.set_xlabel("regions")
    ax.set_ylabel("Log10_(Normalized_counts+1)")
    plt.show()


# In[ ]:




