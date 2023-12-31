#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os

data_type = 'float32'

# this line forces theano to use the GPU and should go before importing cell2location
#os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=' + data_type + ',force_device=True'
# if using the CPU uncomment this:
os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

# silence scanpy that prints a lot of warnings
import warnings
warnings.filterwarnings('ignore')


# In[16]:


# Set paths to data and results used through the document:
sp_data_folder = '/work/project/ext_014/GSE220442/snakerangerout/spaceranger/'
raw_data_path = '/work/project/ext_014/GSE220442/rawdata/'
c2lpath = '/work/project/ext_014/GSE220442/cell2location/'


# In[9]:


def read_and_var_change_name_sp(sample_name, path=sp_data_folder):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    adata = sc.read_visium(path + str(sample_name),
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.set_index('gene_ids', drop=True, inplace=True)
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata

def read_and_var_change_name_paper(sample_name, path=raw_data_path+'counts_and_images/'):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param path: path to data
    """
    adata = sc.read_visium(path + str(sample_name),
                           count_file='filtered_feature_bc_matrix.h5', load_images=True)
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.set_index('gene_ids', drop=True, inplace=True)
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'

    return adata
#######################
# Read the list of spatial experiments
sample_data = pd.read_csv(raw_data_path+'GSE220442visium.csv.gz')

# Read the data into anndata objects
slides_sp = []
for i in sample_data['sample_name']:
    slides_sp.append(read_and_var_change_name_sp(i, path=sp_data_folder))

# Combine anndata objects together
adata_vis = slides_sp[0].concatenate(
    slides_sp[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)
#######################

# Read the data into anndata objects
slides_pp = []
for i in sample_data['sample_name']:
    slides_pp.append(read_and_var_change_name_paper(i, path=paper_data_folder))

# Combine anndata objects together
adata_vis_pp = slides_pp[0].concatenate(
    slides_pp[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)
#######################
#Size matching
adata_vis.obs['in_ori'] = adata_vis.obs.index.isin(adata_vis_pp.obs.index)
adata_vis = adata_vis[adata_vis.obs['in_ori']==True]


# In[10]:


# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


# In[11]:


#change .obsm type
adata_vis.obsm['spatial'] = np.frompyfunc(lambda x: x.replace(',',''),1,1)(adata_vis.obsm['spatial']).astype(np.int64)


# In[12]:


## scRNAseq reference (raw counts)
adata_scrna_raw = sc.read_10x_h5(raw_data_path+'mtg_nature_pure.h5') #mtg_nature_pure.h5 generated by excuting the R code from the original paper
celltype_anno = pd.read_csv(raw_data_path+'MTG_pure_annotation.csv.gz', index_col=0) #MTG_pure_annotation.csv generated by excuting the R code from the original paper
overlap_barcode = np.intersect1d(adata_scrna_raw.obs.index.tolist(), celltype_anno.index.tolist())
celltype_anno = celltype_anno.loc[overlap_barcode, :]
adata_scran_raw = adata_scrna_raw[overlap_barcode, :]
adata_scrna_raw.obs['celltype'] = celltype_anno['celltype'].copy()
adata_scrna_raw.obs['sample_id'] = celltype_anno['sample_id'].copy()
adata_scrna_raw.obs['donor'] = celltype_anno['donor'].copy()
adata_scrna_raw.var['SYMBOL'] = adata_scrna_raw.var_names
adata_scrna_raw.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
adata_scrna_raw.var_names = adata_scrna_raw.var['ENSEMBL']
adata_scrna_raw.var.drop(columns='ENSEMBL', inplace=True)


# In[13]:


from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_scrna_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_scrna_raw = adata_scrna_raw[:, selected].copy()


# In[14]:


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_scrna_raw,
                        # 10X reaction / sample / batch
                        #batch_key='batch',
                        # cell type, covariate used for constructing signatures
                        labels_key='celltype',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['donor']
                       )


# In[15]:


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_scrna_raw)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


# In[17]:


#mod.train is for model traning, the trained model can be loaded with cell2location.models.RegressionModel.load
#Please refer to cell2location_run.sh for slrum job submission
#The regression modle can be reused for cell2location training as long as they share the same reference data.
mod.train(max_epochs=250, use_gpu=False)
#mod = cell2location.models.RegressionModel.load(c2lpath+'combinedvis_scregression_model', adata_scrna_raw)


# In[19]:


mod.plot_history(20)


# In[20]:


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_scrna_raw = mod.export_posterior(
    adata_scrna_raw, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
mod.save('combinedvis_scregression_model', overwrite=True)


# In[21]:


adata_scrna_raw = mod.export_posterior(
    adata_scrna_raw, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'use_gpu': False}
)


# In[22]:


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_scrna_raw.varm.keys():
    inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_scrna_raw.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# In[23]:


#Required for finding shared genes
adata_vis.var['gene_ids'] = adata_vis.var_names
adata_vis.var.set_index('SYMBOL', drop=True, inplace=True)


# In[24]:


# find shared genes and subset both anndata and reference signatures
adata_vis.var_names_make_unique()
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


# In[25]:


#Please refer to cell2location page for the settings of N_cells_per_location and detection_alpha
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5.5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()


# In[26]:

#mod.train is for model traning, the trained model can be loaded with cell2location.models.Cell2location.load
#Please refer to cell2location_run.sh for slrum job submission
#mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
#          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
#          train_size=1,
#          use_gpu=True,
#         )
mod = cell2location.models.Cell2location.load(c2lpath+'combined_model', adata_vis)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);


# In[27]:


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)
# Save model
mod.save(c2lpath+'combined_model', overwrite=True)

# In[28]:


mod.plot_QC()


# In[29]:


#Change the key in .uns slot to match with plot_spatial_QC_across_batches function, if it is already matching, no need to excute this section
def change_dict_key(d, old_key, new_key, default_value=None):
    d[new_key] = d.pop(old_key, default_value)
change_dict_key(adata_vis.uns['spatial'], '1-1', 'CT1')
change_dict_key(adata_vis.uns['spatial'], '18-64', 'CT2')
change_dict_key(adata_vis.uns['spatial'], '2-5', 'CT3')
change_dict_key(adata_vis.uns['spatial'], '2-3', 'AD1')
change_dict_key(adata_vis.uns['spatial'], '2-8', 'AD2')
change_dict_key(adata_vis.uns['spatial'], 'T4857', 'AD3')


# In[30]:


#Hyperparameters checking
#fig = mod.plot_spatial_QC_across_batches()


# In[31]:


# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'AD3')

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=['Astro','Endo','Exc','Inh','Micro','Oligo','OPC'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )


# In[32]:


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['Astro','Endo','Exc','Inh','Micro','Oligo','OPC']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'AD3')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='dark_background',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )


# In[33]:


# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")


# In[34]:


# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata_vis, min_dist = 0.3, spread = 1)

# show regions in UMAP coordinates
with mpl.rc_context({'axes.facecolor':  'white',
                     'figure.figsize': [8, 8]}):
    sc.pl.umap(adata_vis, color=['region_cluster'], size=30,
               color_map = 'RdPu', ncols = 2, legend_loc='on data',
               legend_fontsize=20)
    sc.pl.umap(adata_vis, color=['sample'], size=30,
               color_map = 'RdPu', ncols = 2,
               legend_fontsize=20)


# In[35]:


# plot in spatial coordinates
sample = ['CT1','CT2','CT3','AD1','AD2','AD3']
for names in sample:
    ad = adata_vis[adata_vis.obs['sample'] == names].copy()
    with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': [4.5, 5]}):
        sc.pl.spatial(ad, color=['region_cluster'], size=1.3, img_key='hires',library_id=names,alpha=0.5,title=names)


# In[62]:


adata_vis.write(c2lpath+'saved_h5ad/combined_regional_cluster.h5ad')


# In[ ]:




