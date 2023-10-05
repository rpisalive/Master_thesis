#!/usr/bin/env python
# coding: utf-8

# In[23]:


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


# In[24]:


# Set paths to data and results used through the document:
sp_data_folder = '/work/project/ext_014/EMTAB10972/snakerangerout/spaceranger/'
rawdata_path = '/work/project/ext_014/EMTAB10972/rawdata/'
c2lpath = '/work/project/ext_014/EMTAB10972/cell2location/'


# In[25]:


def read_and_var_change_name(sample_name, path=sp_data_folder):
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

def select_slide(adata, s, s_col='sample'):
    r""" This function selects the data for one slide from the spatial anndata object.

    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """

    slide = adata[adata.obs[s_col].isin([s]), :]
    s_keys = list(slide.uns['spatial'].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]

    slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}

    return slide

#######################
# Read the list of spatial experiments
sample_data = pd.read_csv(rawdata_path+'EMTAB10972visium.csv.gz')

# Read the data into anndata objects
slides = []
for i in sample_data['sample_name']:
    slides.append(read_and_var_change_name(i, path=sp_data_folder))

# Combine anndata objects together
adata_vis = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    uns_merge="unique",
    batch_categories=sample_data['sample_name'],
    index_unique=None
)
#######################
adata_vis_pp = pd.read_csv(rawdata_path+'metadata.tsv.gz',sep='\t')
adata_vis_pp.index = adata_vis_pp["slice"]+ '_' + adata_vis_pp.index.str[:18]
adata_vis_pp.index.name = 'spot_id'
#Size reduce
adata_vis.obs['in_meta'] = adata_vis.obs.index.isin(adata_vis_pp.index)
adata_vis = adata_vis[adata_vis.obs['in_meta'] == True]


# In[26]:


# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]


# In[27]:


#change .obsm type
adata_vis.obsm['spatial'] = np.frompyfunc(lambda x: x.replace(',',''),1,1)(adata_vis.obsm['spatial']).astype(np.int64)


# In[28]:


## scRNAseq reference (raw counts)
adata_scrna_raw = sc.read_10x_h5(rawdata_path+'cerebral_organoidsc.h5')
celltype_anno = pd.read_csv(rawdata_path+'EMTAB10974_annotation.csv.gz', index_col=0)
overlap_barcode = np.intersect1d(adata_scrna_raw.obs.index.tolist(), celltype_anno.index.tolist())
celltype_anno = celltype_anno.loc[overlap_barcode, :]
adata_scrna_raw = adata_scrna_raw[overlap_barcode, :]
adata_scrna_raw.obs['celltype'] = celltype_anno['annot_ct'].copy()
adata_scrna_raw.obs['orig.ident'] = celltype_anno['orig.ident'].copy()
adata_scrna_raw.obs['batch'] = celltype_anno['batch'].copy()
adata_scrna_raw.obs['genotype'] = celltype_anno['genotype'].copy()
adata_scrna_raw.obs['annot_region'] = celltype_anno['annot_region'].copy()
adata_scrna_raw.var['SYMBOL'] = adata_scrna_raw.var_names
adata_scrna_raw.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
adata_scrna_raw.var_names = adata_scrna_raw.var['ENSEMBL']
adata_scrna_raw.var.drop(columns='ENSEMBL', inplace=True)


# In[29]:


# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata_scrna_raw, inplace=True)
adata_scrna_raw.var['mt'] = [gene.startswith('mt-') for gene in adata_scrna_raw.var['SYMBOL']]
adata_scrna_raw.obs['mt_frac'] = adata_scrna_raw[:, adata_scrna_raw.var['mt'].tolist()].X.sum(1).A.squeeze()/adata_scrna_raw.obs['total_counts']


# In[30]:


from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_scrna_raw, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_scrna_raw = adata_scrna_raw[:, selected].copy()


# In[31]:


# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_scrna_raw,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='celltype',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['batch']
                       )


# In[32]:


# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_scrna_raw)

# view anndata_setup as a sanity check
mod.view_anndata_setup()


# In[33]:


#for regression model training, please remove the # from the mod.train line and mute the 2nd line, it is suggested to submit the script to slurm job for training
#The regression modle can be reused for cell2location training as long as they share the same reference data.
#mod.train(max_epochs=250, use_gpu=False)
mod = cell2location.models.RegressionModel.load(c2lpath+'combinedvis_scregression_model', adata_scrna_raw)


# In[34]:


mod.plot_history(20)


# In[35]:


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_scrna_raw = mod.export_posterior(
    adata_scrna_raw, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
#mod.save('combinedvis_scregression_model', overwrite=True)


# In[36]:


adata_scrna_raw = mod.export_posterior(
    adata_scrna_raw, use_quantiles=True,
    # choose quantiles
    add_to_varm=["q05","q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'use_gpu': False}
)


# In[37]:


# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_scrna_raw.varm.keys():
    inf_aver = adata_scrna_raw.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_scrna_raw.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_scrna_raw.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_scrna_raw.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# In[38]:


#Required for the following shared genes finding
adata_vis.var['gene_ids'] = adata_vis.var_names
adata_vis.var.set_index('SYMBOL', drop=True, inplace=True)


# In[39]:


# find shared genes and subset both anndata and reference signatures
adata_vis.var_names_make_unique()
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")


# In[40]:


# create and train the model, please refer to the tutorial for the selection of N_cells_per_location and detection_alpha
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


# In[41]:


#In case of model training, unmute mod.train and mute the mod loading line
#mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
#          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
#          train_size=1,
#          use_gpu=False,
#         )
mod = cell2location.models.Cell2location.load(c2lpath+'combined_model', adata_vis)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);


# In[42]:


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)


# In[43]:


mod.plot_QC()


# In[44]:


fig = mod.plot_spatial_QC_across_batches()


# In[45]:


# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# select one slide
from cell2location.utils import select_slide
slide = select_slide(adata_vis, 'S3')

# plot in spatial coordinates
with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(slide, cmap='magma',
                  # show first 8 cell types
                  color=['choroid_plexus','dien_NPC','dien_neuron','epi','hind_NPC','hind_neuron','mesen','mg',
                         'nc','pns','retina','tel_NPC','tel_neuron'],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )


# In[48]:


# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = ['choroid_plexus','dien_NPC','dien_neuron','epi','hind_NPC','hind_neuron','mesen']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'S3')

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
        colorbar_position='right')


# In[49]:


# select up to 6 clusters
clust_labels = ['mg','nc','pns','retina','tel_NPC','tel_neuron']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'S3')

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


# In[55]:


# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata_vis, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata_vis, resolution=1.1)

# add region as categorical variable
adata_vis.obs["region_cluster"] = adata_vis.obs["leiden"].astype("category")


# In[56]:


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


# In[57]:


# plot in spatial coordinates
sample = ['S1', 'S2', 'S3']
for names in sample:
    ad = adata_vis[adata_vis.obs['sample'] == names].copy()
    with mpl.rc_context({'axes.facecolor':  'black','figure.figsize': [4.5, 5]}):
        sc.pl.spatial(ad, color=['region_cluster'], size=1.3, img_key='hires',library_id=names,alpha=0.5,title=names)


# In[58]:


adata_vis.write(c2lpath+'saved_h5ad/combined_regional_cluster.h5ad')


# In[ ]:




