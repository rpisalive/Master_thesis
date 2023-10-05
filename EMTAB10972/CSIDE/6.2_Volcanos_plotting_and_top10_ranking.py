#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import random


# In[298]:


region_to_plot = 'nc'


# In[299]:


region_dict = {'m':'Mesenchymal',
                'r':'Rhombencephalon',
                't':'Telencephalon',
                'dm':'Diencephalon–mesencephalon',
                'e':'Epithelial',
                'cp':'Choroid Plexus',
                'nc':'Neural crest'}


# In[300]:


palette_dict = {'mesen':'blue',
                'retina':'orange',
                'tel_neuron':'green',
                'tel_NPC':'red',
                'dien_neuron':'purple',
                'hind_NPC':'dodgerblue',
                'pns':'gold',
                'dien_NPC':'lime',
                'epi':'rosybrown',
                'hind_neuron':'darkolivegreen',
                'choroid_plexus':'navy',
                'Insignificant':'lightgrey',
                'nc':'violet'}


# In[301]:


data_path = '/work/project/ext_014/EMTAB10972/CSIDE/vol_plot_data/'
region_list = ['r','nc','dm','t','cp','m','e']


# In[302]:


region = pd.read_csv(data_path+region_to_plot+'/'+region_to_plot+'_combined_im.csv')
region['nlog10'] = -np.log10(region.p)


# In[303]:


#assigning group color
def map_color(a):
    cell_name , significance = a
    
    cell_type = ['mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','mg','retina']
    if significance == 'non_sig':
        return 'Insignificant'
    if significance == 'sig':
        for cell in cell_type:
            if cell_name == cell:
                return cell
region['color'] = region[['cell_name','sig_cat']].apply(map_color, axis = 1)


# In[304]:


#convert nlog10 > 10 
mod_p = []
for pvalues in region.nlog10:
    if pvalues > 10:
        mod_p.append(10)
    else:
        mod_p.append(pvalues)
region['mod_p'] = mod_p


# In[305]:


#Assign shapes for nlog10>10 & >= 10
shape = []
for i in range(len(region)):
    if region.iloc[i].nlog10 <= 10:
        shape.append('nlog10 <= 10')
    else:
        shape.append('nlog10 > 10')
region['shape'] = shape


# In[306]:


cell_type_aval = region['color'].unique().tolist()
shape_aval = region['shape'].unique().tolist()


# In[307]:


palette_list = []
for color in cell_type_aval:
    if color in palette_dict.keys():
        palette_list.append(palette_dict[color])


# In[308]:


plt.figure(figsize = (15,15))

ax = sns.scatterplot(data = region, x = 'log_fc_est', y = 'mod_p', hue = 'color',
                     hue_order = cell_type_aval,
                     palette = palette_list,
                     style = 'shape', style_order = shape_aval,markers = ['o', '^'])
ax.set_title('Differentially expressed genes between '+region_dict[region_to_plot]+' and other regions', size = 20)
#ax.axhline(2.5, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(-0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.set_ylim(-0.5, 10.5)
texts = []
#for i in range(len(region)):
#    if region.iloc[i].nlog10 != 10:
#     if region.iloc[i].sig_cat == 'sig':
#        texts.append(plt.text(x = region.iloc[i].log_fc_est, y = region.iloc[i].mod_p, s = region.iloc[i].genes,fontsize = 7))
#adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'))
for axis in ['bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(width = 2)
plt.xticks(size = 20, weight = 'bold')
plt.yticks(size = 20, weight = 'bold')
plt.xlabel('$log_{2}$ fold change', size = 20)
plt.ylabel('-$log_{10}$ CSIDE p-value', size = 20)
plt.legend(loc = 1, bbox_to_anchor = (1,1), frameon = False,  fontsize = '20')
#plt.savefig('volcano_plots/'+region_to_plot+'_volcano_scTEref_iw.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')


# In[159]:


import gzip

repeat_list = []
unique_list = []
rmsk_path = '/work/project/ext_014/EMTAB10972/rawdata/'

with gzip.open(rmsk_path+'rmsk.txt.gz', 'rt') as file:
    for line in file:
        columns = line.rstrip().split('\t')
        name = columns[10]
        classification = columns[11]
        entry = (name, classification)
        repeat_list.append(entry)  # Add each entry to the repeat_list

        if entry not in unique_list:
            unique_list.append(entry)  # Add unique entries to unique_list


# In[309]:


TE_list = []
TE_type = []
classification_list_te = ['DNA', 'LINE', 'LTR', 'SINE', 'Retroposon']
for te_type in classification_list_te:
    for entry in unique_list:
        name, classification = entry
        if classification == te_type:
            TE_list.append(name)
            TE_type.append(classification)


# In[310]:


#Create a dictionary for the TEs
TE_dict = dict(zip(TE_list, TE_type))


# In[311]:


TE_check = []
for i in range(len(region)):
    if region.iloc[i].genes in TE_list:
        TE_check.append('TE')
    else:
        TE_check.append('non-TE')
region['TE_check'] = TE_check


# In[312]:


family_check = []
for i in range(len(region)):
    if region.iloc[i].genes in TE_list:
        family_check.append(TE_dict[region.iloc[i].genes])
    else:
        family_check.append('non-TE')
region['TE_family'] = family_check


# In[313]:


plt.figure(figsize = (15,15))
ax = sns.scatterplot(data = region, x = 'log_fc_est', y = 'mod_p', hue = 'color',
                     hue_order = cell_type_aval,
                     palette = palette_list,
                     style = 'shape', style_order = shape_aval,markers = ['o', '^'])
ax.set_title('Differentially expressed genes between '+region_dict[region_to_plot]+' and other regions', size = 20)
#ax.axhline(2.5, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(-0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.set_ylim(-0.5, 10.5)
texts = []
for i in range(len(region)):
    if region.iloc[i].TE_check == 'TE':
        if region.iloc[i].sig_cat == 'sig':
            texts.append(plt.text(x = region.iloc[i].log_fc_est, y = region.iloc[i].nlog10, s = region.iloc[i].genes,fontsize = 25))
adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'))
for axis in ['bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(width = 2)
plt.xticks(size = 20, weight = 'bold')
plt.yticks(size = 20, weight = 'bold')
plt.xlabel('$log_{2}$ fold change', size = 20)
plt.ylabel('-$log_{10}$ CSIDE p-value', size = 20)
plt.legend(loc = 1, bbox_to_anchor = (1,1), frameon = False,   fontsize = '20')
plt.savefig('volcano_plots/'+region_to_plot+'_volcano_im.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')


# In[688]:


#List generation for top 10 non-TEs and TEs, ranked by log_fc_est
region_sig = region[region['color']!= 'Insignificant']
region_sig_nonTE = region_sig[region_sig['TE_check'] == 'non-TE']
region_sig_TE = region_sig[region_sig['TE_check'] == 'TE']
region_sig_nonTE = region_sig_nonTE.sort_values(by = 'log_fc_est',ascending=False)
region_sig_TE = region_sig_TE.sort_values(by = 'log_fc_est',ascending=False)


# In[689]:


region_sig_nonTE.to_csv('vol_plot_data/'+region_to_plot+'/'+region_to_plot+'_nonTE_rank.csv',index = False)
region_sig_TE.to_csv('vol_plot_data/'+region_to_plot+'/'+region_to_plot+'_TE_rank.csv',index = False)


# In[ ]:




