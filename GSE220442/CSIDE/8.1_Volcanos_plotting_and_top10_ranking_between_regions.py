#!/usr/bin/env python
# coding: utf-8

# In[22]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import random


# In[23]:


region_to_plot = 'WM'


# In[24]:


region_dict = {'L1':'Layer 1',
                'L2':'Layer 2',
                'L3':'Layer 3',
                'L4':'Layer 4',
                'L5':'Layer 5',
                'L6':'Layer 6',
                'WM':'White matter'}


# In[25]:


palette_dict = {'Astro':'blue',
                'Endo':'orange',
                'Exc':'green',
                'Inh':'red',
                'Micro':'purple',
                'OPC':'dodgerblue',
                'Oligo':'gold',
                'Insignificant':'lightgrey'
                }


# In[26]:


data_path = '/work/project/ext_014/GSE220442/CSIDE/vol_plot_data/'
region_list = ['L1','L2','L3','L4','L5','L6','WM']


# In[27]:


region = pd.read_csv(data_path+region_to_plot+'/combined_across_samples.csv')
region['nlog10'] = -np.log10(region.p)


# In[28]:


#assigning group color
def map_color(a):
    cell_name , significance = a
    
    cell_type = ['Astro','Endo','Exc','Inh','Micro','OPC','Oligo']
    if significance == 'non_sig':
        return 'Insignificant'
    if significance == 'sig':
        for cell in cell_type:
            if cell_name == cell:
                return cell
region['color'] = region[['cell_name','sig_cat']].apply(map_color, axis = 1)


# In[29]:


#convert nlog10 > 10 
mod_p = []
for pvalues in region.nlog10:
    if pvalues > 10:
        mod_p.append(10)
    else:
        mod_p.append(pvalues)
region['mod_p'] = mod_p


# In[30]:


#Assign shapes for nlog10>10 & >= 10
shape = []
for i in range(len(region)):
    if region.iloc[i].nlog10 <= 10:
        shape.append('nlog10 <= 10')
    else:
        shape.append('nlog10 > 10')
region['shape'] = shape


# In[31]:


cell_type_aval = region['color'].unique().tolist()
shape_aval = region['shape'].unique().tolist()


# In[32]:


palette_list = []
for color in cell_type_aval:
    if color in palette_dict.keys():
        palette_list.append(palette_dict[color])


# In[34]:


plt.figure(figsize = (15,15))

ax = sns.scatterplot(data = region, x = 'log_fc_est', y = 'mod_p', hue = 'color',
                     hue_order = cell_type_aval,
                     palette = palette_list,
                     style = 'shape', style_order = shape_aval,markers = ['o', '^'])
ax.set_title('Differentially expressed genes between '+region_dict[region_to_plot]+' and other regions', size = 20)
#ax.axhline(2.5, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(-0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.set_ylim(-0.5, 10.1)
#texts = []
#for i in range(len(region)):
#    if region.iloc[i].nlog10 != 10:
#        if region.iloc[i].mod_p == 10:
#            texts.append(plt.text(x = region.iloc[i].log_fc_est, y = region.iloc[i].mod_p, s = region.iloc[i].genes,fontsize = 20))
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
plt.legend(loc = 1, bbox_to_anchor = (1.15,1), frameon = False, fontsize = '20')
#plt.savefig('volcano_plots/'+region_to_plot+'_volcano_across_samples.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')


# In[35]:


import gzip

repeat_list = []
unique_list = []
rmsk_path = '/work/project/ext_014/GSE220442/rawdata/'

with gzip.open(rmsk_path+'rmsk.txt.gz', 'rt') as file:
    for line in file:
        columns = line.rstrip().split('\t')
        name = columns[10]
        classification = columns[11]
        entry = (name, classification)
        repeat_list.append(entry)  # Add each entry to the repeat_list

        if entry not in unique_list:
            unique_list.append(entry)  # Add unique entries to unique_list


# In[36]:


TE_list = []
TE_type = []
classification_list_te = ['DNA', 'LINE', 'LTR', 'SINE', 'Retroposon']
for te_type in classification_list_te:
    for entry in unique_list:
        name, classification = entry
        if classification == te_type:
            TE_list.append(name)
            TE_type.append(classification)


# In[37]:


#Create a dictionary for the TEs
TE_dict = dict(zip(TE_list, TE_type))


# In[38]:


TE_check = []
for i in range(len(region)):
    if region.iloc[i].genes in TE_list:
        TE_check.append('TE')
    else:
        TE_check.append('non-TE')
region['TE_check'] = TE_check


# In[39]:


family_check = []
for i in range(len(region)):
    if region.iloc[i].genes in TE_list:
        family_check.append(TE_dict[region.iloc[i].genes])
    else:
        family_check.append('non-TE')
region['TE_family'] = family_check


# In[40]:


plt.figure(figsize = (15,15))
ax = sns.scatterplot(data = region, x = 'log_fc_est', y = 'mod_p', hue = 'color',
                     hue_order = cell_type_aval,
                     palette = palette_list,
                     style = 'shape', style_order = shape_aval,markers = ['o', '^'])
ax.set_title('Differentially expressed genes between '+region_dict[region_to_plot]+' and other regions', size = 20)
#ax.axhline(2.5, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
#ax.axvline(-0.8, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.set_ylim(-0.5, 10.1)
texts = []
for i in range(len(region)):
    if region.iloc[i].sig_cat == 'sig':
        if region.iloc[i].TE_check == 'TE':
            texts.append(plt.text(x = region.iloc[i].log_fc_est, y = region.iloc[i].nlog10, s = region.iloc[i].genes,fontsize = 20))
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
plt.legend(loc = 1, bbox_to_anchor = (1.15,1), frameon = False, fontsize = '20')
plt.savefig('volcano_plots/'+region_to_plot+'_volcano_across_samples.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')


# In[20]:


#List generation for top 10 non-TEs and TEs, ranked by log_fc_est
region_sig = region[region['color']!= 'Insignificant']
region_sig_nonTE = region_sig[region_sig['TE_check'] == 'non-TE']
region_sig_TE = region_sig[region_sig['TE_check'] == 'TE']
region_sig_nonTE = region_sig_nonTE.sort_values(by = 'log_fc_est',ascending=False)
region_sig_TE = region_sig_TE.sort_values(by = 'log_fc_est',ascending=False)


# In[21]:


region_sig_nonTE.to_csv('vol_plot_data/'+region_to_plot+'/'+region_to_plot+'_nonTE_rank.csv',index = False)
region_sig_TE.to_csv('vol_plot_data/'+region_to_plot+'/'+region_to_plot+'_TE_rank.csv',index = False)


# In[ ]:




