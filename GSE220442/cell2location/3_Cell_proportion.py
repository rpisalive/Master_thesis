#!/usr/bin/env python
# coding: utf-8

# In[5]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# In[8]:


sample = ['CT1','CT2','CT3','AD1','AD2','AD3']
cell_type = ['Astro','Endo','Inh','Exc','Micro','OPC','Oligo']
Oligo = []
Astro = []
Endo = []
Inh = []
Exc = []
Micro = []
OPC = []
Total = []
sp_te = sc.read_h5ad('spatial_te_with_cluster_anno.h5ad')
#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names] = sp_te[sp_te.obs['sample'] == names]


# In[9]:


i = 0
for names in sample:
    Total_temp = 0
    for cell in cell_type:
        globals()[cell+'_temp'] = globals()[names].obs[cell].sum()
        globals()[cell].append(globals()[cell+'_temp'])
        Total_temp += globals()[cell+'_temp']
    Total.append(Total_temp)
    for cell in cell_type:
        globals()[cell+'_percentage'] = globals()[cell][i]/Total[i]
        globals()[cell][i] = globals()[cell+'_percentage']
    i += 1


# In[10]:


for cell in cell_type:
    globals()[cell] = np.array(globals()[cell])


# In[11]:


cell_type_temp = cell_type.copy()
plt.bar(sample, Astro, width=0.5)
plotted = Astro
cell_type_temp.remove(cell_type_temp[0])
for cell in cell_type_temp:
    plt.bar(sample, globals()[cell], bottom=plotted,width=0.5)
    plotted +=globals()[cell]
plt.xlabel("SLICE",fontsize=10)
plt.ylabel("Cell type proportion")
plt.legend(cell_type, loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.title("Cell distribution per slice")
plt.show()


# In[35]:


for names in sample:
    globals()[names] = globals()[names][globals()[names].obs['Layer'] != 'Noise']


# In[36]:


region = []
for names in sample:
    globals()[names+'region'] = globals()[names].obs['Layer'].unique().tolist()
    globals()[names+'region'].sort()


# In[37]:


for names in sample:
    globals()[names+'Total'] = pd.Series()
    for i in range(len(globals()[names+'region'])):
        globals()[names+'Total'][i] = 0
    globals()[names+'Total'].index = globals()[names+'region']   
    for cell in cell_type:
        globals()[names+'_'+cell] = globals()[names].obs.groupby('Layer')[cell].sum()
        globals()[names+'Total']+=globals()[names+'_'+cell]
    for cell in cell_type:
        globals()[names+'_'+cell+'_percentage_per_region'] = globals()[names+'_'+cell]/globals()[names+'Total']


# In[38]:


for names in sample:
    for cell in cell_type:
        globals()[names+'_'+cell+'_percentage_per_region'] = np.array(globals()[names+'_'+cell+'_percentage_per_region'])


# In[39]:


for names in sample:
    cell_type_temp = cell_type.copy()
    plt.bar(globals()[names+'region'], globals()[names+'_Astro_percentage_per_region'],width=0.5)
    plotted = globals()[names+'_Astro_percentage_per_region']
    cell_type_temp.remove(cell_type_temp[0])
    for cell in cell_type_temp:
        plt.bar(globals()[names+'region'], globals()[names+'_'+cell+'_percentage_per_region'], bottom=plotted,width=0.5)
        plotted +=globals()[names+'_'+cell+'_percentage_per_region']
    plt.xticks(range(len(globals()[names+'region'])), globals()[names+'region'], rotation='vertical')
#plt.xlabel("SLICE",fontsize
    plt.ylabel("Cell type proportion")
    plt.legend(cell_type, loc='center left',bbox_to_anchor=(1.0, 0.5))
    plt.title(names+" - Cell distribution per region")
    plt.show()


# In[ ]:





# In[ ]:




