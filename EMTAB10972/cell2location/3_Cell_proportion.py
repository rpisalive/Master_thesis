#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# In[16]:


h5ad_path = '/work/project/ext_014/EMTAB10972/cell2location/saved_h5ad/'


# In[17]:


sample = ['S1','S2','S3']
cell_type = ['choroid_plexus','dien_NPC','dien_neuron','epi','hind_NPC','hind_neuron','mesen','mg'
            , 'nc','pns','retina','tel_NPC','tel_neuron']
choroid_plexus = []
dien_NPC = []
dien_neuron = []
epi = []
hind_NPC = []
hind_neuron = []
mesen = []
mg = []
nc = []
pns = []
retina = []
tel_NPC = []
tel_neuron = []
Total = []
sp_te = sc.read_h5ad(h5ad_path+'spatial_te_with_cluster_anno.h5ad')
#Splitting the whole object into individual tissue object
for names in sample:
    globals()[names] = sp_te[sp_te.obs['sample'] == names]


# In[18]:


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


# In[19]:


for cell in cell_type:
    globals()[cell] = np.array(globals()[cell])


# In[20]:


cell_type_temp = cell_type.copy()
plt.bar(sample, choroid_plexus, color = 'black',width=0.5)
plotted = choroid_plexus
cell_type_temp.remove(cell_type_temp[0])
plt.bar(sample, dien_NPC, color = 'purple',bottom=plotted,width=0.5)
plotted +=dien_NPC
cell_type_temp.remove(cell_type_temp[0])
plt.bar(sample, dien_neuron, color = 'yellow',bottom=plotted,width=0.5)
plotted +=dien_neuron
cell_type_temp.remove(cell_type_temp[0])
for cell in cell_type_temp:
    plt.bar(sample, globals()[cell], bottom=plotted,width=0.5)
    plotted +=globals()[cell]
plt.xlabel("SLICE",fontsize=10)
plt.ylabel("Cell type proportion")
plt.legend(cell_type, loc='center left',bbox_to_anchor=(1.0, 0.5))
plt.title("Cell distribution per tissue slice")
plt.show()


# In[21]:


region = []
for names in sample:
    globals()[names+'region'] = globals()[names].obs['final_cluster'].unique().tolist()
    globals()[names+'region'].sort()


# In[22]:


for names in sample:
    globals()[names+'Total'] = pd.Series()
    for i in range(len(globals()[names+'region'])):
        globals()[names+'Total'][i] = 0
    globals()[names+'Total'].index = globals()[names+'region']   
    for cell in cell_type:
        globals()[names+'_'+cell] = globals()[names].obs.groupby('final_cluster')[cell].sum()
        globals()[names+'Total']+=globals()[names+'_'+cell]
    for cell in cell_type:
        globals()[names+'_'+cell+'_percentage_per_region'] = globals()[names+'_'+cell]/globals()[names+'Total']
#        print(globals()[names+'_'+cell+'_percentage_per_region'])


# In[23]:


for names in sample:
    for cell in cell_type:
        globals()[names+'_'+cell+'_percentage_per_region'] = np.array(globals()[names+'_'+cell+'_percentage_per_region'])


# In[24]:


for names in sample:
    cell_type_temp = cell_type.copy()
    plt.bar(globals()[names+'region'], globals()[names+'_choroid_plexus_percentage_per_region'], color = 'black',width=0.5)
    plotted = globals()[names+'_choroid_plexus_percentage_per_region']
    cell_type_temp.remove(cell_type_temp[0])
    plt.bar(globals()[names+'region'], globals()[names+'_dien_NPC_percentage_per_region'], color = 'purple',bottom=plotted,width=0.5)
    plotted +=globals()[names+'_dien_NPC_percentage_per_region']
    cell_type_temp.remove(cell_type_temp[0])
    plt.bar(globals()[names+'region'], globals()[names+'_dien_neuron_percentage_per_region'], color = 'yellow',bottom=plotted,width=0.5)
    plotted +=globals()[names+'_dien_neuron_percentage_per_region']
    cell_type_temp.remove(cell_type_temp[0])
    for cell in cell_type_temp:
        plt.bar(globals()[names+'region'], globals()[names+'_'+cell+'_percentage_per_region'], bottom=plotted,width=0.5)
        plotted +=globals()[names+'_'+cell+'_percentage_per_region']
    plt.xticks(range(len(globals()[names+'region'])), globals()[names+'region'], rotation=90)
#plt.xlabel("SLICE",fontsize
    plt.ylabel("Cell type proportion", rotation = 90)
    plt.legend(cell_type, loc='center left',bbox_to_anchor=(1.0, 0.5))
    plt.title(names+" - Cell distribution per region")
    plt.show()


# In[ ]:





# In[ ]:




