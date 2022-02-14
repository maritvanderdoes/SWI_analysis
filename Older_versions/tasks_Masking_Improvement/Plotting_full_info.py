#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:57:35 2020

@author: Marit
"""

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# dir of document
# dir='D:\Analysis_HBL1\Results_quant_hist.csv'
dir='D:\Analysis_HBL1\Results_quant_hist_3_5.csv'
# dir='D:\Analysis_HBL1\Results_quant_hist_4_Sorting.csv'

#open resultfile and save it in different variables
result_file = pd.read_csv(dir)
# Values
result_file['norm_data'] = result_file['mean_curved'] - result_file['min_curved']
# result_file['norm_data'] = result_file['mean_curved'] - result_file['min_curved']
# result_file['norm_data'] = result_file['median_curved']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['median_curved_sp']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['min_curved_sp']
# result_file['norm_data'] = result_file['mean_bg']
# result_file['norm_data'] = result_file['max_curved'] - result_file['min_curved']
# result_file['norm_data'] = result_file['mean_curved'] - result_file['mean_surr']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['mean_surr_sp']
# result_file['norm_data'] = result_file['mean_curved_sp'] -result_file['median_bg_sp']


time_vals = np.arange(0,np.max(result_file['Frame'])+1,1)/4

controls = result_file['Position']>=40

result_ctrl = result_file[controls]
result_smpl = result_file[np.invert(controls)]

sum_values = np.zeros([np.max(result_file['Frame'])+1,2])
avg_values = np.zeros([np.max(result_file['Frame'])+1,2])

for k_sel in np.unique(result_ctrl['Position']):
    n_vals = np.sum(result_ctrl['Position']==k_sel)
    sum_values[0:n_vals,0] = sum_values[0:n_vals,0] + result_ctrl[result_ctrl['Position']==k_sel]['norm_data']

for k_sel in np.unique(result_smpl['Position']):
    n_vals = np.sum(result_smpl['Position']==k_sel)
    sum_values[0:n_vals,1] = sum_values[0:n_vals,1] + result_smpl[result_smpl['Position']==k_sel]['norm_data']

avg_values[:,0] = sum_values[:,0]/np.unique(result_ctrl['Position']).shape[0]
avg_values[:,1] = sum_values[:,1]/np.unique(result_smpl['Position']).shape[0]

#%% BG_test = avg_values
rat = BG_test[:,0]-BG_test[:,1]
avg_values[:,1] = avg_values[:,1]+rat

# %%
plt.figure()
plt.plot(time_vals,avg_values)
# plt.ylim([100, 160])
plt.grid()

plt.figure()
plt.plot(time_vals,np.diff(avg_values))
plt.grid()
# %%
position_sel = 40
position_cmp = 4
# Columns of interest
columns = result_ctrl.columns[77:-3]
# Histogram plotting
hist_vals = result_file[result_file['Position']==position_sel][columns]
hist_vals = hist_vals.div(np.sum(hist_vals,axis=1), axis = 'rows')

hist_comp = result_file[result_file['Position']==position_cmp][columns]
hist_comp = hist_comp.div(np.sum(hist_comp,axis=1), axis = 'rows')

ax = sns.heatmap((hist_vals.T+8))
ax.invert_yaxis()

plt.figure()
plt.plot(columns.values.astype(int),(hist_vals.T+8));

plt.figure()
plt.plot(columns.values.astype(int),np.log2(hist_vals.iloc[40]),'x-') 
plt.plot(columns.values.astype(int),np.log2(hist_comp.iloc[40]),'o-') 
plt.xlabel('Intensity values')
plt.ylabel('Log density')
plt.grid()

#%% Summing histograms
hist_ctrl = np.zeros([np.max(result_file['Frame'])+1,columns.shape[0]])
hist_smpl = np.zeros([np.max(result_file['Frame'])+1,columns.shape[0]])

for k_sel in np.unique(result_ctrl['Position']):
    hist_vals = result_ctrl[result_ctrl['Position']==k_sel][columns] 
    hist_vals = hist_vals.div(np.sum(hist_vals,axis=1), axis = 'rows')
    hist_ctrl[0:hist_vals.shape[0],:] = hist_ctrl[0:hist_vals.shape[0],:] + hist_vals

for k_sel in np.unique(result_smpl['Position']):
    hist_vals = result_smpl[result_smpl['Position']==k_sel][columns] 
    hist_vals = hist_vals.div(np.sum(hist_vals,axis=1), axis = 'rows')
    hist_smpl[0:hist_vals.shape[0],:] = hist_smpl[0:hist_vals.shape[0],:] + hist_vals

hist_ctrl = hist_ctrl/np.unique(result_ctrl['Position']).shape[0]
hist_smpl = hist_smpl/np.unique(result_smpl['Position']).shape[0]
# %%
plt.figure()
ax = sns.heatmap((hist_ctrl.T),yticklabels=columns)
ax.invert_yaxis()

plt.figure()
plt.plot(columns.values.astype(int),np.log2(hist_ctrl.T+8));

plt.figure()
ax = sns.heatmap((hist_smpl.T))
ax.invert_yaxis()


plt.figure()
plt.plot(columns.values.astype(int),np.log2(hist_smpl.T+8));

hist_comp = hist_smpl-hist_ctrl
plt.figure()
ax = sns.heatmap(((hist_smpl-hist_ctrl).T))
ax.invert_yaxis()

plt.figure()
plt.plot(columns.values.astype(int),np.log2((hist_smpl-hist_ctrl).T+8));
# %%
fr_sel = 10
plt.figure()
plt.plot(columns.values.astype(int),hist_ctrl[fr_sel,:])
plt.plot(columns.values.astype(int),hist_smpl[fr_sel,:])

plt.figure()
plt.plot(columns.values.astype(int),hist_comp[fr_sel,:])

plt.figure()
plt.plot(columns.values.astype(int),hist_ctrl[fr_sel,:]/np.max(hist_ctrl[fr_sel,:]))
plt.plot(columns.values.astype(int),hist_smpl[fr_sel,:]/np.max(hist_smpl[fr_sel,:]))

plt.figure()
plt.plot(columns.values.astype(int),hist_smpl[fr_sel,:]/np.max(hist_smpl[fr_sel,:])-hist_ctrl[fr_sel,:]/np.max(hist_ctrl[fr_sel,:]))

# plt.figure()
# plt.plot(columns.values.astype(int),hist_ctrl[fr_sel,:]/np.sum(hist_ctrl[fr_sel,:]))
# plt.plot(columns.values.astype(int),hist_smpl[fr_sel,:]/np.sum(hist_smpl[fr_sel,:]))

# plt.figure()
# plt.plot(columns.values.astype(int),hist_smpl[fr_sel,:]/np.sum(hist_smpl[fr_sel,:])-hist_ctrl[fr_sel,:]/np.sum(hist_ctrl[fr_sel,:]))

# %%
quants = result_ctrl.columns[26:77]

# Summing quantiles
quant_ctrl = np.zeros([np.max(result_file['Frame'])+1,quants.shape[0]])
quant_smpl = np.zeros([np.max(result_file['Frame'])+1,quants.shape[0]])

for k_sel in np.unique(result_ctrl['Position']):
    hist_vals = result_ctrl[result_ctrl['Position']==k_sel][quants] 
    # hist_vals = hist_vals.div(np.sum(hist_vals,axis=1), axis = 'rows')
    quant_ctrl[0:hist_vals.shape[0],:] = quant_ctrl[0:hist_vals.shape[0],:] + hist_vals

for k_sel in np.unique(result_smpl['Position']):
    hist_vals = result_smpl[result_smpl['Position']==k_sel][quants] 
    # hist_vals = hist_vals.div(np.sum(hist_vals,axis=1), axis = 'rows')
    quant_smpl[0:hist_vals.shape[0],:] = quant_smpl[0:hist_vals.shape[0],:] + hist_vals

quant_ctrl = quant_ctrl/np.unique(result_ctrl['Position']).shape[0]
quant_smpl = quant_smpl/np.unique(result_smpl['Position']).shape[0]


#%%
quant_vals = np.linspace(0,1,51)
tp_sel = 10
plt.figure()
plt.plot(quant_ctrl[tp_sel,:],quant_vals)
plt.plot(quant_smpl[tp_sel,:],quant_vals)
plt.xlabel('Pixel value')
plt.ylabel('Quantile')
plt.legend(['Control', 'Sample' ])
plt.grid()
plt.title('Timepoint '+str(tp_sel))

plt.figure()
plt.plot(quant_ctrl[tp_sel,:],quant_ctrl[tp_sel,:],'k:')
plt.plot(quant_ctrl[tp_sel,:],quant_smpl[tp_sel,:],'rx-')
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()
plt.title('Timepoint '+str(tp_sel))

plt.figure()
plt.plot(quant_vals,quant_smpl[tp_sel,:]-quant_ctrl[tp_sel,:],'r')
plt.xlabel('Quantile')
plt.ylabel('$\Delta$')
plt.grid()
plt.title('Timepoint '+str(tp_sel))

#%%
tp_sel = 10
plt.plot(quant_smpl[tp_sel,:],'x')
plt.plot(quant_ctrl[tp_sel,:],'x')

plt.figure()
plt.plot(1/np.diff(quant_smpl[tp_sel,:]),'-x')
plt.plot(1/np.diff(quant_ctrl[tp_sel,:]),'-x')

#%%
plt.figure()
ax = sns.heatmap(quant_smpl[0:150,0:-1].T-quant_ctrl[0:150,0:-1].T)
ax.invert_yaxis()


# %%
plt.figure()
plt.plot(quant_ctrl[0:150,0])
plt.plot(quant_smpl[0:150,0])
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[0]))

plt.figure()
plt.plot(quant_ctrl[0:150,25])
plt.plot(quant_smpl[0:150,25])
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[25]))

plt.figure()
plt.plot(quant_ctrl[0:150,-1])
plt.plot(quant_smpl[0:150,-1])
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[-1]))

plt.figure()
plt.plot(quant_ctrl[0:150,-2])
plt.plot(quant_smpl[0:150,-2])
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[-2]))

#%%
qt_sel = 0
plt.figure()
plt.plot(quant_ctrl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'k:')
plt.plot(quant_smpl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'rx')
plt.xlabel('Sample')
plt.ylabel('Control')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))

qt_sel = 12
plt.figure()
plt.plot(quant_ctrl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'k:')
plt.plot(quant_smpl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'rx')
plt.xlabel('Sample')
plt.ylabel('Control')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))

qt_sel = 25
plt.figure()
plt.plot(quant_ctrl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'k:')
plt.plot(quant_smpl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'rx')
plt.xlabel('Sample')
plt.ylabel('Control')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))

qt_sel = -1
plt.figure()
plt.plot(quant_ctrl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'k:')
plt.plot(quant_smpl[0:140,qt_sel],.8*quant_smpl[0:140,qt_sel]+10,'b:')
plt.plot(quant_smpl[0:140,qt_sel],quant_ctrl[0:140,qt_sel],'rx')
plt.xlabel('Sample')
plt.ylabel('Control')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))

# %%
qt_sel = -1
m = .8
n = 28
sg = -1

plt.figure()
plt.plot(sg*(quant_ctrl[0:140,qt_sel]-m*quant_smpl[0:140,qt_sel]-n))
plt.xlabel('Time')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))


plt.figure()
plt.plot(-sg*quant_ctrl[0:140,qt_sel])
plt.plot(-sg*(m*quant_smpl[0:140,qt_sel]+n))
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[qt_sel]))

# %%
qt_sel = -2
m = .8
n = 28
sg = 1

plt.figure()
plt.plot(sg*(quant_ctrl[0:140,qt_sel]-m*quant_smpl[0:140,qt_sel]-n))
plt.xlabel('Time')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))


plt.figure()
plt.plot(-sg*quant_ctrl[0:140,qt_sel])
plt.plot(-sg*(m*quant_smpl[0:140,qt_sel]+n))
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[qt_sel]))

#%%
qt_sel = 30
m = 1
n = 0
sg = -1

plt.figure()
plt.plot(sg*(quant_ctrl[0:140,qt_sel]-m*quant_smpl[0:140,qt_sel]-n))
plt.xlabel('Time')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))


plt.figure()
plt.plot(-sg*quant_ctrl[0:140,qt_sel])
plt.plot(-sg*(m*quant_smpl[0:140,qt_sel]+n))
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[qt_sel]))

# %%
qt_sel = -1
m = .8
n = 0
sg = -1

plt.figure()
plt.plot(sg*(quant_ctrl[0:140,qt_sel]-m*quant_smpl[0:140,qt_sel]-n))
plt.xlabel('Time')
plt.grid()
plt.title('Quantile selected: '+str(quants[qt_sel]))


plt.figure()
plt.plot(-sg*quant_ctrl[0:140,qt_sel])
plt.plot(-sg*(m*quant_smpl[0:140,qt_sel]+n))
plt.xlabel('Time')
plt.grid()
plt.title(str(quants[qt_sel]))

# %%
plt.figure()
ax = sns.heatmap((quant_ctrl[0:140,:-1].T))
ax.invert_yaxis()

plt.figure()
ax = sns.heatmap((quant_smpl[0:140,:-1].T))
ax.invert_yaxis()

# %%
qt_sel = -1
plt.figure()
for k_sel in np.unique(result_ctrl['Position']):
    hist_sel = result_ctrl[result_ctrl['Position']==k_sel][quants[qt_sel]]
    plt.plot(np.arange(0,hist_sel.shape[0]),hist_sel)

plt.plot(np.arange(0,140),quant_ctrl[0:140,qt_sel],'k')

plt.figure()
for k_sel in np.unique(result_smpl['Position']):
    hist_sel = result_smpl[result_smpl['Position']==k_sel][quants[qt_sel]]
    plt.plot(np.arange(0,hist_sel.shape[0]),hist_sel)

plt.plot(np.arange(0,140),quant_smpl[0:140,qt_sel],'k')

plt.figure()
plt.plot(np.arange(0,140),quant_ctrl[0:140,qt_sel])
plt.plot(np.arange(0,140),quant_smpl[0:140,qt_sel])

plt.figure()
plt.plot(quant_ctrl[0:140,qt_sel],quant_smpl[0:140,qt_sel],'x')


# %% Background correlations
feature = '1'
plt.figure()
plt.plot(result_file['mean_bg'],result_file[feature],'x')
plt.plot(result_ctrl['mean_bg'],result_ctrl[feature],'ro')
plt.ylabel(feature)
plt.xlabel('Background')

tp_sel = 10
plt.figure()
plt.plot(result_file[result_file['Frame']==tp_sel]['mean_bg'],result_file[result_file['Frame']==tp_sel][feature],'x')
plt.plot(result_ctrl[result_ctrl['Frame']==tp_sel]['mean_bg'],result_ctrl[result_ctrl['Frame']==tp_sel][feature],'ro')
plt.ylabel(feature)
plt.xlabel('Background')

plt.figure()
plt.plot(result_file[result_file['Position']==4]['mean_bg'],result_file[result_file['Position']==4][feature],'x')
plt.plot(result_ctrl[result_ctrl['Position']==40]['mean_bg'],result_ctrl[result_ctrl['Position']==40][feature],'ro')
plt.ylabel(feature)
plt.xlabel('Background')

# %%
plt.figure()
plt.plot(result_file[result_file['Position']==41]['mean_bg'],result_file[result_file['Position']==4][feature],'x')
plt.plot(result_ctrl[result_ctrl['Position']==40]['mean_bg'],result_ctrl[result_ctrl['Position']==40][feature],'ro')
plt.ylabel(feature)
plt.xlabel('Background')

#%%

#%%
sp_sel = 40
qt_sel = -1
plt.plot((result_file[result_file['Position']==sp_sel]['Frame']-result_file[result_file['Position']==sp_sel]['Frame'].iloc[0])/4,
    result_file[result_file['Position']==sp_sel][quants[qt_sel]])

# %%
sp_sel = 40
qt_sel = -3
for sp_sel in [1,2,4]:
    plt.plot((result_file[result_file['Position']==sp_sel]['Frame']-result_file[result_file['Position']==sp_sel]['Frame'].iloc[0])/4,
        result_file[result_file['Position']==sp_sel][quants[qt_sel]])

plt.figure()
for sp_sel in [40,41,44]:
    plt.plot((result_file[result_file['Position']==sp_sel]['Frame']-result_file[result_file['Position']==sp_sel]['Frame'].iloc[0])/4,
        result_file[result_file['Position']==sp_sel][quants[qt_sel]])
# %%
