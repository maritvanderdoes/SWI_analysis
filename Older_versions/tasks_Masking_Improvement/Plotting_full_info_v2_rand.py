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
import scipy.signal as sig
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf


# dir of document
# dir='D:\Analysis_HBL1\Results_quant_hist.csv'
# dir='D:\Analysis_HBL1\Results_quant_hist_3_5.csv'
dir='D:\Analysis_HBL1\Results_quant_hist_4_Sorting.csv'

#open resultfile and save it in different variables
result_file = pd.read_csv(dir)
# Values
result_file['norm_data'] = result_file['mean_curved'] - result_file['min_curved']
result_file['norm_data'] = result_file['mean_curved'] - result_file['totest']
result_file['norm_data'] = result_file['mean_curved'] - result_file['mean_bg']
# result_file['norm_data'] = (result_file['mean_curved'] - result_file['totest'])/result_file['mean_curved']
# result_file['norm_data'] = result_file['mean_curved'] - result_file['mean_surr']
# result_file['norm_data'] = result_file['0.5']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['median_curved_sp']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['min_curved_sp']
# result_file['norm_data'] = result_file['mean_bg']
# result_file['norm_data'] = result_file['mean_curved']
# result_file['norm_data'] = result_file['max_curved'] - result_file['min_curved']
# result_file['norm_data'] = result_file['mean_curved'] - result_file['mean_surr']
# result_file['norm_data'] = result_file['mean_curved_sp'] - result_file['mean_surr_sp']
# result_file['norm_data'] = result_file['mean_curved_sp'] -result_file['median_bg_sp']


time_vals = np.arange(0,np.max(result_file['Frame'])+1,1)/4

# Creating adjusted times
result_file['Hours'] = result_file['Frame']

for k_v, k_sel in enumerate(np.unique(result_file['Position'])):
    times_vals = result_file[result_file['Position']==k_sel]['Frame'].values
    times_vals = (times_vals - times_vals[0])/4
    result_file.loc[result_file[result_file['Position']==k_sel].index,'Hours'] = times_vals

# result_file = result_file[result_file['Position']>39]
result_file = result_file[result_file['Position']<40]

new_pos = np.random.permutation(np.unique(result_file['Position']))

result_file['Type'] = 0
result_file.loc[result_file[result_file['Position'].isin(new_pos[0:int(new_pos.shape[0]/2)])].index,'Type'] = 1

    

#%%

# Defining the samples
result_file['Type'] = (result_file['Position']<40).astype(int)

# D
result_ctrl = result_file[result_file['Type']==0]
result_smpl = result_file[result_file['Type']==1]

# Positions
positions_ctrl = np.unique(result_ctrl['Position'])
positions_smpl = np.unique(result_smpl['Position'])

# Values for the measurement
values_ctrl = np.zeros([np.max(result_file['Frame'])+1,positions_ctrl.shape[0]])*np.nan
values_smpl = np.zeros([np.max(result_file['Frame'])+1,positions_smpl.shape[0]])*np.nan

# Values for the background
bg_ctrl = np.zeros([np.max(result_file['Frame'])+1,positions_ctrl.shape[0]])*np.nan
bg_smpl = np.zeros([np.max(result_file['Frame'])+1,positions_smpl.shape[0]])*np.nan

avg_values = np.zeros([np.max(result_file['Frame'])+1,2])*np.nan
avg_bg = np.zeros([np.max(result_file['Frame'])+1,2])*np.nan

for k_v, k_sel in enumerate(positions_ctrl):
    n_vals = np.sum(result_ctrl['Position']==k_sel)
    # values_ctrl[0:n_vals,k_v] = sig.medfilt(result_ctrl[result_ctrl['Position']==k_sel]['norm_data'],5)
    values_ctrl[0:n_vals,k_v] = result_ctrl[result_ctrl['Position']==k_sel]['norm_data']

for k_v, k_sel in enumerate(positions_smpl):
    n_vals = np.sum(result_smpl['Position']==k_sel)
    # values_smpl[0:n_vals,k_v] = sig.medfilt(result_smpl[result_smpl['Position']==k_sel]['norm_data'],5)
    values_smpl[0:n_vals,k_v] = result_smpl[result_smpl['Position']==k_sel]['norm_data']

avg_values[:,0] = np.nanmean(values_ctrl,axis=1)
avg_values[:,1] = np.nanmean(values_smpl,axis=1)

values_bg = np.zeros([np.max(result_file['Frame'])+1,positions_ctrl.shape[0]])*np.nan
for k_v, k_sel in enumerate(positions_ctrl):
    n_vals = np.sum(result_ctrl['Position']==k_sel)
    # values_ctrl[0:n_vals,k_v] = sig.medfilt(result_ctrl[result_ctrl['Position']==k_sel]['norm_data'],5)
    values_bg[0:n_vals,k_v] = result_ctrl[result_ctrl['Position']==k_sel]['mean_surr']

avg_values_bg = np.nanmean(values_bg,axis=1)

# %%
plt.figure()
plt.plot(time_vals,avg_values)
# plt.ylim([100, 160])
plt.xlabel('Time')
plt.ylabel('Fluorescence levels')
plt.grid()

plt.figure()
plt.plot(time_vals,np.diff(avg_values))
plt.xlabel('Time')
plt.ylabel('Difference')
plt.grid()

plt.figure()
plt.plot(avg_values[0:100,0],avg_values[0:100,1],'x')
# plt.ylim([100, 160])
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()

#%%

plt.figure()
plt.plot(result_file[result_file['Type']==0]['totest'],result_file[result_file['Type']==0]['mean_curved'],'x')
plt.xlabel('Median')
plt.ylabel('Mean')
plt.grid()

plt.figure()
plt.plot(result_file[result_file['Type']==1]['totest'],result_file[result_file['Type']==1]['mean_curved'],'o')
plt.xlabel('Median')
plt.ylabel('Mean')
plt.grid()

plt.figure()
plt.plot(result_file[result_file['Type']==0]['totest'],result_file[result_file['Type']==0]['mean_curved'],'x')
plt.plot(result_file[result_file['Type']==1]['totest'],result_file[result_file['Type']==1]['mean_curved'],'o')
plt.xlabel('Median')
plt.ylabel('Mean')
plt.grid()

# %%
plt.plot(values_smpl[:,3])
# %% Checking background
tp_sel = 0
result_sel = result_file[result_file['Hours']==tp_sel]
plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_bg'],result_sel[result_sel['Type']==0]['mean_curved'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_bg'],result_sel[result_sel['Type']==1]['mean_curved'],'o')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved'],'o')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved']-result_sel[result_sel['Type']==0]['mean_surr'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved']-result_sel[result_sel['Type']==1]['mean_surr'],'o')
plt.ylabel('Difference')
plt.xlabel('Average Surroundings')

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved']/result_sel[result_sel['Type']==0]['mean_surr'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved']/result_sel[result_sel['Type']==1]['mean_surr'],'o')
plt.ylabel('Ratio')
plt.xlabel('Average Surroundings')

# %%
plt.figure()
for sp_sel in np.arange(0,40,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel['mean_bg'],result_sel['mean_curved'],'o')
    plt.ylabel('Mean Value')
    plt.xlabel('Average Surroundings')

plt.figure()
for sp_sel in np.arange(40,67,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel['mean_bg'],result_sel['mean_curved'],'x')
    plt.ylabel('Mean Value')
    plt.xlabel('Average Surroundings')

# %%
tp_th = 10
plt.figure()
for sp_sel in np.arange(0,40,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel[result_sel['Hours']<tp_th]['mean_bg'],result_sel[result_sel['Hours']<tp_th]['mean_curved'],'o')
    plt.ylabel('Mean Value')
    plt.xlabel('Average Surroundings')

plt.figure()
for sp_sel in np.arange(40,67,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel[result_sel['Hours']<tp_th]['mean_bg'],result_sel[result_sel['Hours']<tp_th]['mean_curved'],'x')
    plt.ylabel('Mean Value')
    plt.xlabel('Average Surroundings')


#%%
tp_th = 20
sp_sel = 41
result_sel = result_file[result_file['Position']==sp_sel]
plt.figure()
plt.plot(result_sel[result_sel['Hours']<tp_th]['mean_bg'],result_sel[result_sel['Hours']<tp_th]['mean_curved'],'-o')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')

plt.figure()
plt.plot(result_sel[result_sel['Hours']<tp_th]['mean_bg'],result_sel[result_sel['Hours']<tp_th]['mean_curved']-result_sel[result_sel['Hours']<tp_th]['mean_bg'],'-o')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')

#%%
tp_sel = 10
result_sel = result_file[result_file['Hours']==tp_sel]

set_sel = 0
X = np.array([result_sel[result_sel['Type']==set_sel]['mean_surr']])
X = np.transpose(X)
X = sm.add_constant(X)
Y = result_sel[result_sel['Type']==set_sel]['mean_curved']

est = sm.OLS(Y, X)
est0 = est.fit()
print('Control')
print(est0.params)

set_sel = 1
X = np.array([result_sel[result_sel['Type']==set_sel]['mean_surr']])
X = np.transpose(X)
X = sm.add_constant(X)
Y = result_sel[result_sel['Type']==set_sel]['mean_curved']

est = sm.OLS(Y, X)
est1 = est.fit()
print('Sample')
print(est1.params)

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved'],'o')
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_surr']*est0.params[1]+est0.params[0],'k:')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_surr']*est1.params[1]+est1.params[0],'k:')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')



plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved']-result_sel[result_sel['Type']==0]['mean_surr'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved']-result_sel[result_sel['Type']==1]['mean_surr'],'o')
plt.ylabel('Difference')
plt.xlabel('Average Surroundings')

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved']/result_sel[result_sel['Type']==0]['mean_surr'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved']/result_sel[result_sel['Type']==1]['mean_surr'],'o')
plt.ylabel('Ratio')
plt.xlabel('Average Surroundings')
# %%


tp_sel = 10
param_mat = np.zeros([150, 2, 2])

for k_tp, tp_sel in enumerate(np.arange(0,150,)/4):
    result_sel = result_file[result_file['Hours']==tp_sel]

    for set_sel in [0,1]:

        X = np.array([result_sel[result_sel['Type']==set_sel]['mean_surr']])
        X = np.transpose(X)
        X = sm.add_constant(X)
        Y = result_sel[result_sel['Type']==set_sel]['mean_curved']

        est = sm.OLS(Y, X)
        est2 = est.fit()

        param_mat[k_tp,set_sel,:] = est2.params[None,None, :]

# %%
plt.figure()
plt.plot(param_mat[:,1,1].squeeze())
plt.plot(param_mat[:,0,1].squeeze())

plt.figure()
plt.plot(param_mat[:,1,0].squeeze())
plt.plot(param_mat[:,0,0].squeeze())

plt.figure()
plt.plot(param_mat[:,1,0].squeeze()-param_mat[:,0,0].squeeze())

# %%
md = smf.mixedlm("mean_curved ~ mean_bg", result_sel, groups=result_sel["Type"])
mdf = md.fit()

print(mdf.summary())
# %%
vcf = {"scenario": "0 + C(scenario)"}                                                         

md = sm.MixedLM.from_formula("mean_curved ~ mean_bg", result_sel, groups=result_sel["Type"])
mdf = md.fit()
print(mdf.summary())

# %%

tp_sel = 20
result_sel = result_file[result_file['Hours']==tp_sel]

X = np.array([result_sel['mean_surr'],result_sel['Type']])
X = np.transpose(X)
X = sm.add_constant(X)
Y = result_sel['mean_curved']

est = sm.OLS(Y, X)
est2 = est.fit()
fit_vals = est2.fittedvalues

print(est2.summary())

plt.figure()
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved'],'x')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],result_sel[result_sel['Type']==1]['mean_curved'],'o')
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],fit_vals[result_sel['Type']==0],'k:')
plt.plot(result_sel[result_sel['Type']==1]['mean_surr'],fit_vals[result_sel['Type']==1],'k:')
plt.ylabel('Mean Value')
plt.xlabel('Average Surroundings')
plt.grid()
#%%
plt.plot(result_sel[result_sel['Type']==0]['mean_surr'],result_sel[result_sel['Type']==0]['mean_curved'],'x')


# %%
n_tps = 170
param_mat = np.zeros([n_tps, 3])
stder_mat = np.zeros([n_tps, 3])
result_file['signal_fitted'] = 0
result_file['signal_bg_removed'] = 0

for k_tp, tp_sel in enumerate(time_vals[0:n_tps]):
    result_sel = result_file[result_file['Hours']==tp_sel]
    # result_sel = result_file[result_file['Hours'].isin(np.arange(tp_sel,tp_sel+1,0.25))]

    X = np.array([result_sel['mean_surr'],result_sel['Type']])
    # X = np.array([result_sel['mean_bg'],result_sel['Type']])
    X = np.transpose(X)
    X = sm.add_constant(X)
    Y = result_sel['mean_curved']
    # Y = result_sel['norm_data']

    est = sm.OLS(Y, X)
    est2 = est.fit()

    param_mat[k_tp,:] = est2.params
    stder_mat[k_tp,:] = est2.bse
    result_file.loc[result_sel.index,'signal_fitted'] = est2.fittedvalues
    result_file.loc[result_sel.index,'signal_bg_removed'] = Y - result_sel['mean_surr']*est2.params[1]


# %%

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,0])
plt.grid()

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,1])
plt.grid()

# plt.figure()
# plt.plot(time_vals[0:n_tps],param_mat[:,2],'b')
# plt.plot(time_vals[0:n_tps],param_mat[:,2]-2*stder_mat[:,2],'b:')
# plt.plot(time_vals[0:n_tps],param_mat[:,2]+2*stder_mat[:,2],'b:')
# plt.grid()
# plt.xlabel('Time (hours)')

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,2],'-')
# plt.plot(time_vals[0:n_tps],np.convolve(param_mat[:,2],np.ones(5),'same')/5,'k')
plt.grid()
plt.xlabel('Time (hours)')
plt.ylabel('Hbl-1 signal')


# %%

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,0])
plt.xlabel('Time (hours)')
plt.ylabel('Common Intercept')
plt.grid()

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,1])
plt.xlabel('Time (hours)')
plt.ylabel('Slope')
plt.grid()

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,2],'b')
plt.plot(time_vals[0:n_tps],param_mat[:,2]-2*stder_mat[:,2],'b:')
plt.plot(time_vals[0:n_tps],param_mat[:,2]+2*stder_mat[:,2],'b:')
plt.grid()
plt.xlabel('Time (hours)')

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,2],'-')
# plt.plot(time_vals[0:n_tps],np.convolve(param_mat[:,2],np.ones(5),'same')/5,'k')
plt.grid()
plt.xlabel('Time (hours)')
plt.ylabel('Hbl-1 signal')

plt.figure()
plt.semilogy(result_file['Hours'],result_file['volume'],'x')

plt.figure()
plt.plot(param_mat[:,0],param_mat[:,1],'x')
plt.xlabel('Intercept')
plt.ylabel('Slope')
plt.grid()

plt.figure()
plt.plot(param_mat[:,0],param_mat[:,2],'x-')
plt.xlabel('Intercept')
plt.ylabel('Signal')
plt.grid()

plt.figure()
plt.plot(param_mat[:,1],param_mat[:,2],'x-')
plt.xlabel('Slope')
plt.ylabel('Signal')
plt.grid()

#%%
plt.figure()
plt.plot(avg_values_bg[0:param_mat.shape[0]],param_mat[:,1],'x')
plt.xlabel('Time (hours)')
plt.ylabel('Slope')
plt.grid()

plt.figure()
plt.plot(avg_values_bg[0:param_mat.shape[0]])
plt.xlabel('Time (hours)')
plt.ylabel('Slope')
plt.grid()

# %%
tp_th = 30
plt.figure()
for sp_sel in np.arange(0,5,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel[result_sel['Hours']<tp_th]['Hours'],result_sel[result_sel['Hours']<tp_th]['signal_fitted'],'-o')
    plt.ylabel('Mean Value')
    plt.xlabel('Time')

plt.figure()
for sp_sel in np.arange(40,42,1):
    result_sel = result_file[result_file['Position']==sp_sel]
    plt.plot(result_sel[result_sel['Hours']<tp_th]['Hours'],result_sel[result_sel['Hours']<tp_th]['signal_fitted'],'-x')
    plt.ylabel('Mean Value')
    plt.xlabel('Time')
# %%


# %%

t_its = 10
n_tps = 160
param_mat = np.zeros([n_tps, 3, t_its, 2])
stder_mat = np.zeros([n_tps, 3, t_its, 2])

dir='D:\Analysis_HBL1\Results_quant_hist_4_Sorting.csv'

#open resultfile and save it in different variables
for typ in [0, 1]:
    for its in np.arange(0,t_its):
        result_file = pd.read_csv(dir)

        time_vals = np.arange(0,np.max(result_file['Frame'])+1,1)/4

        # Creating adjusted times
        result_file['Hours'] = result_file['Frame']

        for k_v, k_sel in enumerate(np.unique(result_file['Position'])):
            times_vals = result_file[result_file['Position']==k_sel]['Frame'].values
            times_vals = (times_vals - times_vals[0])/4
            result_file.loc[result_file[result_file['Position']==k_sel].index,'Hours'] = times_vals

        if typ == 0:
            result_file = result_file[result_file['Position']>39]
        else:
            result_file = result_file[result_file['Position']<40]

        new_pos = np.random.permutation(np.unique(result_file['Position']))

        result_file['Type'] = 0
        result_file.loc[result_file[result_file['Position'].isin(new_pos[0:int(new_pos.shape[0]/2)])].index,'Type'] = 1

        for k_tp, tp_sel in enumerate(time_vals[0:n_tps]):
            result_sel = result_file[result_file['Hours']==tp_sel]
            # result_sel = result_file[result_file['Hours'].isin(np.arange(tp_sel,tp_sel+1,0.25))]

            X = np.array([result_sel['mean_surr'],result_sel['Type']])
            # X = np.array([result_sel['mean_bg'],result_sel['Type']])
            X = np.transpose(X)
            X = sm.add_constant(X)
            Y = result_sel['mean_curved']
            # Y = result_sel['norm_data']

            est = sm.OLS(Y, X)
            est2 = est.fit()

            param_mat[k_tp,:,its,typ] = est2.params
            stder_mat[k_tp,:,its,typ] = est2.bse
# %%
typ_sel = 0
plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,0,:,typ_sel])
plt.grid()

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,1,:,typ_sel])
plt.grid()

# plt.figure()
# plt.plot(time_vals[0:n_tps],param_mat[:,2],'b')
# plt.plot(time_vals[0:n_tps],param_mat[:,2]-2*stder_mat[:,2],'b:')
# plt.plot(time_vals[0:n_tps],param_mat[:,2]+2*stder_mat[:,2],'b:')
# plt.grid()
# plt.xlabel('Time (hours)')

plt.figure()
plt.plot(time_vals[0:n_tps],param_mat[:,2,:,typ_sel],'-')
# plt.plot(time_vals[0:n_tps],np.convolve(param_mat[:,2],np.ones(5),'same')/5,'k')
plt.grid()
plt.xlabel('Time (hours)')
plt.ylabel('Hbl-1 signal')
# %%
plt.figure()
plt.plot(time_vals[0:n_tps],np.mean(param_mat[:,0,:,:],axis=1))
plt.grid()

plt.figure()
plt.plot(time_vals[0:n_tps],np.mean(param_mat[:,1,:,:],axis=1))
plt.grid()
plt.figure()
plt.plot(time_vals[0:n_tps],np.mean(param_mat[:,2,:,:],axis=1))
plt.grid()

# %%
result_sel = result_file[result_file['Hours']==14.5]
sns.scatterplot(x="max_surr",y="mean_curved",data=result_sel,hue="Type")

# %%
sns.lineplot(data=result_file, x="Hours",y="totest", hue = "Type")
# %%

result_file["totest"]=result_file["mean_curved"]-result_file["min_curved"]
# %%
