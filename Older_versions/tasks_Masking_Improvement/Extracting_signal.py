#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:57:35 2020

@author: Marit
"""

#%%
import numpy as np
from numpy.lib.function_base import diff
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as sig
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf

# dir='D:\Analysis_HBL1\Results_quant_hist.csv'
# dir='D:\Analysis_HBL1\Results_quant_hist_3_5.csv'
dir='D:\Analysis_HBL1\Results_quant_hist_4_Sorting.csv'

result_file = pd.read_csv(dir)

# Timepoints
time_vals = np.arange(0,np.max(result_file['Frame'])+1,1)/4

# Creating adjusted times
result_file['Time'] = result_file['Frame']

# Adjusting for timepoints
for k_v, k_sel in enumerate(np.unique(result_file['Position'])):
    times_vals = result_file[result_file['Position']==k_sel]['Frame'].values
    times_vals = (times_vals - times_vals[0])/4
    result_file.loc[result_file[result_file['Position']==k_sel].index,'Time'] = times_vals

# Defining the samples
result_file['Type'] = (result_file['Position']<40).astype(int)


# %%
n_tps = 173

# Creating parameters dataframe
param_df = pd.DataFrame(
    index   = np.arange(0, n_tps),
    columns = ['Time',  'Internal', 'External_Factor', 'Signal', 
                        'Internal SE', 'External_Factor SE', 'Signal SE',
                        'Internal pval', 'External_Factor pval', 'Signal pval'])
param_df['Time'] = result_file['Time']

# Creating the extracted values
result_file[['signal_fitted','signal_bg_removed','bacterial_contribution','noise_res']] = 0

for k_tp, tp_sel in enumerate(param_df['Time'][0:n_tps]):
    # Selection of the time
    # result_sel = result_file[result_file['Time']==tp_sel]
    result_sel = result_file[result_file['Time'].isin(np.arange(tp_sel,tp_sel+.75,0.25))]

    # Creating the variables
    X = result_sel[['mean_surr','Type']]
    X = sm.add_constant(X)
    Y = result_sel['mean_curved']

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    # Saving the parameters
    param_df.loc[param_df[param_df['Time']==tp_sel].index,param_df.columns[1:]] = \
        np.concatenate([est2.params,est2.bse,est2.pvalues])

    # Extracting the features
    result_file.loc[result_sel.index,'signal_fitted'] = est2.fittedvalues
    result_file.loc[result_sel.index,'signal_bg_removed'] = Y - result_sel['mean_surr']*est2.params[1] - est2.params[0]
    result_file.loc[result_sel.index,'bacterial_contribution'] = result_sel['mean_surr']*est2.params[1]
    result_file.loc[result_sel.index,'noise_res'] = Y - est2.fittedvalues

#%%
sns.lineplot(data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Type', legend=True)
# plt.plot(param_df['Time'],param_df['Signal']*(param_df['Signal pval']<.05),'k:')
# plt.legend(['Control','hbl-1 mean','computed signal'])
plt.grid()

#%%
sns.lineplot(data = result_file, x = 'Time', y = 'noise_res', hue = 'Type')
plt.grid()


#%%
# df_m = result_file[result_file['Type']==1][['Position', 'Time', 'signal_bg_removed']]
df_m = result_file[['Position', 'Time', 'signal_bg_removed']]
# df_m = df_m.groupby(['Position', 'Time']).mean()
df_m = df_m.groupby(['Time', 'Position']).mean()
df_m = df_m.unstack(level=0)
df_m.columns = df_m.columns.droplevel()
indices = (-df_m.mean(axis = 1)).argsort()
df_m = df_m.iloc[indices]

f, axs = plt.subplots(2,1,figsize=(15,15))
sns.lineplot(ax = axs[0], data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Type', legend=False)
sns.heatmap(ax = axs[1], data=df_m, robust=True, cbar = False)
# axs.margins(x=0)

# #%%
# signal_bg_removed = result_file[result_file['Type']==1].groupby('Time', as_index=False)['signal_bg_removed'].mean()

# plt.plot(signal_bg_removed['Time'], signal_bg_removed['signal_bg_removed'])
# #%%
# plt.plot(signal_bg_removed['signal_bg_removed'],param_df['Signal'],'x')

# #%%
# sns.lineplot(data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Position')

# #%%
# sns.lineplot(data = result_file, x = 'Time', y = 'signal_fitted', hue = 'Type')
# plt.grid()
# #%%
# sns.lineplot(data = result_file, x = 'Time', y = 'mean_curved', hue = 'Type')
# plt.grid()
#%%
plt.plot(param_df['Time'],param_df['Signal pval']<.05,'k:')
plt.legend(['Control','hbl-1 mean','computed signal'])
plt.xlabel('Time')
plt.ylabel('Significance')
plt.grid()
#%%
sns.lineplot(data = result_file[result_file['Type']==1], x = 'Time', y = 'signal_bg_removed', hue = 'Position', palette=sns.color_palette("husl",pd.unique(result_file[result_file['Type']==1]['Position']).shape[0]), legend = False)
plt.plot(param_df['Time'],param_df['Signal']*(param_df['Signal pval']<.05),'k:')
plt.grid()

#%%
s_sel = 19
plt.plot(result_file[result_file['Position']==s_sel]['Time'], result_file[result_file['Position']==s_sel]['signal_bg_removed'])
plt.plot(param_df['Time'], param_df['Signal'],'k:')
plt.grid()


#%%
# %%
# dir='D:\Analysis_HBL1\Results_quant_hist_4_Sorting.csv'

# results = pd.read_csv(dir)

# dir = 'D:\Analysis_HBL1\Results_parallel_new.csv'

# result_mip = pd.read_csv(dir)

# result_file = results.merge(result_mip[['Position','Frame','mean_mip'] ], on=['Position','Frame'], how = 'inner')


# # Timepoints
# time_vals = np.arange(0,np.max(result_file['Frame'])+1,1)/4

# # Creating adjusted times
# result_file['Time'] = result_file['Frame']

# for k_v, k_sel in enumerate(np.unique(result_file['Position'])):
#     times_vals = result_file[result_file['Position']==k_sel]['Frame'].values
#     times_vals = (times_vals - times_vals[0])/4
#     result_file.loc[result_file[result_file['Position']==k_sel].index,'Hours'] = times_vals

# # Defining the samples
# result_file['Type'] = (result_file['Position']<40).astype(int)

# # %%

#%%
sns.lineplot(data = result_file, x = 'Time', y = 'bacterial_contribution', hue = 'Type')
plt.grid()


bac = result_file.groupby(['Type','Time'], as_index=False)['bacterial_contribution'].mean()
bac = bac[bac['Time']<43]

plt.figure()
plt.semilogy(bac[bac['Type']==0]['bacterial_contribution'].values)
plt.semilogy(bac[bac['Type']==1]['bacterial_contribution'].values)


plt.figure()
plt.semilogy(bac[bac['Type']==0]['bacterial_contribution'].values - bac[bac['Type']==1]['bacterial_contribution'].values)


#%%
# Subset for volume
volume_file = result_file[['Position','Time','volume']]

s_sel = 4
# Extract ordered times
sample_file = volume_file[volume_file['Position']==s_sel].reset_index()
time_idx = sample_file['Time'].argsort()
volume_file = sample_file.loc[time_idx]

# Compute volume
volume = volume_file[volume_file['Position']==s_sel]['volume']
log_volume = np.log2(volume)

plt.plot(log_volume)

# Computing differences
diff_kern = [-1,2,-1]
# diff_kern = [-1,8,0,-8,1]
diff_log_volume = np.convolve(log_volume,diff_kern,'same')

# diff_log_volume[diff_log_volume**2>4] = 0

plt.figure()
plt.plot(diff_log_volume**2>8)

#%%
med_log_volume = sig.medfilt(log_volume,5)

plt.figure()
# plt.plot(log_volume)
plt.plot(med_log_volume)

#%%
avg_kern = np.ones(3)
avg_log_volume = np.convolve(med_log_volume,avg_kern,'same')/np.sum(avg_kern)

plt.figure()
plt.plot(avg_log_volume)


# %%
med_filt_size = 5
avg_kern = np.ones(1)
diff_kern = [-1,0,1]
diff_kern = [-1,8,0,-8,1]

med_log_volume = sig.medfilt(log_volume,med_filt_size)

plt.figure()
# plt.plot(log_volume)
plt.plot(med_log_volume)

avg_log_volume = np.convolve(med_log_volume,avg_kern,'same')/np.sum(avg_kern)

plt.figure()
plt.plot(avg_log_volume)

# Computing differences
diff_log_volume = np.convolve(avg_log_volume,diff_kern,'same')

diff_log_volume[diff_log_volume**2>4] = 0

plt.figure()
plt.plot(diff_log_volume)
# %%
# Removing outliers by computing differences
outlier_th = 1
med_filt_size = 5
# Computing differences
diff_kern = [1,0,-1]
# diff_kern = [-1,8,0,-8,1]
diff_median = np.convolve(log_volume,diff_kern,'same')
median_rep = sig.medfilt(log_volume,med_filt_size)


# Extract extremes
outliers = diff_median**2>outlier_th**2
log_volume[outliers] = median_rep[outliers]

plt.figure()
plt.plot(log_volume)

# %%
avg_kern = np.ones(5)
mlt_kern = np.ones(5)
# diff_kern = [1,0,-1]
diff_kern = [-1,8,0,-8,1]
adj = 0.1



avg_log_volume = np.convolve(log_volume,avg_kern,'same')/np.sum(avg_kern)

plt.figure()
plt.plot(avg_log_volume)
plt.ylim([15, 21])

# Computing differences
diff_log_volume = np.convolve(avg_log_volume,diff_kern,'same')

# Extracting average growth rate
avg_growth = np.mean(diff_log_volume[(diff_log_volume**2<4)&(diff_log_volume>0)])

diff_log_volume[diff_log_volume**2>4] = avg_growth


molts = np.convolve(diff_log_volume<avg_growth*adj,mlt_kern,'same')/np.sum(mlt_kern)>.5


plt.figure()
plt.plot(diff_log_volume)
plt.plot(molts)
plt.grid()

#%%
# Order timepoints (if not done before)
result_file = result_file.sort_values(by = ['Position','Time'])

s_sel = 1 # 25,42
log_volume = np.log2(result_file[result_file['Position']==s_sel]['volume'])
time = result_file[result_file['Position']==s_sel]['Time']

plt.figure
plt.plot(time,log_volume)
# plt.xlim([2,5])



# % FULL ALGORITHM

# PARAMETERS
# Threshold for molt
molt_th = .4
# Threshold for outliers
outlier_th = .5
# Median filter to replace values
med_filt_size = 5
med_filt_size_2 = 9
# Mean filter to remove noise after outlier removal
avg_kern = np.ones(7)
# Mean filter to fix the derivative
avg_kern_diff = np.ones(7)
# Fixing molt filter
mlt_kern = np.ones(5)
# Finite difference scheme for outliers (First order)
outl_kern = [-1,2,-1]
# Finite difference scheme for derivative (Second order)
# diff_kern = [-1,8,0,-8,1]
# diff_kern = [-1/12,2/3,0,-2/3,1/12]
# diff_kern = [1/60,-3/20,3/4,0,-3/4,3/20,-1/60]
diff_kern = [3,1,-4]

# REMOVING OUTLIERS
# Computing how far the values are by finite differences
diff_median = np.convolve(log_volume,outl_kern,'same')
# Values to replace
median_rep = sig.medfilt(log_volume,med_filt_size)
# Detect outliers
outliers = diff_median**2>outlier_th**2
# Remove outliers
log_volume[outliers] = median_rep[outliers]
# Second median
log_volume_2 = sig.medfilt(log_volume,med_filt_size_2)

plt.figure()
plt.plot(time,log_volume)
# plt.xlim([2,5])


# COMPUTING DERIVATIVES
# Removing noise
avg_log_volume = np.convolve(log_volume_2,avg_kern,'same')/np.sum(avg_kern)
# Computing differences
diff_log_volume = np.convolve(avg_log_volume,diff_kern,'same')
# Smoothing
diff_log_volume = np.convolve(diff_log_volume,avg_kern_diff,'same')
# Computing average growth
avg_growth = np.mean(diff_log_volume[(diff_log_volume**2<4)&(diff_log_volume>0)])
# Removing the extremes
diff_log_volume[diff_log_volume**2>4] = avg_growth


# COMPUTING THE MOLTS
molts_pre = diff_log_volume<avg_growth*molt_th
molts = np.convolve(diff_log_volume<avg_growth*molt_th,mlt_kern,'same')/np.sum(mlt_kern)>.5

# Plotting results
plt.figure()
plt.plot(avg_log_volume)

# Plotting results
plt.figure()
plt.plot(outliers)

plt.figure()
plt.plot(diff_log_volume)
plt.plot(molts_pre*avg_growth)
plt.plot(molts*avg_growth)
plt.grid()
plt.ylim([-avg_growth*2, avg_growth*2])


plt.figure(figsize=(8,15))
plt.plot(time, log_volume)
plt.plot(time[molts_pre],log_volume[molts_pre],'x')
plt.plot(time[molts],log_volume[molts],'x')

#%%
from scipy.signal import savgol_filter
yhat = savgol_filter(log_volume, 11, 3) # window size 51, polynomial order 3

plt.figure(figsize=(15,8))
plt.plot(time, yhat)

plt.figure(figsize=(15,8))
plt.plot(time[1:], np.diff(yhat))
plt.grid()



#%%
molts, _ = molt_detection(log_volume)
molts_df, n_molts = molt_analysis(molts, time)

plt.figure(figsize=(15,8))
plt.plot(time, log_volume)
plt.plot(time[molts],log_volume[molts],'x')
# plt.plot(time,log_volume*molts,'x')



# %%

def molt_detection(log_volume, outlier_th = 0.5, molt_th = 0.1, 
                   median_filt = 5, noise_filt = 5, der_filt = 3, molt_filt = 3):
    # Molt detection v1.1
    # Note: Takes the last molt if molt_th>0 
    # (so it slows down the growth, but it does not decreases)
    #
    # Requires:
    #
    # Numpy, Scipy

    # KERNELS
    # Mean filter to remove noise after outlier removal
    avg_kern = np.ones(noise_filt)
    # Mean filter to remove noise after outlier removal
    avg_kern_diff = np.ones(der_filt)
    # Fixing molt filter
    mlt_kern = np.ones(molt_filt)
    # Finite difference scheme for outliers (First order, central, Second Der)
    outl_kern = [-1,2,-1]
    # # Finite difference scheme for derivative (Second order, central)
    # diff_kern = [-1,8,0,-8,1]    
    # # Finite difference scheme for derivative (Third order, central)
    # diff_kern = [1/60,-3/20,3/4,0,-3/4,3/20,-1/60]

    # Finite difference scheme for derivative (Second order, backwards)
    diff_kern = [3,1,-4]


    # REMOVING OUTLIERS
    # Computing how far the values are by finite differences
    diff_median = np.convolve(log_volume,outl_kern,'same')
    # Values to replace
    median_rep = sig.medfilt(log_volume,median_filt)
    # Detect outliers
    outliers = diff_median**2>outlier_th**2
    # Fix extremes
    outliers[[0,-1]] = False
    # Remove outliers
    log_volume[outliers] = median_rep[outliers]

    # COMPUTING DERIVATIVES
    # Removing noise
    avg_log_volume = np.convolve(log_volume,avg_kern,'same')/np.sum(avg_kern)
    # Computing differences
    diff_log_volume = np.convolve(avg_log_volume,diff_kern,'same')
    # Smoothing
    diff_log_volume = np.convolve(diff_log_volume,avg_kern_diff,'same')/np.sum(avg_kern_diff)
    # Computing average growth
    avg_growth = np.mean(diff_log_volume[(diff_log_volume**2<4)&(diff_log_volume>0)])
    # Removing the extremes
    diff_log_volume[diff_log_volume**2>4] = avg_growth

    # COMPUTING THE MOLTS
    molts = np.convolve(diff_log_volume<avg_growth*molt_th,mlt_kern,'same')/np.sum(mlt_kern)>.5

    return molts, (outliers, log_volume, diff_log_volume)

def molt_analysis(molts, time = None, bwd = 0, fwd = 0, time_th = 0, T_correct = 0):
    # Compute the indices
    idx = np.arange(0,molts.shape[0])

    # Check if time is given, otherwise use index
    if np.sum(time == None):
        time = idx

    # Remove possible molts at earlier timepoints
    molts[time<time_th] = 0

    # Finding the changes
    changes = np.diff(molts.astype(int), prepend=0)

    # Computing molt entry and exit times
    molt_entry = time[idx[changes ==  1]+bwd]
    molt_exit  = time[idx[changes == -1]+fwd]

    n_molts = molt_entry.shape[0]

    # Checking last molt exit
    if molt_exit.shape[0]<molt_entry.shape[0]:
        molt_exit = np.append(molt_exit,time[-1])

    # Correcting gaps if necessary
    correction_metric = (molt_entry[1:]-molt_exit[:-1])>T_correct
    idx_correction = np.append(True, correction_metric)
    molt_entry = molt_entry[idx_correction]
    molt_exit  = molt_exit[idx_correction]

    # Computing the molt duration
    molt_duration = molt_exit-molt_entry


    # Creating pandas dataframe column names
    str_molt_entry     = ['Molt_entry_'+str(x+1) for x in range(molt_entry.shape[0])]
    str_molt_exit      = ['Molt_exit_' +str(x+1) for x in range(molt_exit.shape[0])]
    str_molt_duration  = ['Molt_duration_' +str(x+1) for x in range(molt_duration.shape[0])]

    # Creating pandas dataframe
    molts_df = pd.DataFrame(data = [np.append(molt_entry, molt_exit)],
                      columns = np.append(str_molt_entry, str_molt_exit))

    # Ordering the dataframe
    molts_df = molts_df.sort_values(by = 0, ascending=True, axis=1)

    # Adding the duration
    molts_df = molts_df.join(pd.DataFrame([dict(zip(str_molt_duration, molt_duration))]))

    return molts_df, n_molts
    



# %%
#%%
# Order timepoints (if not done before)
result_file = result_file.sort_values(by = ['Position','Time'])

s_sel = 25 # Check 41 and 44
log_volume = np.log2(result_file[result_file['Position']==s_sel]['volume']).values
time = result_file[result_file['Position']==s_sel]['Time'].values

plt.figure
plt.plot(time,log_volume)

molts, (outliers, log_volume_avg, diff_log_volume) = molt_detection(log_volume)

plt.figure(figsize=(15,8))
plt.plot(time, log_volume)
# plt.plot(time, log_volume_avg)
plt.plot(time[molts],log_volume[molts],'x')

molts_df, n_molts = molt_analysis(molts, time)

molts_df 

# %%
#%% Iterative
# molts_results = pd.DataFrame
result_file['Log_Volume'] = 0

# Order timepoints (if not done before)
result_file = result_file.sort_values(by = ['Position','Time'])

for k_sel, s_sel in enumerate(pd.unique(result_file['Position'])):
    # Compute the volume and times in numpy form
    log_volume = np.log2(result_file[result_file['Position']==s_sel]['volume']).values
    time = result_file[result_file['Position']==s_sel]['Time'].values

    # Detect molts
    molts, (outliers, log_volume_avg, diff_log_volume) = molt_detection(log_volume, molt_th=.2, noise_filt=5, molt_filt= 3)

    # Compute molt times
    molts_df, n_molts = molt_analysis(molts, time, time_th = 10, T_correct = 3)

    # Adding the position information
    molts_df = pd.DataFrame({'Position': [s_sel]}).join(molts_df)

    # Iterating
    if k_sel == 0:
        molts_results = molts_df
    else:
        molts_results = molts_results.append(molts_df)

    # Saving cleaned volume
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Log_Volume'] = log_volume_avg


# Final adjustments
molts_results = molts_results.reset_index()
molts_results = molts_results.drop('index', axis = 1)
molts_results

# %%

molts_results.mean()
# %%

molts_results.std()

# %% Adjusting curves
average_events = molts_results[molts_results.columns[1:7]].mean()



# %%
plt.plot(average_events,molts_results.loc[0, molts_results.columns[1:7] ],'x')
# %%
k_sel = 3
# Creating the variables
X = average_events
X = sm.add_constant(X)
Y = molts_results.loc[k_sel, molts_results.columns[1:7] ]

# Fitting the graph
est = sm.OLS(Y, X)
est2 = est.fit()

print(est2.summary())

plt.figure()
plt.plot(average_events, est2.fittedvalues,'o')
plt.plot(average_events, Y,'x')
# %%
# Compute average worm dynamics
column_names = molts_results.columns[1:7]
average_events = molts_results[column_names].mean()

# Storing values
result_file['Adj_Time'] = 0
molts_results['Growth_rate'] = 0
molts_results['Adj_Hatching'] = 0

# Creating matrix
Y = average_events

# Correcting for timepoints
for k_sel, s_sel in enumerate(np.unique(result_file['Position'])):    
    X = molts_results[molts_results['Position'] == s_sel][column_names].T
    X = sm.add_constant(X)  

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    # Adjusting the times
    times_vals = result_file[result_file['Position']==s_sel]['Time'].values
    times_vals = times_vals*est2.params.iloc[1]+est2.params.iloc[0]

    # Saving
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Adj_Time'] = times_vals

    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Adj_Hatching'] = est2.params.iloc[0]
    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Growth_rate'] = est2.params.iloc[1]


# %% No constant
# Compute average worm dynamics
column_names = molts_results.columns[1:7]
average_events = molts_results[column_names].mean()

# Storing values
result_file['Adj_Time'] = 0
molts_results['Growth_rate'] = 0
molts_results['Adj_Hatching'] = 0

# Creating matrix
Y = average_events

# Correcting for timepoints
for k_sel, s_sel in enumerate(np.unique(result_file['Position'])):    
    X = molts_results[molts_results['Position'] == s_sel][column_names].T

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    # Adjusting the times
    times_vals = result_file[result_file['Position']==s_sel]['Time'].values
    times_vals = times_vals*est2.params.iloc[0]

    # Saving
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Adj_Time'] = times_vals

    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Adj_Hatching'] = 0
    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Growth_rate'] = est2.params.iloc[0]



#%%
sns.histplot(data = molts_results, x = 'Growth_rate', bins = 20)

# %%
# sns.lineplot(data = result_file, x = 'Adj_Time', y = 'volume', hue = 'Position', legend=True)
sns.lineplot(data = result_file[result_file['Position'] <7], x = 'Adj_Time', y = 'Log_Volume', hue = 'Position', legend=True)

#%%
sns.lineplot(data = result_file[result_file['Position'] <7], x = 'Time', y = 'Log_Volume', hue = 'Position', legend=True)

#%%
plt.plot(result_file['Adj_Time'],result_file['volume'],'x')
plt.figure()
plt.plot(result_file['Time'],result_file['volume'],'x')
#%%
n_tps = 173

# Ranges for time
tp_range = .2

# Creating parameters dataframe
param_df = pd.DataFrame(
    index   = np.arange(0, n_tps),
    columns = ['Time',  'Internal', 'External_Factor', 'Signal', 
                        'Internal SE', 'External_Factor SE', 'Signal SE',
                        'Internal pval', 'External_Factor pval', 'Signal pval'])
param_df['Time'] = result_file['Time']

# Creating the extracted values
result_file[['signal_fitted','signal_bg_removed','bacterial_contribution','noise_res']] = 0

for k_tp, tp_sel in enumerate(param_df['Time'][0:n_tps]):
    # Selection of the time
    # result_sel = result_file[result_file['Time']==tp_sel]
    # result_sel = result_file[result_file['Adj_Time'].isin(np.arange(tp_sel,tp_sel+.75,0.25))]
    result_sel = result_file[(result_file['Adj_Time']<(tp_sel+tp_range))&(result_file['Adj_Time']>(tp_sel-tp_range))]

    # Creating the variables
    X = result_sel[['mean_surr','Type']]
    X = sm.add_constant(X)
    Y = result_sel['mean_curved']

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    # Saving the parameters
    param_df.loc[param_df[param_df['Time']==tp_sel].index,param_df.columns[1:]] = \
        np.concatenate([est2.params,est2.bse,est2.pvalues])

    # Extracting the features
    result_file.loc[result_sel.index,'signal_fitted'] = est2.fittedvalues
    result_file.loc[result_sel.index,'signal_bg_removed'] = Y - result_sel['mean_surr']*est2.params[1] - est2.params[0]
    result_file.loc[result_sel.index,'bacterial_contribution'] = result_sel['mean_surr']*est2.params[1]
    result_file.loc[result_sel.index,'noise_res'] = Y - est2.fittedvalues

#%%
sns.lineplot(data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Type', legend=False)
plt.plot(param_df['Time'],param_df['Signal']*(param_df['Signal pval']<.05),'k:')
plt.legend(['Control','hbl-1 mean','computed signal'])
plt.grid()

#%%
sns.lineplot(data = result_file, x = 'Time', y = 'noise_res', hue = 'Type')
plt.grid()


#%%
# df_m = result_file[result_file['Type']==1][['Position', 'Time', 'signal_bg_removed']]
df_m = result_file[['Position', 'Time', 'signal_bg_removed']]
# df_m = df_m.groupby(['Position', 'Time']).mean()
df_m = df_m.groupby(['Time', 'Position']).mean()
df_m = df_m.unstack(level=0)
df_m.columns = df_m.columns.droplevel()
indices = (-df_m.mean(axis = 1)).argsort()
df_m = df_m.iloc[indices]

f, axs = plt.subplots(2,1,figsize=(15,15))
sns.lineplot(ax = axs[0], data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Type', legend=False)
sns.heatmap(ax = axs[1], data=df_m, robust=True, cbar = False)
# axs.margins(x=0)

# #%%
# %%
