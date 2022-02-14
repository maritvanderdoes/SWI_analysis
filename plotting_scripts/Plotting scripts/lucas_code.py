#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 

@author: Lucas J. Morales
"""


#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as sig
import seaborn as sns
import statsmodels.api as sm

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

    #outliers
    #log_volume
    #diff_log_volume
    return molts

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
    

#%%
# Reading CSV file
# dir='D:\Analysis_HBL1\Results_quant_hist.csv'
# dir='D:\Analysis_HBL1\Results_quant_hist_3_5.csv'
dir='/Volumes/ggrossha/Marit/HetP_Quant/20201224_HBL1_analysis/Marits_analysis/V3/Results_HBL1_V3.csv'

result_file = pd.read_csv(dir)

# Defining the samples
result_file['Type'] = (result_file['Position']<40).astype(int)

# Creating adjusted metrics
result_file['Time'] = result_file['Frame']
result_file['Adj_Time'] = result_file['Frame']
result_file['Log_Volume'] = np.log2(result_file['volume_curved'])

# Cleaning the file
result_file = result_file[['Position', 'Type', 'Frame','Time', 'Adj_Time', 'mean_curved','mean_background', 'volume_curved', 'Log_Volume']]

# Creating new columns
result_file[['Volume_Outlier', 'Signal_fitted','Signal_bg_removed','Signal_bg_removed_med',
    'Bacterial_contribution','Noise_res','Sig2Noise_Ratio']] = np.nan

# Adjusting for timepoints
for k_v, k_sel in enumerate(np.unique(result_file['Position'])):
    times_vals = result_file[result_file['Position']==k_sel]['Frame'].values
    times_vals = (times_vals - times_vals[0])/4
    result_file.loc[result_file[result_file['Position']==k_sel].index,'Time'] = times_vals


# Order timepoints
result_file = result_file.sort_values(by = ['Position','Time'])

#%%
# Molt calculation
for k_sel, s_sel in enumerate(pd.unique(result_file['Position'])):
    # Compute the volume and times in numpy form
    log_volume = np.log2(result_file[result_file['Position']==s_sel]['volume_curved']).values
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

    # Saving cleaned metrics
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Log_Volume'] = log_volume_avg
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Volume_Outlier'] = outliers


# Final adjustments
molts_results = molts_results.reset_index()
molts_results = molts_results.drop('index', axis = 1)
molts_results

# Adding the average worm
average_worm = molts_results.mean()
average_worm['Position'] = 'X'
molts_results = molts_results.append(average_worm, ignore_index=True)

# Creating the average worm dynamics
n_molts = 3
time_vec = pd.unique(result_file['Time'])
average_molts = np.sum([(time_vec >= average_worm[2*k+1])&(time_vec <= average_worm[2*k+2]) for k in range(0,n_molts)],axis = 0)

# %%
# Force Hatch to be zero?
hatch_zero = True

# Compute average worm dynamics
column_names = molts_results.columns[1:7]
average_events = molts_results[column_names].mean()

# Storing values
molts_results['Growth_rate'] = 1
molts_results['Adj_Hatching'] = 0

# Creating matrix
Y = average_events

# Correcting for timepoints
for k_sel, s_sel in enumerate(np.unique(result_file['Position'])):    
    X = molts_results[molts_results['Position'] == s_sel][column_names].T
    if not hatch_zero:
        X = sm.add_constant(X,prepend=False)  

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    if hatch_zero:
        coordinate_at_origin = 0
    else:
        coordinate_at_origin = est2.params.iloc[1]

    # Adjusting the times
    times_vals = result_file[result_file['Position']==s_sel]['Time'].values
    times_vals = times_vals*est2.params.iloc[0]+coordinate_at_origin

    # Saving
    result_file.loc[result_file[result_file['Position']==s_sel].index,'Adj_Time'] = times_vals

    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Adj_Hatching'] = coordinate_at_origin
    molts_results.loc[molts_results[molts_results['Position']==s_sel].index,'Growth_rate'] = est2.params.iloc[0]


# Rounding Adj_time
result_file['Adj_Time_Round'] = np.round(result_file['Adj_Time']*4)/4

#%%
sns.histplot(data = molts_results, x = 'Growth_rate', bins = 20)

#%%
sns.lineplot(data = result_file, x = 'Adj_Time', y = 'Log_Volume', hue = 'Position', legend=True)
plt.figure()
sns.lineplot(data = result_file, x = 'Time', y = 'Log_Volume', hue = 'Position', legend=True)


# %%
n_tps = 40*4

# Time grouping
time_var = 'Adj_Time'

# Time range (Average Time ± tp_range)
tp_range = .25

# Creating parameters dataframe
param_df = pd.DataFrame(
    index   = np.arange(0, n_tps),
    columns = ['Time',  'Internal', 'External_Factor', 'Signal', 
                        'Internal SE', 'External_Factor SE', 'Signal SE',
                        'Internal pval', 'External_Factor pval', 'Signal pval'])
param_df['Time'] = result_file['Time']

for k_tp, tp_sel in enumerate(param_df['Time'][0:n_tps]):
    # Selection of the reference time
    result_sel = result_file[(result_file[time_var]<=(tp_sel+tp_range))&(result_file['Adj_Time']>=(tp_sel-tp_range))]

    # Remove volume outliers
    # result_sel = result_sel[result_sel['Volume_Outlier']==False]

    # Selection of the target times
    target_idx = result_sel[result_sel['Adj_Time_Round'] == tp_sel].index
    positions = [np.where(result_sel.index == k_sel)[0][0] for k_sel in target_idx]

    # Creating the variables
    X = result_sel[['mean_background','Type']]
    X = sm.add_constant(X)
    Y = result_sel['mean_curved']

    # Fitting the graph
    est = sm.OLS(Y, X)
    est2 = est.fit()

    # Saving the parameters
    param_df.loc[param_df[param_df['Time']==tp_sel].index,param_df.columns[1:]] = \
        np.concatenate([est2.params,est2.bse,est2.pvalues])

    # Extracting teh points
    idx_sel = result_sel.index[positions]

    # Computing the parameters
    external_bg = result_sel.loc[idx_sel,'mean_surr']*est2.params[1]
    internal_bg = est2.params[0]
    background = external_bg + internal_bg
    fitted_fluorescence = est2.fittedvalues[idx_sel]
    signal = fitted_fluorescence - background
    noise = result_sel.loc[idx_sel,'mean_curved'] - fitted_fluorescence
    sig2noise = signal/noise

    # Saving the features
    result_file.loc[idx_sel,'Signal_fitted'] = fitted_fluorescence
    result_file.loc[idx_sel,'Sig2Noise_ratio'] = sig2noise
    result_file.loc[idx_sel,'Bacterial_contribution'] = external_bg
    result_file.loc[idx_sel,'Noise_res'] = noise
    # result_file.loc[idx_sel,'Signal_bg_removed'] = fitted_fluorescence + (signal/fitted_fluorescence)*noise
    result_file.loc[idx_sel,'Signal_bg_removed'] = result_sel.loc[idx_sel,'mean_curved'] - background

#%% Removing outliers from signal_bg_removed
for s_sel in pd.unique(result_file['Position']):
    vals2filt = result_file[result_file['Position']==s_sel]['Signal_bg_removed']
    (idx_vals, vals2filt) = (vals2filt.index,vals2filt.values)
    med_comp = sig.medfilt(vals2filt,5)
    act_comp = vals2filt-med_comp
    pk_sel = act_comp>2*np.nanstd(act_comp)
    vals2filt[pk_sel] = med_comp[pk_sel]
    result_file.loc[idx_vals,'Signal_bg_removed_med'] = vals2filt


#%%
# sns.lineplot(data = result_file, x = 'Time', y = 'signal_bg_removed', hue = 'Type', legend=True)
# # plt.plot(param_df['Time'],param_df['Signal']*(param_df['Signal pval']<.05),'k:')
# # plt.legend(['Control','hbl-1 mean','computed signal'])
# plt.grid()
# plt.gca().set_ylim(bottom=-1.5)

# %%
sns.lineplot(data = result_file, x = 'Adj_Time_Round', y = 'Signal_bg_removed', hue = 'Type', legend=False)
plt.fill_between(time_vec, 1.1*average_molts*param_df['Signal'].max(), -1.5*average_molts, facecolor="red", alpha = 0.25)
# plt.plot(param_df['Time'],param_df['Signal']*(param_df['Signal pval']<.05),'k:')
# plt.legend(['Control','hbl-1 mean','computed signal'])
plt.grid()
plt.gca().set_ylim(bottom=-1.5)

# %%
sns.lineplot(data = result_file, x = 'Adj_Time_Round', y = 'Log_Volume', hue = 'Type', legend=False)
plt.fill_between(time_vec, 21*average_molts, 15*average_molts, facecolor="red", alpha = 0.25)
plt.grid()
plt.gca().set_ylim(bottom=15)

#%%
# df_m = result_file[result_file['Type']==1][['Position', 'Time', 'signal_bg_removed']]
df_m = result_file[['Position', 'Adj_Time_Round', 'Signal_bg_removed_med']]
# df_m = df_m.groupby(['Position', 'Time']).mean()
df_m = df_m.groupby(['Adj_Time_Round', 'Position']).mean()
df_m = df_m.unstack(level=0)
df_m.columns = df_m.columns.droplevel()
indices = (-df_m.mean(axis = 1)).argsort()
df_m = df_m.iloc[indices]

result_plot = result_file.copy()
result_plot['Time_Plot'] = 4*result_plot['Adj_Time_Round']

f, axs = plt.subplots(2,1,sharex = True, figsize=(15,15))
sns.lineplot(ax = axs[0], data = result_plot, x = 'Time_Plot', y = 'Signal_bg_removed_med', hue = 'Type', legend=False)
axs[0].fill_between(time_vec*4, 1.1*average_molts*param_df['Signal'].max(), -1.5*average_molts, facecolor="red", alpha = 0.25)
axs[0].grid()
sns.heatmap(ax = axs[1], data=df_m, robust=True, cbar = False)

#%%
sns.lineplot(data = result_file[result_file['Type']==1], x = 'Adj_Time_Round', y = 'Signal_bg_removed', hue = 'Position', legend=False)
plt.grid()

# %%
denoise = 5
krn_siz = 21
trend_signal = np.convolve(param_df['Signal'],np.ones(krn_siz)/krn_siz,'same')
smooth_signal = np.convolve(param_df['Signal'],np.ones(denoise)/denoise,'same')
plt.plot(param_df['Time'],smooth_signal-trend_signal)
plt.fill_between(time_vec, average_molts, -average_molts, facecolor="red", alpha = 0.25)
# plt.plot(time_vec, average_molts)


plt.figure()
plt.plot(param_df['Time'],smooth_signal)
plt.plot(param_df['Time'],trend_signal)
# %%
sns.scatterplot(data=result_file[result_file['Adj_Time']==32],x = 'mean_background', y = 'mean_curved', hue = 'Type')
# %%
nan_vals = result_file[np.isnan(result_file['Signal_bg_removed'])]
# %%
avg_vol = result_file[result_file['Adj_Time_Round'].isin(pd.unique(param_df['Time']))]
avg_vol = avg_vol[['Adj_Time_Round', 'Log_Volume']]
avg_vol = avg_vol.groupby(['Adj_Time_Round']).mean()


avg_vol =avg_vol['Log_Volume'].values
plt.plot(param_df['Time'],param_df['Signal']+avg_vol)

plt.figure()
plt.plot(param_df['Time'],avg_vol)
# %%
