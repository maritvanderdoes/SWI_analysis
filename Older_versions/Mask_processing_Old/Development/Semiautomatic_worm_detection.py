#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from utils import molt_detection, molt_analysis
from utils import multiple_regression

def escape_determination(df, threshold = 8, tps = 'Frame', volume_diff = 'volume_diff', area = 'mask_area', add_end = True):
    exit_likelihood = (-(df[volume_diff]/df[area])>threshold) + \
                      (  df[volume_diff].isnull().astype(float).diff()>0)
    
    escape_estimates = df.loc[exit_likelihood,tps].values

    if add_end or (escape_estimates.size==0):
        escape_estimates = np.append(np.array([df[tps].max()-1]),escape_estimates)

    return escape_estimates

def hatching_determination_motion(df, threshold = 8, window = 3, eccentricity_th = .9, tps = 'Frame', displacement = 'displacement', eccentricity = 'mask_eccentricity', add_start = True):
    # On motion
    df["egg_likelihood_motion"] = ((df[displacement]<threshold).astype(float).rolling(window,center = True).mean().values>.5).astype(float)
    hatching_estimates_motion = df["egg_likelihood_motion"].diff()<0

    hatching_estimates = df.loc[hatching_estimates_motion,tps].values

    if add_start or (df.loc[df[tps]==1,eccentricity].values > eccentricity_th):
        hatching_estimates = np.append([2],hatching_estimates)

    return hatching_estimates

def hatching_determination_eccentricity(df, threshold = 2, window = 1, eccentricity_th = .9, tps = 'Frame', eccentricity = 'mask_eccentricity', add_start = True):

    fun = np.log(df[eccentricity]).rolling(window, center = True).mean().diff()
    hatching_estimates = df.loc[fun<-threshold*fun.std(),tps].values

    if add_start or (df.loc[df[tps]==1,eccentricity].values > eccentricity_th):
        hatching_estimates = np.append([2],hatching_estimates)

    return hatching_estimates

# %%
dir='D:/GitHub/SWI_analysis/Analysis V3/'
filename = 'cleanup_Results_LJ_Lin28_Project_210305.csv'

# Opening the file
preanalysis_csv = pd.read_csv(dir+filename)
preanalysis_csv = preanalysis_csv.sort_values(['Position','Frame'])

# Checking NAs
print(preanalysis_csv.isnull().sum())

# %%
position_sel = 71 # 50, 68-74
motion_level = 7
motion_window = 3
ec_th = 0.9
vol_change_th = 8

preanalysis_csv.groupby('Position').median()['mask_eccentricity'].sort_values()

selection_csv = preanalysis_csv[preanalysis_csv['Position']==position_sel]
print(selection_csv.isnull().sum())

# New metrics
selection_csv['displacement'] = \
    np.sqrt((selection_csv['mask_centroid_x'].diff())**2+\
            (selection_csv['mask_centroid_y'].diff())**2)
selection_csv['volume_diff'] = selection_csv['mask_area'].diff()

# Estimations
selection_csv['egg_likelihood'] = selection_csv['displacement']<motion_level
selection_csv['egg_likelihood'] = (selection_csv['egg_likelihood'].astype(float).rolling(motion_window,center = True).mean().values>.5).astype(float)
selection_csv['exit_likelihood'] = (-selection_csv['volume_diff'] /selection_csv['mask_area']>vol_change_th ) + (selection_csv['volume_diff'].isnull().astype(float).diff()>0)

fig, axs = plt.subplots(5,1,figsize = (6,12),sharex = True)
axs[0].semilogy(selection_csv['Frame'],selection_csv['volume_curved'])
axs[0].semilogy(selection_csv['Frame'],selection_csv['mask_area'])
# axs[0].semilogy(selection_csv['Frame'],selection_csv['mask_total'])
axs[1].plot(selection_csv['Frame'],selection_csv['volume_diff'] )
axs[2].plot(selection_csv['Frame'],selection_csv['mask_mean'] )

axs[3].plot(selection_csv['Frame'],selection_csv['mask_eccentricity'])
axs[3].plot(selection_csv['Frame'],selection_csv['egg_likelihood'],'r')


axs[4].plot(selection_csv['Frame'],selection_csv['displacement'])
axs[4].plot(selection_csv['Frame'],selection_csv['egg_likelihood']*selection_csv['displacement'].max(),'r')
axs[4].plot(selection_csv['Frame'],selection_csv['exit_likelihood']*selection_csv['displacement'].max(),'c')


# %% Selecting egg
# goodworms_pre = pd.DataFrame(columns = ['Position','Hat'])
hatching_estimates = selection_csv.loc[selection_csv['egg_likelihood'].diff()<0,'Frame'].values
selection_csv.loc[selection_csv['egg_likelihood'].diff()<0,'mask_eccentricity'].values
if selection_csv.loc[selection_csv['Frame']==1,'mask_eccentricity'].values >ec_th:
    hatching_estimates = np.append([2],hatching_estimates)
plt.plot(selection_csv['Frame'],selection_csv['egg_likelihood'].diff())

new_hatches_m = hatching_determination_motion(selection_csv, motion_level)
new_hatches_e = hatching_determination_eccentricity(selection_csv, threshold = 2)

hatches_merged = np.unique(np.concatenate([new_hatches_m,new_hatches_e]))

# %%
escape_estimates = selection_csv.loc[selection_csv['exit_likelihood'],'Frame'].values
if escape_estimates.shape[0]==0:
    escape_estimates = np.array([2])

new_escapes = escape_determination(selection_csv, vol_change_th)
# worm_entry = selection_csv.iloc[np.arange(0,-2),'volume_diff'] /selection_csv.iloc[np.arange(1,-1),'mask_area']>100
plt.figure()
plt.plot(selection_csv['Frame'],selection_csv['exit_likelihood'])
# plt.plot(selection_csv['Frame'],-selection_csv['volume_diff'] /selection_csv['mask_area'])
# plt.plot(selection_csv['Frame'],-selection_csv['volume_diff'] /selection_csv['mask_area'])
# plt.ylim([0,20])
# plt.plot(selection_csv.iloc[np.arange(0,-2),'Frame'],worm_entry)
# %%
plt.semilogy(selection_csv['Frame'],selection_csv['volume_curved'])
plt.axvline(hatching_estimates[0],color = 'k')
if hatching_estimates.shape[0]>1:
    plt.axvline(hatching_estimates[1],color = 'gray')
if (selection_csv.loc[selection_csv['Frame']==1,'mask_eccentricity'].values >ec_th)&(hatching_estimates.shape[0]>1):
    plt.axvline(hatching_estimates[2],color = 'gray')
plt.axvline(escape_estimates[0],color = 'r')
if escape_estimates.shape[0]>1:
    plt.axvline(escape_estimates[1],color = 'orange')

plt.xlabel('Frame')
plt.ylabel('Volume')
# plt.savefig('D:/GitHub/SWI_analysis/output/'+basename+'_Feature_analysis_s'+str(position_sel)+'.pdf')

#%% 
plt.plot(-selection_csv['volume_diff'] /selection_csv['mask_area'])

#%%%
plt.plot(selection_csv['Frame'],selection_csv['mask_eccentricity'].rolling(4).mean().diff()**2)

#%% Recomputing the metrics
ecc_window = 5
hatching_prob_displacement = selection_csv['displacement']<motion_level
hatching_prob_eccentricity = selection_csv['mask_eccentricity'].rolling(ecc_window, center = True).mean().diff()

print(selection_csv.loc[hatching_prob_eccentricity<-hatching_prob_eccentricity.std(),'Frame'])
plt.plot(selection_csv['Frame'],hatching_prob_eccentricity)
# plt.plot(selection_csv['Frame'],selection_csv['mask_eccentricity'])

#%% Refining metrics
mask_window = 3
selection_csv['mask_area_log'] = np.log(selection_csv['mask_area'])
selection_csv['mask_area_w'] = selection_csv['mask_area'].rolling(mask_window,center = True).mean()
selection_csv['mask_area_w_log'] = np.log(selection_csv['mask_area_w'])
selection_csv['mask_area_log_w'] = selection_csv['mask_area_log'].rolling(mask_window,center = True).mean()


fig, axs = plt.subplots(4,1,figsize = (6,12),sharex = True)
axs[0].plot(selection_csv['Frame'],selection_csv['mask_area'])
axs[0].plot(selection_csv['Frame'],selection_csv['mask_area_w'])

axs[1].plot(selection_csv['Frame'],selection_csv['mask_area_log'])

axs[2].plot(selection_csv['Frame'],selection_csv['mask_area_w_log'])
axs[2].plot(selection_csv['Frame'],selection_csv['mask_area_log_w'])

axs[3].plot(selection_csv['Frame'],selection_csv['mask_area_log'].diff())
axs[3].plot(selection_csv['Frame'],selection_csv['mask_area_w_log'].diff())
axs[3].plot(selection_csv['Frame'],selection_csv['mask_area_log_w'].diff())

std_th = 1
param = 'mask_area_log' 
param = 'mask_eccentricity' 
fun = selection_csv[param].diff()
fun = selection_csv[param].rolling(mask_window,center = True).mean().diff()
plt.figure()
# plt.plot(selection_csv['Frame'],fun<-std_th*fun.std())
plt.plot(selection_csv['Frame'],fun>std_th*fun.std())

print(hatching_estimates)
print(escape_estimates)
# print(selection_csv.loc[fun<-fun.std(),'Frame'])
print(selection_csv.loc[fun>std_th*fun.std(),'Frame'])

#%%
std_th = 1
mask_window = 3

param = 'displacement' 
fun = selection_csv[param].rolling(mask_window,center = True).mean().diff()
hatching_prob = fun>std_th*fun.std()
hatching_events = selection_csv.loc[hatching_prob,'Frame']
print(hatching_events)

#%%
std_th = 1

plt.figure()
param = 'mask_area_log' 
fun = selection_csv[param].diff()
plt.plot(selection_csv['Frame'],fun<-std_th*fun.std())

print(hatching_estimates)
print(escape_estimates)
print(selection_csv.loc[fun<-fun.std(),'Frame'])

#%%
std_th = 1
std_window = 3

param = 'mask_eccentricity' 
fun = selection_csv[param].diff()
fun_std = selection_csv[param].rolling(std_window, center = True).std()
plt.figure()
plt.plot(selection_csv['Frame'],fun<-std_th*fun.std())


plt.figure()
plt.plot(selection_csv['Frame'],fun)
plt.figure()
plt.plot(selection_csv['Frame'],fun_std)

plt.figure()
plt.plot(selection_csv['Frame'],np.sign(fun)*fun>std_th*fun_std)
plt.plot(selection_csv['Frame'],np.sign(fun)*fun>std_th*fun.std(),'x:')

print(hatching_estimates)
print(escape_estimates)
print(selection_csv.loc[fun<-fun.std(),'Frame'])

#%%
fun_anal = selection_csv['displacement'].rolling(3,center = True).mean()
plt.semilogy(selection_csv['Frame'],fun_anal)

#%%
diff_anal = np.log(fun_anal).diff()
diff_std = diff_anal.std()
plt.plot(selection_csv['Frame'],diff_anal>3*diff_std)
selection_csv.loc[(diff_anal>3*diff_std),'Frame']

# %%

locidx = selection_csv['Frame'].isin(np.arange(0,hatching_estimates[0]))
plt.figure()
plt.plot(selection_csv.loc[locidx,'mask_centroid_x'],selection_csv.loc[locidx,'mask_centroid_y'],'rx-')
plt.xlim([0, 1200])
plt.ylim([0, 1200])
# %%

# %%
from utils import image_lists, single_image_lists, read_image_and_metadata
import os
import io
import imageio
import scipy.misc as misc

# load parameters
dirpath         = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210305_Compressed'
channel_GFP     = 'w1Lucas-sim-488-561'
channel_mcherry = 'w2Lucas-sim-561-488'

# # Additional parameters
# n_workers = 12
# data_format = 'ts'

# s_sel = position_sel
# t_sel = hatching_estimates[0]

# # list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
# (list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel2 = channel_GFP, s_sel = s_sel, t_sel = t_sel, channel3 = None, data_format = data_format)

#%%
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)

#%% Iterating files
colset = ['Frame', 'Position', 'Condition', 'Name']
list_df = pd.DataFrame(columns = colset)

for k_f in np.arange(0,len(list_mcherry)):
    if np.mod(k_f,200)==0:
        print(k_f)
    (directory, filename) = list_mcherry[k_f].split('\\')
    (basename, rest) = filename.split('_t')
    (frame_set, rest) = rest.split('_s')
    (position_set, rest) = rest.split('_l')
    (condition_set, rest) = rest.split('p')

    list_df=list_df.append(dict(zip(colset,[int(frame_set), int(position_set), condition_set, list_mcherry[k_f]])), ignore_index=True)

#%%

n_plots_th = 3
n_plots = max(min(hatching_estimates.shape[0],n_plots_th),
              min(escape_estimates.shape[0],n_plots_th))

# Plotting hatching estimates
fig, axs = plt.subplots(n_plots,6,figsize = (12,n_plots*3),sharex = True, sharey= True, squeeze=False)
subset = list_df[list_df['Position']==position_sel]

for prob_set in range(0,min(hatching_estimates.shape[0],n_plots_th)):
    t_sel = hatching_estimates[prob_set]

    print('Plotting hatching: '+str(prob_set+1))
    for k_t, t_set in enumerate([t_sel-1,t_sel,t_sel+1]):
        filename_it = subset.loc[subset['Frame']==t_set,'Name'].values[0]
        img_mcherry = np.array(imageio.mimread(filename_it))

        axs[prob_set,k_t].imshow(img_mcherry.max(axis = 0))
        axs[prob_set,k_t].set_title('f: '+str(t_set))
        axs[prob_set,k_t].axis('off')

for prob_set in range(0,min(escape_estimates.shape[0],n_plots_th)):
    t_sel = escape_estimates[prob_set]

    print('Plotting escaping: '+str(prob_set+1))
    
    for k_t, t_set in enumerate([t_sel-1,t_sel,t_sel+1]):
        filename_it = subset.loc[subset['Frame']==t_set,'Name'].values[0]
        img_mcherry = np.array(imageio.mimread(filename_it))

        axs[prob_set,k_t+3].imshow(img_mcherry.max(axis = 0))
        axs[prob_set,k_t+3].set_title('f: '+str(t_set))
        axs[prob_set,k_t+3].axis('off')

fig.suptitle('Sample: '+str(position_sel)+', Condition: '+subset.loc[subset['Frame']==t_set,'Condition'].values[0]+\
    ', NAs: '+str(selection_csv.isnull().sum()['mask_mean'])+ ', mean_ecc: '+str(np.round(selection_csv['mask_eccentricity'].mean(),2)))

line = plt.Line2D((.51,.51),(.1,.9), color="k", linewidth=3)
fig.add_artist(line)


fig.savefig('D:/GitHub/SWI_analysis/output/'+basename+'_Sample_analysis_s'+str(position_sel)+'.pdf')

# %%
# FUll code

dirsel = 'D:/GitHub/SWI_analysis/output/'
dirsel = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_HetPathQuant_Results'

motion_level = 8
motion_window = 5
ec_th = 0.9
vol_change_th = 10
t_th = 255

# iterating_vec = np.arange(4,preanalysis_csv['Position'].max()+1)

iterating_vec = preanalysis_csv.groupby('Position').median()['mask_eccentricity'].sort_values(ascending=False).index

for position_sel in [42]:#iterating_vec[0:]:
    print(position_sel)
    selection_csv = preanalysis_csv[preanalysis_csv['Position']==position_sel]

    # New Metrics
    selection_csv['displacement'] = \
        np.sqrt((selection_csv['mask_centroid_x'].diff())**2+\
                (selection_csv['mask_centroid_y'].diff())**2)
    selection_csv['volume_diff'] = selection_csv['mask_area'].diff()

    # Estimations
    selection_csv['egg_likelihood'] = selection_csv['displacement']<motion_level
    selection_csv['egg_likelihood'] = (selection_csv['egg_likelihood'].astype(float).rolling(motion_window,center = True).mean().values>.5).astype(float)
    selection_csv['exit_likelihood'] = (-selection_csv['volume_diff'] /selection_csv['mask_area']>vol_change_th) + (selection_csv['volume_diff'].isnull().astype(float).diff()>0)

    # Hatching
    hatching_estimates = selection_csv.loc[selection_csv['egg_likelihood'].diff()<0,'Frame'].values

    if selection_csv.loc[selection_csv['Frame']==1,'mask_eccentricity'].values >ec_th:
        hatching_estimates = np.append([2],hatching_estimates)

    if hatching_estimates.size==0:
        hatching_estimates = np.array([2])

    # Exit
    escape_estimates = selection_csv.loc[selection_csv['exit_likelihood'],'Frame'].values
    if escape_estimates.shape[0]==0:
        escape_estimates = np.array([2])

    if escape_estimates.size==0:
        escape_estimates = np.array([selection_csv['Frame'].max()-1])

    plt.figure(figsize=(10,6))
    plt.semilogy(selection_csv['Frame'],selection_csv['volume_curved'])
    plt.axvline(hatching_estimates[0],color = 'k')
    if hatching_estimates.shape[0]>1:
        plt.axvline(hatching_estimates[1],color = 'gray')
    if (selection_csv.loc[selection_csv['Frame']==1,'mask_eccentricity'].values >ec_th)&(hatching_estimates.shape[0]>1):
        plt.axvline(hatching_estimates[1],color = 'gray')
    plt.axvline(escape_estimates[0],color = 'r')
    if escape_estimates.shape[0]>1:
        plt.axvline(escape_estimates[1],color = 'orange')

    plt.xlabel('Frame')
    plt.ylabel('Volume')
    plt.savefig(dirsel+'/worm_detect/'+basename+'_Feature_analysis_s'+str(position_sel)+'.pdf')
    plt.close()

    # Check images
    n_plots_th = 3
    n_plots = max(min(hatching_estimates.shape[0],n_plots_th),
                min(escape_estimates.shape[0],n_plots_th))

    # Plotting hatching estimates
    fig, axs = plt.subplots(n_plots,6,figsize = (12,n_plots*3),sharex = True, sharey= True, squeeze=False)
    subset = list_df[list_df['Position']==position_sel]

    for prob_set in range(0,min(hatching_estimates.shape[0],n_plots_th)):
        t_sel = hatching_estimates[prob_set]

        print('Plotting hatching: '+str(prob_set+1))
        for k_t, t_set in enumerate([t_sel-1,t_sel,t_sel+1]):
            if (t_set>0)&(t_set<=subset['Frame'].max()):
                filename_it = subset.loc[subset['Frame']==t_set,'Name'].values[0]
                img_mcherry = np.array(imageio.mimread(filename_it))

                axs[prob_set,k_t].imshow(img_mcherry.max(axis = 0))
                axs[prob_set,k_t].set_title('f: '+str(t_set))
                axs[prob_set,k_t].axis('off')
            elif (t_set<subset['Frame'].max()):
                axs[prob_set,k_t+3].imshow(img_mcherry.max(axis = 0))
                axs[prob_set,k_t+3].set_title('f: '+str(t_set))
                axs[prob_set,k_t+3].axis('off')  

    for prob_set in range(0,min(escape_estimates.shape[0],n_plots_th)):
        t_sel = escape_estimates[prob_set]

        print('Plotting escaping: '+str(prob_set+1))
        
        for k_t, t_set in enumerate([t_sel-1,t_sel,t_sel+1]):
            if (t_set>0)&(t_set<=subset['Frame'].max()):
                filename_it = subset.loc[subset['Frame']==t_set,'Name'].values[0]
                img_mcherry = np.array(imageio.mimread(filename_it))

                axs[prob_set,k_t+3].imshow(img_mcherry.max(axis = 0))
                axs[prob_set,k_t+3].set_title('f: '+str(t_set))
                axs[prob_set,k_t+3].axis('off')
            elif (t_set<subset['Frame'].max()):
                axs[prob_set,k_t+3].imshow(img_mcherry.max(axis = 0))
                axs[prob_set,k_t+3].set_title('f: '+str(t_set))
                axs[prob_set,k_t+3].axis('off')  

    fig.suptitle('Sample: '+str(position_sel)+', Condition: '+subset.loc[subset['Frame']==t_sel,'Condition'].values[0]+\
        ', NAs: '+str(selection_csv.isnull().sum()['mask_mean'])+ ', mean_ecc: '+str(np.round(selection_csv['mask_eccentricity'].mean(),2)))

    line = plt.Line2D((.515,.515),(.1,.9), color="k", linewidth=3)
    fig.add_artist(line)


    fig.savefig(dirsel+'/worm_detect/'+basename+'_Sample_analysis_s'+str(position_sel)+'.pdf', dpi = 200)
    plt.close()

# %%
