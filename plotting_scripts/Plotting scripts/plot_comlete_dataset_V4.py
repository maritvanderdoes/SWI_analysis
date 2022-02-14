#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas.core.frame import DataFrame
import seaborn as sns

# %% Import CSV files
def select_timepoints(frame,hatch,escape):
    if (frame>=hatch)and(frame<escape):
        return (frame-hatch)*min_image_taken/60
    else:
        return None

dir='/Users/Marit/Documents/work/IWM_presentation/'

goodworms=['goodworms_HBL1.csv','goodworms_LIN14.csv','goodworms_lin41.csv','goodworms_LIN42.csv']
results=['Results_parallel_HBL1.csv', 'Results_parallel_LIN14.csv','Results_parallel_lin41.csv','Results_parallel_LIN42.csv']
min_image_taken=15

CSV_quantification_merged=pd.DataFrame()
nr_samples=np.zeros(np.size(goodworms))

for i in range(np.size(goodworms)):
    CSV_quantification = pd.read_csv(dir+results[i])
    CSV_annotation= pd.read_csv(dir+goodworms[i])

    CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
    CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1

    #Add time
    CSV_quantification['Time'] = CSV_quantification.apply(lambda row : select_timepoints(row['Frame'],row['Hatch'],row['Escape']), axis = 1) #add time column
    CSV_quantification= CSV_quantification.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape
    #CSV_quantification= CSV_quantification[CSV_quantification['Time']<60] #selects only timepoints smaller then 40hours

    
    #get mean of controls to substract
    control_only=CSV_quantification[CSV_quantification.Condition=='control']
    to_substract=control_only.groupby('Time').mean()
    to_substract= pd.DataFrame( to_substract["mean_curved"])
    to_substract['Time'] = to_substract.index
    to_substract.index.names = ['Index']
    to_substract2 = to_substract.rename(columns={'mean_curved':'background_signal'},errors="raise")

    #add column to substract, and substract
    CSV_quantification=CSV_quantification.merge(to_substract2,on=['Time'], how='inner') #merges two datasets
    CSV_quantification["final intensity (a.u.)"] = CSV_quantification['mean_curved']-CSV_quantification['background_signal']

    #delete controls lines
    sample_only=CSV_quantification[CSV_quantification.Condition!='control']
    CSV_quantification_merged = pd.concat([CSV_quantification_merged, sample_only])



# filtering dataset for quality and Time

#CSV_quantification_merged= CSV_quantification_merged.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape
#CSV_quantification_merged= CSV_quantification_merged[CSV_quantification_merged['Time']<40] #selects only worms in between hatch and excape
#delete rows that are not good
#CSV_quantification_merged= CSV_quantification_merged[CSV_quantification !=0]
#CSV_quantification_merged= CSV_quantification_merged.dropna(axis='index') 


# %% plots whole dataset
plt.figure(1)
sns.set_theme(style="whitegrid")
fig, axs = plt.subplots(1, 2, figsize=(20,8),sharex=False)
sns.lineplot(ax=axs[0],x='Time',y="final intensity (a.u.)",hue='Condition',data=CSV_quantification_merged,legend=True)
axs[0].axvline(12.5,color='gray')
axs[0].axvline(21,color='gray')
axs[0].axvline(30,color='gray')
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',hue='Condition',data=CSV_quantification_merged,legend=False)
axs[1].axvline(12.5,color='gray')
axs[1].axvline(21,color='gray')
axs[1].axvline(30,color='gray')
axs[1].semilogy()
# %% select the conditions you want to plot


#CSV_selected_condition=CSV_quantification_merged[(CSV_quantification_merged.Condition=='LIN_14')
#(CSV_quantification_merged.Condition=='LIN_28') ]

CSV_selected_condition=CSV_quantification_merged[(CSV_quantification_merged.Condition=='LIN_28')]

plt.figure(2)
#sns.set_theme(style="whitegrid")
fig, axs = plt.subplots(2, 1, figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y="final intensity (a.u.)",hue='Condition',data=CSV_selected_condition,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',hue='Condition',data=CSV_selected_condition,legend=False).set_title('volume')
axs[1].semilogy()

axs[0].axvline(12.5,color='gray')
axs[0].axvline(20,color='gray')
axs[0].axvline(28,color='gray')
axs[0].axvline(40,color='gray')
axs[1].axvline(12.5,color='gray')
axs[1].axvline(20,color='gray')
axs[1].axvline(28,color='gray')
axs[1].axvline(40,color='gray')


# %% Plot individual worms or conditions
CSV_selected_condition=CSV_quantification_merged[(CSV_quantification_merged.Condition=='HBL_1')
                                    | (CSV_quantification_merged.Condition=='LIN_29')
                                    | (CSV_quantification_merged.Condition=='LIN_42')
                                                                                ]


plt.figure(3)
sns.set_theme(style="whitegrid")
fig, axs = plt.subplots(2, 1, figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y="final intensity (a.u.)",hue='Condition',data=CSV_selected_condition,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',hue='Condition',data=CSV_selected_condition,legend=False).set_title('volume')
axs[1].semilogy()

# %% curved 3D
plt.figure(1)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1, figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='mean_curved',hue='Condition',data=CSV_quantification2,legend=True).set_title('properties curved 3D')
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',hue='Condition',data=CSV_quantification2,legend=False)
axs[1].semilogy()

# %% volume comparision
plt.figure(2)
sns.lineplot(x='Time',y='volume_curved',data=CSV_quantification2,legend=True).set_title('volume comparison 2D vs 3D')
sns.lineplot(x='Time',y='volume_curved_2D',data=CSV_quantification2,legend=True)
sns.lineplot(x='Time',y='area_curved_2D',data=CSV_quantification2,legend=True)
plt.legend(['volume_curved', 'volume_curved_2D','area_curved_2D'])
plt.semilogy()


# %% curved with head vs curved without head
plt.figure(3)
sns.lineplot(x='Time',y='mean_curved',data=CSV_quantification2[CSV_quantification2['Condition']=="HBL_1"],legend=True).set_title('intenisty comparision minus head/tail')
sns.lineplot(x='Time',y='mean_curved_ht',data=CSV_quantification2[CSV_quantification2['Condition']=="HBL_1"],legend=True)
sns.lineplot(x='Time',y='mean_straightened',data=CSV_quantification2[CSV_quantification2['Condition']=="HBL_1"],legend=True)
sns.lineplot(x='Time',y='mean_straightened_ht',data=CSV_quantification2[CSV_quantification2['Condition']=="HBL_1"],legend=True)
plt.legend(['mean_curved', 'mean_curved_ht','mean_straightened','mean_straightened_ht'])


# %%
tp20=CSV_quantification2[CSV_quantification2['Time']==2]
volume_curved=np.mean(tp20['volume_curved'])
print(volume_curved)
area_curved = np.mean(tp20['area_curved_2D'])
hight=volume_curved/area_curved
Length=(np.pi*area_curved**2)/(4*volume_curved)


# %%
np.where(CSV_annotation.Quality[]==1)

# %%
np.size(np.where((CSV_annotation.Quality[(CSV_annotation.Condition!='control')]==1)))
# %%(CSV_annotation.Condition!='control')
