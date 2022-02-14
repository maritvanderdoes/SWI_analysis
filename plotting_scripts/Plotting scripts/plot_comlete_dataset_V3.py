#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# %% Import CSV files
dir='/Users/Marit/Documents/programming/python/'

CSV_quantification = pd.read_csv(dir+"Results_parallel_HBL1.csv")
CSV_annotation= pd.read_csv(dir+"goodworms_HBL1.csv")
min_image_taken=15

# filtering dataset for quality and Time
def select_timepoints(frame,hatch,escape):
    if (frame>=hatch)and(frame<escape):
        return (frame-hatch)*min_image_taken/60
    else:
        return None

CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1
CSV_quantification['Time'] = CSV_quantification.apply(lambda row : select_timepoints(row['Frame'],row['Hatch'],row['Escape']), axis = 1) #add time column
CSV_quantification= CSV_quantification.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape

#delete rows that are not good
CSV_quantification2= CSV_quantification[CSV_quantification !=0]
CSV_quantification2= CSV_quantification2.dropna(axis='index') 



# %% select the conditions you want to plot
CSV_singlecondition=CSV_quantification[(CSV_quantification.Condition=='control') 
                                    | (CSV_quantification.Condition=='lin14')
                                                                                ]

# oselect the worm you want to plot
wormtoplot=3
CSV_singleworm=CSV_quantification[CSV_quantification.Position==wormtoplot]

# %% plots whole dataset
plt.figure(1)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1, figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='final_intensity',hue='Condition',data=CSV_quantification,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume',hue='Condition',data=CSV_quantification,legend=False).set_title('volume')
axs[1].semilogy()

# %% plots whole dataset
plt.figure()
sns.set_theme(style="darkgrid")
sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=CSV_quantification,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume',hue='Condition',data=CSV_quantification,legend=False).set_title('volume')
#axs[1].semilogy()
plt.show()
# %% Plot individual worms or conditions
plt.figure(2)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1,figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='final_intensity',hue='Condition',data=CSV_singleworm,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume',hue='Condition',data=CSV_singleworm,legend=False).set_title('volume')
axs[1].semilogy()

plt.figure(3)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1,figsize=(10,15),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='final_intensity',hue='Condition',data=CSV_singlecondition,legend=True).set_title('final intensity')
sns.lineplot(ax=axs[1],x='Time',y='volume',hue='Condition',data=CSV_singlecondition,legend=False).set_title('volume')
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
