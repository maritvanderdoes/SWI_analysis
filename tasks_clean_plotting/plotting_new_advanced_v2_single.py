#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from utils import molt_detection, molt_analysis
from utils import multiple_regression

def select_timepoints(frame,hatch,escape):
    if (frame>=hatch)and(frame<escape): return (frame-hatch)*min_image_taken/60
    else: return None

def removing_outliers(volume):
    return volume

# %% Import CSV files and merge datasets
dir='D:/GitHub/SWI_analysis/Analysis V3/'
goodworms=['goodworms_LIN29_NHR23i.csv']
results=['Results_NRH23i.csv']

min_image_taken=15


# %% Import CSV files and merge datasets
CSV_quantification_merged=pd.DataFrame()
control_only_merged=pd.DataFrame()
nr_samples=np.zeros(np.size(goodworms))

for i in range(np.size(goodworms)):
    CSV_quantification = pd.read_csv(dir+results[i])
    CSV_annotation= pd.read_csv(dir+goodworms[i])

    CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
    CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1

    #Add time
    CSV_quantification['Time'] = CSV_quantification.apply(lambda row : select_timepoints(row['Frame'],row['Hatch'],row['Escape']), axis = 1) #add time column
    CSV_quantification= CSV_quantification.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape
    # CSV_quantification= CSV_quantification[CSV_quantification['Time']<42] #selects only timepoints smaller then 40hours

    #merges all dataset
    CSV_quantification_merged = pd.concat([CSV_quantification_merged, CSV_quantification])

# REmove lin-14
# CSV_quantification_merged = CSV_quantification_merged[CSV_quantification_merged["Condition"]!='LIN-14']

#apply regression for all controls, and save results in: signal
controls=CSV_quantification_merged[CSV_quantification_merged["Condition"]=='control'] #get control worms to fit data on


#%%
sample_sel = CSV_quantification_merged[\
            (CSV_quantification_merged['Position']==10)&\
            (CSV_quantification_merged['Condition']=='LIN-29')].sort_values(by='Time')

CSV_sing = CSV_quantification_merged.copy()
CSV_sing.loc[CSV_sing['old_condition'].isnull(),'old_condition'] = CSV_sing.loc[CSV_sing['old_condition'].isnull(),'Condition']+'_smpl'
print(pd.unique(CSV_sing['old_condition']))

CSV_sing = CSV_sing.groupby(['old_condition','Position','Time']).mean().reset_index(2)

# %% Volume processing
# print(pd.unique(CSV_sing['Condition']))
from  scipy.ndimage import median_filter

median_volumes = CSV_quantification_merged.reset_index().groupby(['Condition','Time']).max()
median_volumes = CSV_quantification_merged.reset_index().groupby(['Time']).max()

th = 0.8
sample = 97
chamber_idx_vec = pd.unique(CSV_quantification_merged.index)
chamber_volume = CSV_quantification_merged.loc[chamber_idx_vec[sample],['Time','volume_curved']].reset_index().groupby('Time').mean()

chamber_volume['diff_vol'] = chamber_volume['volume_curved'] - median_volumes.loc[chamber_idx_vec[50][0],'volume_curved']
chamber_volume['volume_curved_log'] = np.log(chamber_volume['volume_curved'])
chamber_volume['diff_vol_log'] = np.log(chamber_volume['volume_curved']) - np.log(median_volumes.loc[chamber_idx_vec[50][0],'volume_curved'])
chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_volume'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['volume_curved']

g_results = sns.lineplot(data = median_volumes['volume_curved'].reset_index(), x = 'Time', y='volume_curved', hue = 'old_condition')
g_results.set(yscale='log')

CSV_quantification_merged.loc[chamber_idx_vec[sample],'actual_volume'] = chamber_volume['actual_volume']

plt.figure()
plt.semilogy(chamber_volume['volume_curved'])
plt.semilogy(chamber_volume['actual_volume'],'-o')

# chamber_volume.apply(median_filter, size = 3)
plt.figure()
plt.semilogy(chamber_volume['volume_curved'])
# plt.semilogy(chamber_volume[chamber_volume['diff_vol']<-5e4]['volume_curved'],'rx')
plt.semilogy(chamber_volume[chamber_volume['diff_vol_log']<-th]['volume_curved'],'rx')
# plt.semilogy(chamber_median)
plt.semilogy(median_volumes.loc[chamber_idx_vec[sample][0],'volume_curved'])

# plt.figure()
# plt.plot(chamber_volume['diff_vol'])
# plt.plot(chamber_volume[chamber_volume['diff_vol']<-5e4]['diff_vol'],'rx')

plt.figure()
plt.plot(chamber_volume['diff_vol_log'])
plt.plot(chamber_volume[chamber_volume['diff_vol_log']<-th]['diff_vol_log'],'x')

# plt.plot(chamber_volume[chamber_volume['diff_vol']<-5e4]['diff_vol'],'rx')
#%% Iterate
CSV_quantification_merged = CSV_quantification_merged.groupby(['Condition','Position','Time']).mean().reset_index(2)
median_volumes = CSV_quantification_merged.reset_index().groupby(['Condition','Time']).max()
CSV_processed = CSV_quantification_merged.copy()
th = 0.8
sample = 97
chamber_idx_vec = pd.unique(CSV_quantification_merged.index)
for sample, chamber_idx in enumerate(pd.unique(CSV_quantification_merged.index)):
    chamber_volume = CSV_quantification_merged.loc[chamber_idx,['Time','volume_curved','mean_curved', 'mean_background']].reset_index().groupby('Time').mean()

    chamber_volume['diff_vol'] = chamber_volume['volume_curved'] - median_volumes.loc[chamber_idx[0],'volume_curved']
    chamber_volume['volume_curved_log'] = np.log(chamber_volume['volume_curved'])
    chamber_volume['diff_vol_log'] = np.log(chamber_volume['volume_curved']) - np.log(median_volumes.loc[chamber_idx[0],'volume_curved'])
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_volume'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['volume_curved']
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_mean'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['mean_curved']
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_background'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['mean_background']

    CSV_processed.loc[chamber_idx,'volume_curved'] = chamber_volume['actual_volume'].values
    CSV_processed.loc[chamber_idx,'mean_curved'] = chamber_volume['actual_mean'].values
    CSV_processed.loc[chamber_idx,'mean_background'] = chamber_volume['actual_background'].values
    plt.semilogy(chamber_volume['actual_volume'].values)

CSV_processed = CSV_processed.reset_index().dropna()

print(CSV_processed['volume_curved'].isnull().sum())

#%% Iterate
CSV_quantification_merged = CSV_quantification_merged.groupby(['Condition','Position','Time']).mean().reset_index(2)
median_volumes = CSV_quantification_merged.reset_index().groupby(['Time']).max()
CSV_processed = CSV_quantification_merged.copy()
th = 0.8
sample = 97
chamber_idx_vec = pd.unique(CSV_quantification_merged.index)
for sample, chamber_idx in enumerate(pd.unique(CSV_quantification_merged.index)):
    chamber_volume = CSV_quantification_merged.loc[chamber_idx,['Time','volume_curved','mean_curved', 'mean_background']].reset_index().groupby('Time').mean()

    chamber_volume['diff_vol'] = chamber_volume['volume_curved'] - median_volumes['volume_curved']
    chamber_volume['volume_curved_log'] = np.log(chamber_volume['volume_curved'])
    chamber_volume['diff_vol_log'] = np.log(chamber_volume['volume_curved']) - np.log(median_volumes['volume_curved'])
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_volume'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['volume_curved']
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_mean'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['mean_curved']
    chamber_volume.loc[chamber_volume['diff_vol_log']>-th,'actual_background'] = chamber_volume[chamber_volume['diff_vol_log']>-th]['mean_background']

    CSV_processed.loc[chamber_idx,'volume_curved'] = chamber_volume['actual_volume'].values
    CSV_processed.loc[chamber_idx,'mean_curved'] = chamber_volume['actual_mean'].values
    CSV_processed.loc[chamber_idx,'mean_background'] = chamber_volume['actual_background'].values
    plt.semilogy(chamber_volume['actual_volume'].values)

CSV_processed = CSV_processed.reset_index().dropna()

print(CSV_processed['volume_curved'].isnull().sum())

#%%
# dynamics_to_plot = CSV_processed.reset_index().groupby(['Condition','Time']).mean().reset_index()
nhr23_df = CSV_processed.reset_index().groupby(['Condition','Time','Position']).mean().reset_index()
sns.lineplot(data = nhr23_df, x = 'Time', y = 'mean_curved', hue = 'Condition')

# nhr23_df = dynamics_to_plot.loc[dynamics_to_plot['Condition'].isin(['LIN29_GFP_smpl','LIN42_GFP_smpl','control_smpl'])]
# nhr23_df = dynamics_to_plot.loc[np.invert(dynamics_to_plot['old_condition'].isin(['LIN29_GFP_smpl','LIN42_GFP_smpl','control_smpl']))]
nhr23_df = nhr23_df[['Condition','Time','Position','mean_curved','mean_background','volume_curved']].dropna()
plt.figure()
sns.lineplot(data = nhr23_df, x = 'Time', y = 'mean_curved', hue = 'Condition')

plt.figure()
sns.lineplot(data = nhr23_df, x = 'Time', y = 'mean_background', hue = 'Condition')

#%%
nhr23_df['QuickComp'] = nhr23_df['mean_curved']-nhr23_df['mean_background']
sns.lineplot(data = nhr23_df, x = 'Time', y = 'QuickComp', hue = 'Condition')


#%%
dep_var = 'mean_curved'
condition_variable = 'Condition'
reference_conditon =  'control'
reference_conditon =  'LIN42_GFP'
indep_var = ['mean_background',condition_variable]

Regression_results=nhr23_df.groupby('Time').apply(multiple_regression, dep_var, indep_var, condition_variable, reference_conditon) #do linear regression per timepoint

nhr23_df = nhr23_df.merge(Regression_results.reset_index(), on = ['Time', 'Condition'])
nhr23_df['Background'] = nhr23_df['bg_factor']*nhr23_df['mean_background'] + nhr23_df['internal_noise']
nhr23_df['signal_ind'] = nhr23_df['mean_curved'] - nhr23_df['Background'] 

#%%
sns.lineplot(data = nhr23_df, x= 'Time', y = 'signal_ind', hue = 'Condition')
# plt.ylim([-2,30])

#%%
comp = nhr23_df.groupby(['Condition','Time']).mean().reset_index()
comp[comp['Condition'].isin(['LIN29_GFP','control'])].groupby('Condition').diff()
plt.plot((comp[comp['Condition']=='LIN29_GFP']['signal_ind']-comp[comp['Condition']=='control']['signal_ind']).values)


#%%
snsplot = sns.lineplot(data = nhr23_df, x= 'Time', y = 'volume_curved', hue = 'Condition')
snsplot.set(yscale='log')

#%%
sns.lineplot(data = Regression_results.reset_index(), x= 'Time', y = 'signal', hue = 'Condition')

#%%
sns.lineplot(data = Regression_results.reset_index(), x= 'Time', y = 'signal_pval', hue = 'Condition')
plt.axhline(0.05,color='gray')

#%%
sns.lineplot(data = Regression_results.reset_index(), x= 'Time', y = 'r_squared', hue = 'Condition')

#%%
sns.lineplot(data = Regression_results[Regression_results['r_squared']>0.55].reset_index(), x= 'Time', y = 'signal', hue = 'Condition')

#%%
controls = nhr23_df[nhr23_df['old_condition']=='control_smpl']
#%%
Regression_results=controls.groupby('Time').apply(multiple_regression, 'mean_curved', ['mean_background']) #do linear regression per timepoint
Regression_results = Regression_results.droplevel(1)
nhr23_df=nhr23_df.merge(Regression_results,on=['Time'],how='inner') #merges quantification and annotation. 
nhr23_df["signal"]=nhr23_df["mean_curved"]-CSV_quantification_merged["bg_factor"]*CSV_quantification_merged["mean_background"]-CSV_quantification_merged["internal_noise"] #calculate signal

#get only the samples
Samples=nhr23_df[nhr23_df.old_condition!='control_smpl'] #delete the control lines

#scale signal to 1/100 to get the: final intensity 
scale_to="signal" 
maximum=Samples.groupby(["old_condition","Time",]).mean().groupby("old_condition").max()[scale_to]
maximum=maximum.rename("maximum")    
Samples=Samples.merge(maximum,on=['old_condition'], how='inner') #merges two datasets

minimum=Samples.groupby(["old_condition","Time",]).mean().groupby("old_condition").min()[scale_to]
mimimum=minimum.rename("minimum")    
Samples=Samples.merge(mimimum,on=['old_condition'], how='inner') #merges two datasets

Samples["final intensity"] = (Samples[scale_to]-Samples['minimum'])/(Samples['maximum']-Samples['minimum'])*100

# %%
sns.set_theme(style="whitegrid")
for y_plot in ['signal', 'final intensity']:
    plt.figure()
    snsplot = sns.lineplot(x='Time',y=y_plot,data=Samples, hue='old_condition', err_style="band", legend=True)
    snsplot.set_title("raw data")
    snsplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    [plt.axvline(k,color='gray') for k in [12.5,21,30]]
    #plt.savefig("raw_data.pdf")

# %% Individual plots
condition_vec = ['LIN-29', 'LIN-28', 'LIN42BC', 'LIN-14', 'LIN-41', 'HBL-1']
fig, axs = plt.subplots(len(condition_vec), 1, figsize=(5,18),sharex=True)

sns.set_theme(style="whitegrid")
for k_c, condition_sel in enumerate(condition_vec):
    toplot=Samples[Samples["Condition"]==condition_sel]
    sns.lineplot(ax=axs[k_c],x='Time',y="signal",data=toplot).set_title(condition_sel)
    [axs[k_c].axvline(k,color='gray') for k in [12.5,21,30]]

plt.tight_layout()
plt.savefig('Individual_tracs_perexp.pdf')

#%% Plotting the regression for an individual timepoint
tp = 22
controls_tp = nhr23_df[nhr23_df['Time']==tp]
# controls_tp = controls_tp[controls_tp['old_condition']=='control_smpl']

sns.regplot(data = controls_tp,x = 'mean_background', y = 'mean_curved')

# Check all together by condition
dep_var = 'mean_curved'
condition_variable = 'Condition'
reference_conditon =  'control_smpl'
indep_var = ['mean_background',condition_variable]
regression_multiple = multiple_regression(controls_tp, dep_var, indep_var, condition_variable, reference_conditon)

# Check individually
indep_var = ['mean_background']
regression_controls_new = controls_tp.groupby('Condition').apply(multiple_regression, dep_var, indep_var) #do linear regression per timepoint
regression_controls_new = regression_controls_new.droplevel(1)

# Check all together directly
regression_all_new = multiple_regression(controls_tp, dep_var, indep_var) 

#% Plotting
plt.figure()
sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="Condition").set_title("Controls for timepoint "+str(tp)+ " (individually)")

x_fit = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)
y_fit = regression_controls_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_controls_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit)

y_fit = regression_all_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_all_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit,'k')

plt.figure()
sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="Condition").set_title("Controls for timepoint "+str(tp)+ " (by condition)")

x_fit = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)
y_fit = regression_multiple["bg_factor"].values[None,:]*x_fit[:,None] + regression_multiple["internal_noise"].values[None,:] +  regression_multiple["signal"].values[None,:]
plt.plot(x_fit, y_fit)

y_fit = regression_all_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_all_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit,'k')

regression_controls_new

#%% Checking per sample
fig, axs = plt.subplots(pd.unique(controls_tp['old_condition']).shape[0],1,\
            figsize = (10,30), sharex = False)

for k_c, condition in enumerate(pd.unique(controls_tp['old_condition'])):
    subset = controls_tp[controls_tp['old_condition']==condition]
    # sns.regplot(data=subset,x="mean_background", y="mean_curved", ax=axs[k_c])
    sns.scatterplot(data=subset,x="mean_background", y="mean_curved", hue = 'old_condition', ax=axs[k_c])

    x_fit = np.linspace(np.min(subset["mean_background"])-10,np.max(subset["mean_background"])+10,100)
    y_fit = regression_multiple.loc[condition,"bg_factor"]*x_fit[:,None] + regression_multiple.loc[condition,"internal_noise"]+  regression_multiple.loc[condition,"signal"]
    axs[k_c].plot(x_fit, y_fit, 'b:')

    y_fit = regression_controls_new.loc[condition,"bg_factor"]*x_fit[:,None] + regression_controls_new.loc[condition,"internal_noise"]
    axs[k_c].plot(x_fit, y_fit,'c:')

    y_fit = regression_all_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_all_new["internal_noise"].values[None,:]
    axs[k_c].plot(x_fit, y_fit,'k')


# %%
dep_var = 'mean_curved'
condition_variable = 'old_condition'
reference_conditon =  'control_HBL1'
reference_conditon =  'control_Lin28'
indep_var = ['mean_background',condition_variable]

Regression_results=controls.groupby('Time').apply(multiple_regression, dep_var, indep_var, condition_variable, reference_conditon) #do linear regression per timepoint

# %%
snsplot = sns.lineplot(data = Regression_results.reset_index(), y = 'signal', x = 'Time', hue = 'old_condition')
snsplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
[plt.axvline(k,color='gray') for k in [12.5,21,30]];


# %% COMPARING BETWEEN MULTIPLE AND INDIVIDUAL REGRESSIONS
dep_var = 'mean_curved'
condition_variable = 'Condition'
reference_conditon =  'control'
indep_var = ['mean_background',condition_variable]

Regression_results = CSV_quantification_merged.groupby('Time').apply(multiple_regression, dep_var, indep_var, condition_variable, reference_conditon) #do linear regression per timepoint
Independent_results = Samples.groupby(['Time','Condition']).mean()

# %%
snsplot = sns.lineplot(data = Regression_results.reset_index(), y = 'signal', x = 'Time', hue = 'Condition')
snsplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
[plt.axvline(k,color='gray') for k in [12.5,21,30]];

plt.figure()
snsplot = sns.lineplot(data = Independent_results.reset_index(), y = 'signal', x = 'Time', hue = 'Condition')
snsplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
[plt.axvline(k,color='gray') for k in [12.5,21,30]];

# Plotting correlations in a single condition
mult_reg = Regression_results.reset_index()
indp_reg = Independent_results.reset_index()
cndt = 'HBL-1'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt]['signal'],
         indp_reg[indp_reg['Condition']==cndt]['signal'],'x')
plt.xlabel('Multiple regression')
plt.ylabel('Independent regression')

# Plotting correlations in all conditions
plt.figure()
for cndt in pd.unique(indp_reg['Condition']):
    plt.plot(mult_reg[mult_reg['Condition']==cndt]['signal'],
             indp_reg[indp_reg['Condition']==cndt]['signal'],'x')
plt.legend(pd.unique(indp_reg['Condition']))
plt.xlabel('Multiple regression')
plt.ylabel('Independent regression')
# %% Phase spaces
w_size = 1*4
mult_reg = Regression_results.reset_index()[['Time', 'Condition', 'signal']]
time_reg = pd.unique(mult_reg['Time'])
cndt_A = 'HBL-1'
cndt_B = 'LIN-14'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean(),mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean(),'-o')
plt.xlabel(cndt_A)
plt.ylabel(cndt_B)

plt.figure()
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean())
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean())
plt.xlabel('Time')
plt.ylabel('Levels')

cndt_A = 'LIN-28'
cndt_B = 'LIN-14'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean(),mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean(),'-o')
plt.xlabel(cndt_A)
plt.ylabel(cndt_B)

plt.figure()
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean())
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean())
plt.xlabel('Time')
plt.ylabel('Levels')

cndt_A = 'LIN-28'
cndt_B = 'LIN-29'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean(),mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean(),'-o')
plt.xlabel(cndt_A)
plt.ylabel(cndt_B)

plt.figure()
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean())
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean())
plt.xlabel('Time')
plt.ylabel('Levels')


cndt_A = 'LIN-29'
cndt_B = 'LIN-14'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean(),mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean(),'-o')
plt.xlabel(cndt_A)
plt.ylabel(cndt_B)

plt.figure()
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean())
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean())
plt.xlabel('Time')
plt.ylabel('Levels')

cndt_A = 'HBL-1'
cndt_B = 'LIN-29'
plt.figure()
plt.plot(mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean(),mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean(),'-o')
plt.xlabel(cndt_A)
plt.ylabel(cndt_B)

plt.figure()
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_A]['signal'].rolling(w_size).mean())
plt.plot(time_reg,mult_reg[mult_reg['Condition']==cndt_B]['signal'].rolling(w_size).mean())
plt.xlabel('Time')
plt.ylabel('Levels')
# %%
