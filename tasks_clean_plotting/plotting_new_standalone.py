#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
from datetime import time
from numpy.core.defchararray import endswith, index
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas.core.frame import DataFrame
import seaborn as sns
import statsmodels.api as sm

def select_timepoints(frame,hatch,escape):
    if (frame>=hatch)and(frame<escape): return (frame-hatch)*min_image_taken/60
    else: return None

# Defining the function
def multiple_regression(data, dep_var, indep_var,  condition_variable = 'Condition_default', reference_conditon = 'Control_default', print_outliers = False):
    # Subsetting
    subset = data[[dep_var]+indep_var]

    # Check whether you run this for more than one condition
    if len(indep_var) == 1:
        subset[condition_variable] = reference_conditon
        indep_var = indep_var + [condition_variable]

    # Preparing
    subset['internal_noise'], subset['signal'] = [1., 0.]

    # Turning into dummy variables
    types = pd.unique(subset[condition_variable])
    # Selecting the reference
    if (not reference_conditon) or (reference_conditon == 0): type_vec = types[1:] # Set the first as reference
    else: type_vec = types[types!=reference_conditon]

    # Selecting the variables
    for type_sel in type_vec: subset[type_sel] = (subset[condition_variable] == type_sel).astype(float)
    subset_reg = subset.drop([condition_variable,'signal'], axis = 1)

    # Computing the sort
    result = sm.OLS(subset_reg.iloc[:,0], subset_reg.iloc[:,1:]).fit()

    # Storing the results
    # Creating
    param_df = subset[['internal_noise','signal']+indep_var].groupby(indep_var[1]).mean()
    #Parameters
    param_df[indep_var[0]] = result.params[indep_var[0]]
    param_df['internal_noise'] = result.params['internal_noise']
    # For each sample
    for type_sel in type_vec: param_df.loc[type_sel,'signal'] = result.params[type_sel]
    # Std
    param_df[indep_var[0]+'_stderr'] = result.bse[indep_var[0]]
    param_df['internal_noise_stderr'] = result.bse['internal_noise']
    for type_sel in type_vec: param_df.loc[type_sel,'signal_stderr'] = result.bse[type_sel]
    # Pvalues
    param_df[indep_var[0]+'_pval'] = result.pvalues[indep_var[0]]
    param_df['internal_noise_pval'] = result.pvalues['internal_noise']
    for type_sel in type_vec: param_df.loc[type_sel,'signal_pval'] = result.pvalues[type_sel]

    # Model statistics
    param_df['fval'] = result.f_pvalue
    param_df['r_squared'] = result.rsquared

    # Renaming
    param_df.rename(columns={"mean_background":"bg_factor"},inplace=True)
    param_df.rename(columns={"mean_background_stderr":"bg_factor_stderr"},inplace=True)
    param_df.rename(columns={"mean_background_pval":"bg_factor_pval"},inplace=True)

    # Give an alert of outliers
    if (result.outlier_test()['bonf(p)'].sum() != subset.shape[0]) and print_outliers:
        print('Outliers found at: ')
        print(data.loc[result.outlier_test()['bonf(p)']<1])

    return param_df

# %%
#  Import CSV files and merge datasets
dir='D:/GitHub/SWI_analysis/Analysis V3/'
goodworms=['goodworms_LIN14_LIN29.csv','goodworms_LIN28.csv',"goodworms_HBL1.csv", 'goodworms_LIN42BC.csv','goodworms_LIN41.csv','goodworms_LIN42AB.csv']
results=['Results_lin14_lin29_V3.csv','Results_lin28_V3.csv',"Results_HBL1_V3.csv",'Results_LIN42BC_V3.csv','Results_LIN41_V3.csv','Results_lin42AB_V3.csv']

min_image_taken=15


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
    CSV_quantification= CSV_quantification[CSV_quantification['Time']<42] #selects only timepoints smaller then 40hours

    #merges all dataset
    CSV_quantification_merged = pd.concat([CSV_quantification_merged, CSV_quantification])

#apply regression for all controls, and save results in: signal
controls=CSV_quantification_merged[CSV_quantification_merged["Condition"]=='control'] #get control worms to fit data on

Regression_results=controls.groupby('Time').apply(multiple_regression, 'mean_curved', ['mean_background']) #do linear regression per timepoint
Regression_results = Regression_results.droplevel(1)
CSV_quantification_merged=CSV_quantification_merged.merge(Regression_results,on=['Time'],how='inner') #merges quantification and annotation. 
CSV_quantification_merged["signal"]=CSV_quantification_merged["mean_curved"]-CSV_quantification_merged["bg_factor"]*CSV_quantification_merged["mean_background"]-CSV_quantification_merged["internal_noise"] #calculate signal

#get only the samples
Samples=CSV_quantification_merged[CSV_quantification_merged.Condition!='control'] #delete the control lines

#scale signal to 1/100 to get the: final intensity 
scale_to="signal" 
maximum=Samples.groupby(["Condition","Time",]).mean().groupby("Condition").max()[scale_to]
maximum=maximum.rename("maximum")    
Samples=Samples.merge(maximum,on=['Condition'], how='inner') #merges two datasets

minimum=Samples.groupby(["Condition","Time",]).mean().groupby("Condition").min()[scale_to]
mimimum=minimum.rename("minimum")    
Samples=Samples.merge(mimimum,on=['Condition'], how='inner') #merges two datasets

Samples["final intensity"] = (Samples[scale_to]-Samples['minimum'])/(Samples['maximum']-Samples['minimum'])*100



# %%
for y_plot in ['signal', 'final intensity']:
    plt.figure()
    sns.set_theme(style="whitegrid")
    snsplot = sns.lineplot(x='Time',y=y_plot,data=Samples, hue='Condition', err_style="band", legend=True)#.set_title("raw data")
    [plt.axvline(k,color='gray') for k in [12.5,21,30]]
    snsplot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.savefig("raw_data.pdf")


#%% Plotting the regression for an individual timepoint
tp = 20
controls_tp = CSV_quantification_merged[CSV_quantification_merged['Time']==tp]
controls_tp = controls_tp[controls_tp['Condition']=='control']

# Check all together by condition
dep_var = 'mean_curved'
condition_variable = 'old_condition'
reference_conditon =  'control_HBL1'
indep_var = ['mean_background',condition_variable]
regression_multiple = multiple_regression(controls_tp, dep_var, indep_var, condition_variable, reference_conditon)

# Check individually
indep_var = ['mean_background']
regression_controls_new = controls_tp.groupby('old_condition').apply(multiple_regression, dep_var, indep_var) #do linear regression per timepoint
regression_controls_new = regression_controls_new.droplevel(1)

# Check all together directly
regression_all_new = multiple_regression(controls_tp, dep_var, indep_var) 

#%% Plotting
sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="old_condition").set_title("Controls for timepoint "+str(tp)+ " (individually)")

x_fit = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)
y_fit = regression_controls_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_controls_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit)

y_fit = regression_all_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_all_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit,'k')

plt.figure()
sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="old_condition").set_title("Controls for timepoint "+str(tp)+ " (by condition)")

x_fit = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)
y_fit = regression_multiple["bg_factor"].values[None,:]*x_fit[:,None] + regression_multiple["internal_noise"].values[None,:] +  regression_multiple["signal"].values[None,:]
plt.plot(x_fit, y_fit)

y_fit = regression_all_new["bg_factor"].values[None,:]*x_fit[:,None] + regression_all_new["internal_noise"].values[None,:]
plt.plot(x_fit, y_fit,'k')

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
