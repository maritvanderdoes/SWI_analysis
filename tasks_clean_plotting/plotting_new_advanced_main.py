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


# %% Import CSV files and merge datasets
dir='D:/GitHub/SWI_analysis/Analysis V3/'
# goodworms=['goodworms_LIN14_LIN29.csv','goodworms_LIN28.csv',"goodworms_HBL1.csv",'goodworms_LIN41.csv','goodworms_LIN42BC.csv']
# results=['Results_lin14_lin29_V3.csv','Results_lin28_V3.csv',"Results_HBL1_V3.csv",'Results_LIN41_V3.csv','Results_lin42BC_V3.csv']

# goodworms=['goodworms_LIN28.csv']
# results=['Results_lin28_V3.csv']

goodworms=['goodworms_HBL1.csv']
results=['Results_HBL1_V3.csv']


# mRNA adding
filename_rna =  "D:/Datasets/TL_all_fused.csv"
filename_miRNA =  "D:/Datasets/miRNA_profiles.csv"

#add pseudotime
molts_SWI = np.array([12.5, 21, 30, 42])
molts_RNA = np.array([12, 18, 25, 35])
scal_factor = np.mean(molts_SWI/molts_RNA)

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
    CSV_quantification= CSV_quantification[CSV_quantification['Time']<42] #selects only timepoints smaller then 40hours

    #merges all dataset
    CSV_quantification_merged = pd.concat([CSV_quantification_merged, CSV_quantification])

#apply regression for all controls, and save results in: signal
controls=CSV_quantification_merged[CSV_quantification_merged["Condition"]=='control'] #get control worms to fit data on

# %% Molt aligment (V1)

# %% Quantification
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
sns.set_theme(style="whitegrid")
for y_plot in ['signal', 'final intensity']:
    plt.figure()
    snsplot = sns.lineplot(x='Time',y=y_plot,data=Samples, hue='Condition', err_style="band", legend=True)
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

#%%

#reading filename
mRNA= pd.read_csv(filename_rna,index_col=0)
mRNA.columns = np.arange(len(mRNA.columns))+1
mRNA= pd.melt(mRNA,ignore_index=False ,var_name="Time")

miRNA=pd.read_csv(filename_miRNA,index_col=0)
miRNA.columns = np.arange(len(miRNA.columns))+5
miRNA= pd.melt(miRNA,ignore_index=False ,var_name="Time")

RNA=pd.concat([mRNA,miRNA])
RNA["Time_scaled"]=RNA["Time"]*scal_factor

#select certain times, and scale
RNA= RNA[RNA['Time_scaled']<42]

#make RNA in linear scale
RNA["linear"]=2**RNA["value"]

#add genename
genename={"Samples":['WBGene00003003','WBGene00003014','WBGene00003015','WBGene00003026','WBGene00018572','WBGene00001824','cel-lin-4-5p','cel-lin-4-3p','cel-let-7-5p','cel-let-7-3p','cel-miR-48-5p','cel-miR-48-3p','cel-miR-84-5p','cel-miR-84-3p','cel-miR-241-5p','cel-miR-241-3p'],
          "Name":['LIN-14','LIN-28','LIN-29','LIN-41','LIN42BC','HBL-1', 'LIN-4 guide','LIN-4 passenger','LET-7 guide', 'LET-7 passenger', 'miR-48 guide','miR-48 passenger','miR-84 guide','miR-84 passenger','miR-241 guide','miR-241 passenger']}
genename=pd.DataFrame(genename).set_index("Samples")
RNA=RNA.merge(genename, how="inner",left_index=True,right_index=True)

#scale to maximum and minimum
scale_to="linear"
maximum=RNA.groupby("Name").max()[scale_to]
maximum=maximum.rename("maximum")  
minimum=RNA.groupby("Name").min()[scale_to]
minimum=minimum.rename("minimum")  
RNA=RNA.merge(maximum,on=['Name'], how='inner') #merges two datasets
RNA=RNA.merge(minimum,on=['Name'], how='inner') #merges two datasets
RNA["final intensity"] = (RNA[scale_to]-RNA['minimum'])/(RNA['maximum']-RNA['minimum'])*100

#%%
fig, axs = plt.subplots(6, 2, figsize=(10,18),sharex=True,sharey=False)
sns.set_theme(style="whitegrid")

#plotting proteins
for k,name in enumerate(genename["Name"][0:6]):

    toplot_protein=Samples[Samples["Condition"]==name]
    sns.lineplot(ax=axs[k,0],x='Time',y="final intensity",data=toplot_protein).set_title(name)
    toplot_mRNA=RNA[RNA["Name"]==name]
    sns.lineplot(ax=axs[k,0],x='Time_scaled',y="final intensity",data=toplot_mRNA)
    axs[k,0].legend(["Protein","mRNA"])

    axs[k,0].axvline(12.5,color='gray')
    axs[k,0].axvline(21,color='gray')
    axs[k,0].axvline(30,color='gray')


title=["Lin-4","Let-7",'MiR-48',"MiR-84","miR-241"]
for k,name in enumerate(genename["Name"][6::]):

    toplot_miRNA=RNA[RNA["Name"]==name]
    sns.lineplot(ax=axs[int(k/2),1],x='Time',y="value",data=toplot_miRNA).set_title(title[int(k/2)])

    
    if np.mod(k,2)==1:
        axs[int(k/2),1].legend(["miRNA-guide","miRNA-passenger"])
    
        axs[int(k/2),1].axvline(12.5,color='gray')
        axs[int(k/2),1].axvline(21,color='gray')
        axs[int(k/2),1].axvline(30,color='gray')

toplot_protein=Samples[Samples["Condition"]=="lin42_ab"]
sns.lineplot(ax=axs[4,0], x='Time',y="final intensity",data=toplot_protein).set_title("LIN-42")


plt.savefig("With_mRNA.pdf")


###############################################
## ANALYSIS REGRESSION ########################
###############################################

#%% Plotting the regression for an individual timepoint
tp = 5
controls_tp = CSV_quantification_merged[CSV_quantification_merged['Time']==tp]
# controls_tp = controls_tp[controls_tp['Condition']=='control']

# Check all together by condition
dep_var = 'mean_curved'
# condition_variable = 'old_condition'
# reference_conditon =  'control_HBL1'
condition_variable = 'Condition'
reference_conditon =  'control'
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

#%%
plt.figure(figsize=(6,4))
controls_tp = controls_tp.sort_values(by = 'Condition')

sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="Condition").set_title("Regression for timepoint "+str(tp)+ " ")

controls_tp

x_fit = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)
y_fit = regression_multiple["bg_factor"].values[None,:]*x_fit[:,None] + regression_multiple["internal_noise"].values[None,:] +  regression_multiple["signal"].values[None,:]
plt.plot(x_fit, y_fit)


plt.savefig('Figure_HBL1_5.pdf', dpi = 200)


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
