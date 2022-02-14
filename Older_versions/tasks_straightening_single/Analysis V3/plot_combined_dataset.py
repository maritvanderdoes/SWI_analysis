#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% import packages
from numpy.core.defchararray import endswith, index
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pandas.core.frame import DataFrame
import seaborn as sns
import statsmodels.api as sm

def select_timepoints(frame,hatch,escape):
    if (frame>=hatch)and(frame<escape):
        return (frame-hatch)*min_image_taken/60
    else:
        return None

def regress(data, yvar, xvars):
    Y = data[yvar]
    X = data[xvars]    
    X['internal_noise'] = 1.
    X.rename(columns={"mean_background":"bg_factor"},inplace=True)
    result = sm.OLS(Y, X).fit()
    return result.params


# %%
#  Import CSV files and merge datasets
dir='/Volumes/ggrossha/Marit/HetP_Quant/Analysis V3/'
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

Regression_results=controls.groupby('Time').apply(regress, 'mean_curved', ['mean_background']) #do linear regression per timepoint
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
#check if controls are similar
tp=20
controls_tp=controls[controls["Time"]==tp]

#plot dots plus regression
plt.figure()
sns.scatterplot(data=controls_tp,x="mean_background", y="mean_curved",hue="old_condition").set_title("controls for tp "+str(tp))
sns.scatterplot(data=Samples[(Samples["Condition"]=="HBL-1")&(Samples["Time"]==tp)],x="mean_background",y="mean_curved",color="black")


Regression_controls=controls_tp.groupby('old_condition').apply(regress, 'mean_curved', ['mean_background']) #do linear regression per timepoint
x = np.linspace(np.min(controls_tp["mean_background"])-10,np.max(controls_tp["mean_background"])+10,100)

for i in range(np.size(Regression_controls,0)):
    y= x*Regression_controls["bg_factor"].values[i]+Regression_controls["internal_noise"].values[i]
    plt.plot(x,y)

Regression_all=regress(controls_tp,'mean_curved', ['mean_background'])
y= x*Regression_all["bg_factor"]+Regression_all["internal_noise"]
plt.plot(x,y,color="black")

controls_to_plot=CSV_quantification_merged[CSV_quantification_merged["Condition"]=="control"]
plt.figure()
sns.lineplot(data=controls_to_plot,x="Time",y="signal",hue="old_condition").set_title("control_experiments together")



# %% 
#plot all samples together
plt.figure(1)
sns.set_theme(style="whitegrid")
sns.lineplot(x='Time',y="signal",data=Samples, hue='Condition', err_style="band", legend=True).set_title("raw data")
plt.axvline(12.5,color='gray')
plt.axvline(21,color='gray')
plt.axvline(30,color='gray')
#plt.savefig("raw_data.pdf")

plt.figure(2)
sns.set_theme(style="whitegrid")
sns.lineplot(x='Time',y="final intensity",data=Samples, hue='Condition', err_style="band", legend=True).set_title("raw data")
plt.axvline(12.5,color='gray')
plt.axvline(21,color='gray')
plt.axvline(30,color='gray')
#plt.savefig("raw_data.pdf")


# %% fig 2 plots individal experiments
fig, axs = plt.subplots(6, 1, figsize=(5,18),sharex=True)
sns.set_theme(style="whitegrid")
toplot=Samples[Samples["Condition"]=="LIN-29"]
sns.lineplot(ax=axs[0],x='Time',y="signal",data=toplot).set_title("LIN-29")
axs[0].axvline(12.5,color='gray')
axs[0].axvline(21,color='gray')
axs[0].axvline(30,color='gray')

toplot=Samples[Samples["Condition"]=="LIN-28"]
sns.lineplot(ax=axs[1],x='Time',y="signal",data=toplot).set_title("LIN-28")
axs[1].axvline(12.5,color='gray')
axs[1].axvline(21,color='gray')
axs[1].axvline(30,color='gray')

toplot=Samples[(Samples["Condition"]=="LIN42BC")| (Samples["Condition"]=="lin42_ab")]
sns.lineplot(ax=axs[2],x='Time',y="signal",data=toplot,hue="Condition").set_title("LIN-42")
axs[2].axvline(12.5,color='gray')
axs[2].axvline(21,color='gray')
axs[2].axvline(30,color='gray')

toplot=Samples[Samples["Condition"]=="LIN-14"]
sns.lineplot(ax=axs[3],x='Time',y="signal",data=toplot).set_title("LIN-14")
axs[3].axvline(12.5,color='gray')
axs[3].axvline(21,color='gray')
axs[3].axvline(30,color='gray')

toplot=Samples[Samples["Condition"]=="LIN-41"]
sns.lineplot(ax=axs[4],x='Time',y="signal",data=toplot).set_title("LIN-41")
axs[4].axvline(12.5,color='gray')
axs[4].axvline(21,color='gray')
axs[4].axvline(30,color='gray')

toplot=Samples[Samples["Condition"]=="HBL-1"]
sns.lineplot(ax=axs[5],x='Time',y="signal",data=toplot).set_title("HBL-1")
axs[5].axvline(12.5,color='gray')
axs[5].axvline(21,color='gray')
axs[5].axvline(30,color='gray')

plt.tight_layout()
plt.savefig('Individual_tracs_alltogether.pdf')

# %% mRNA adding
filename =  "/Volumes/ggrossha/Lucas/RNAseq_Datasets/TL_all_fused.csv"
filename_miRNA =  "/Volumes/ggrossha/Lucas/RNAseq_Datasets/miRNA_profiles.csv"

#reading filename
mRNA= pd.read_csv(filename,index_col=0)
mRNA.columns = np.arange(len(mRNA.columns))+1
mRNA= pd.melt(mRNA,ignore_index=False ,var_name="Time")

miRNA=pd.read_csv(filename_miRNA,index_col=0)
miRNA.columns = np.arange(len(miRNA.columns))+5
miRNA= pd.melt(miRNA,ignore_index=False ,var_name="Time")

RNA=pd.concat([mRNA,miRNA])

#add pseudotime
molts_SWI = np.array([12.5, 21, 30, 42])
molts_RNA = np.array([12, 18, 25, 35])
scal_factor = np.mean(molts_SWI/molts_RNA)
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


# %% 
# fig 3 plotting genes with dataset
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


# %%
