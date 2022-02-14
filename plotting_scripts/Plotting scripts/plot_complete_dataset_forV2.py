#%% import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#%% import files and merge together
CSV_quantification = pd.read_csv("/Volumes/ggrossha/Marit/HetP_Quant/20210507_lin28_analysis/V2_analysis/Results_parallel_lin28.csv")
CSV_annotation= pd.read_csv("/Volumes/ggrossha/Marit/HetP_Quant/20210507_lin28_analysis/V2_analysis/goodworms_lin28.csv")
min_image_taken=15 #time interval in minutes
#color_graphs=["blue","orange"]
include_molt_exit=False
plot_till= "Escape"
saving_directory='/Users/Marit/Documents/'

#%% 
# Basic merging of files
def add_time(frame,start,end):
    if (frame>=start)and(frame<end):
        return (frame-start)*min_image_taken/60
    else:
        return None

def molt_annotation(hatch,molt):
    return (molt-hatch)*min_image_taken/60


#merge arrays for complete dataset
CSV_quantification.head()
CSV_annotation.head()

CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1
CSV_quantification['Time'] = CSV_quantification.apply(lambda row : add_time(row['Frame'],row['Hatch'],row['Escape']), axis = 1) #add overall time
CSV_quantification= CSV_quantification.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape


#extract conditions and number of worms
conditions=pd.unique(CSV_annotation["Condition"])
goodworms=CSV_annotation[CSV_annotation['Quality']==1] #selects only worms where quality is 1
number_of_worms=np.empty(np.shape(conditions))
for i,condition in enumerate(conditions):
    number_of_worms[i]=int(np.sum(goodworms["Condition"]==condition))


#add timing molts
if include_molt_exit==True:
    CSV_quantification['TimeMolt1'] = CSV_quantification.apply(lambda row : molt_annotation(row['Hatch'],row['M1']), axis = 1) #set timing based on molt1
    CSV_quantification['TimeMolt2'] = CSV_quantification.apply(lambda row : molt_annotation(row['Hatch'],row['M2']), axis = 1) #set timing based on molt1
    CSV_quantification['TimeMolt3'] = CSV_quantification.apply(lambda row : molt_annotation(row['Hatch'],row['M3']), axis = 1) #set timing based on molt1
    CSV_quantification['TimeMolt4'] = CSV_quantification.apply(lambda row : molt_annotation(row['Hatch'],row['M4']), axis = 1) #set timing based on molt1


#quick plot to check
sns.lineplot(x='Time',y='mean_curved',hue='Condition',data=CSV_quantification,ci="sd").set_title("quick plot")


#%% 
# Fig 1 Plot individual worms
def filter_data(quantification_df,begin,end,condition):

    #create arrays to plot for complete dataset
    quantification_df['Time_temp'] = quantification_df.apply(lambda row : add_time(row['Frame'],row[begin],row[end]), axis = 1) #add time column
    quantification_df= quantification_df.dropna(subset=['Time_temp'],axis='index') #selects only worms in between begin and and

    track_df=quantification_df[quantification_df["Condition"]==condition].pivot(index="Time",columns="Position",values="mean_curved")
    track=track_df.to_numpy()
    track_index=track_df.index.to_numpy()
    position=track_df.columns.to_numpy()

    if include_molt_exit==True:
        molting_time=quantification_df[quantification_df["Condition"]==condition].groupby("Position").mean()
        TM1=molting_time["TimeMolt1"].to_numpy()
        TM2=molting_time["TimeMolt2"].to_numpy()
        TM3=molting_time["TimeMolt3"].to_numpy()
        TM4=molting_time["TimeMolt4"].to_numpy()
    else:
        TM1=0
        TM2=0
        TM3=0
        TM4=0


    return track,track_index,position,TM1,TM2,TM3,TM4


#make figure
fig, axs = plt.subplots(int(np.max(number_of_worms)),np.size(conditions), figsize=(15,80),sharex=True)

for column,condition in enumerate(conditions):
    [all, idx_all,position,TM1,TM2,TM3,TM4]=filter_data(CSV_quantification,'Hatch',plot_till,condition) #all

    for plotpos in range(np.shape(all)[1]):
        axs[plotpos,column].plot(idx_all,all[:,plotpos],color=color_graphs[column], linewidth=1)
        axs[plotpos,column].plot(idx_all,np.mean(all,1),color="black",linestyle='dashed',linewidth=1) #average

        if include_molt_exit==True:
            axs[plotpos,column].axvline(x=TM1[plotpos], color="black",linestyle=":")
            axs[plotpos,column].axvline(x=TM2[plotpos], color="black",linestyle=":")
            axs[plotpos,column].axvline(x=TM3[plotpos], color="black",linestyle=":")
            axs[plotpos,column].axvline(x=TM4[plotpos], color="black",linestyle=":")

        axs[plotpos,column].title.set_text(condition+"_position "+str(position[plotpos]))


plt.tight_layout()
plt.savefig(saving_directory+"individual_tracks.pdf")
print("done")




#%% fig 2 Plot with molt exit
def calculate_molt(condition,quantification_df):

    quantification_df=quantification_df[quantification_df["Condition"]==condition]
    avg_molting=quantification_df.groupby("Position").mean().mean()
    std_molting=quantification_df.groupby("Position").mean().std()


    return avg_molting, std_molting

plt.figure(2,figsize=[9,5])
sns.lineplot(x='Time',y='mean_curved',hue='Condition',data=CSV_quantification) #palette=color_graphs
plt.xlabel("Time after hatch (hours)")
plt.ylabel("GFP intensities (a.u.)")
plt.legend(labels=[conditions[0]+ " N="+str(number_of_worms[0]),
                    conditions[1]+ " N="+str(number_of_worms[1]),
                    #conditions[2]+ " N="+str(number_of_worms[2])
                    ],
                     bbox_to_anchor=(1.05, 1), loc='upper left')

for color,condition in zip(color_graphs,conditions):

    avg,std=calculate_molt(condition,CSV_quantification)
    if include_molt_exit==True:
        plt.axvspan(avg["TimeMolt1"]-std["TimeMolt1"],avg["TimeMolt1"]+std["TimeMolt1"],facecolor=color,alpha=0.2)
        plt.axvspan(avg["TimeMolt2"]-std["TimeMolt2"],avg["TimeMolt2"]+std["TimeMolt2"],facecolor=color,alpha=0.2)
        plt.axvspan(avg["TimeMolt3"]-std["TimeMolt3"],avg["TimeMolt3"]+std["TimeMolt3"],facecolor=color,alpha=0.2)
        plt.axvspan(avg["TimeMolt4"]-std["TimeMolt4"],avg["TimeMolt4"]+std["TimeMolt4"],facecolor=color,alpha=0.2)


plt.tight_layout()
plt.savefig(saving_directory+"all worms.pdf")



# %% plot with volume side by side, in 2D and 3D
plt.figure(3)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1, figsize=(5,6),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='mean_curved',hue='Condition',data=CSV_quantification,legend=True).set_title('properties curved 3D')
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',hue='Condition',data=CSV_quantification,legend=False)
axs[1].semilogy()
plt.savefig(saving_directory+"3D.pdf")


plt.figure(4)
sns.set_theme(style="darkgrid")
fig, axs = plt.subplots(2, 1, figsize=(5,6),sharex=True)
sns.lineplot(ax=axs[0],x='Time',y='mean_mip',hue='Condition',data=CSV_quantification,legend=True).set_title('properties MIP (2D)')
sns.lineplot(ax=axs[1],x='Time',y='area_mip',hue='Condition',data=CSV_quantification,legend=False)
axs[1].semilogy()
plt.savefig(saving_directory+"2D.pdf")

# %%
toplot=CSV_quantification[CSV_quantification["Condition"]=="LIN_28"]
molts_SWI = np.array([12.5, 21, 30, 42])

sns.set_theme(style="whitegrid")
fig, axs = plt.subplots(1, 2, figsize=(5,2),sharex=False)
sns.lineplot(ax=axs[0],x='Time',y='mean_curved',data=toplot,legend=True,ci='sd').set_title('Intensity of LIN-28::GFP')
axs[0].axvline(12.5,color='gray')
axs[0].axvline(21,color='gray')
axs[0].axvline(30,color='gray')
axs[0].set_ylabel("Intensity (a.u.)")
axs[0].set_xlabel("time (hours)")
sns.lineplot(ax=axs[1],x='Time',y='volume_curved',data=toplot,legend=False,ci="sd").set_title('Volume of worm over time')
axs[1].semilogy()
axs[1].axvline(12.5,color='gray')
axs[1].axvline(21,color='gray')
axs[1].axvline(30,color='gray')
axs[1].set_ylabel("volume (Pixel^3)")
axs[1].set_xlabel("time (hours)")
plt.tight_layout()
plt.savefig(saving_directory+"3D.pdf")



# %%