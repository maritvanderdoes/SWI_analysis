#%% import packages
#plotting SWI of experiment done on 12.2.21 seg1 and 2 vs seg 1 and 2scr
#import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#%% import files and merge together
CSV_quantification = pd.read_csv("kymograph_stats.csv")
CSV_annotation= pd.read_csv("goodworms_GRH1_detailed.csv")
min_image_taken=10 #time interval in minutes

label_condition1="scr2"
label_condition2="control"

def add_time(frame,start,end):
    if (frame>=start)and(frame<end):
        return (frame-start)*min_image_taken/60
    else:
        return None

#merge arrays for complete dataset
CSV_quantification.head()
CSV_annotation.head()

CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1
CSV_quantification['Time'] = CSV_quantification.apply(lambda row : add_time(row['Frame'],row['Hatch'],row['Escape']), axis = 1) #add overall time
CSV_quantification= CSV_quantification.dropna(subset=['Time'],axis='index') #selects only worms in between hatch and excape

goodworms=CSV_annotation[CSV_annotation['Quality']==1] #selects only worms where quality is 1
Nc1=np.sum(goodworms["Condition"]==label_condition1)
Nc2=np.sum(goodworms["Condition"]==label_condition2)

#%%filtering dataset for certain regions in time to plot

def filter_data(quantification_df,begin,end):
    """filter_data is a function that you give the quantification as input, and from where to where you want to plot it

    Args:
        quantification_df (pandas dataframe): Contains the 
        begin (str): Give the name of the column in the dataframe you have annotated the beginning
        end (str): Give the name of the column in the dataframe where you have annotated the end

    Returns:
        cond1, idx_cond1, cond2, idx_cond2 (np arrays): gives the names of the ]
    """
    #create arrays to plot for complete dataset
    quantification_df['Time_temp'] = quantification_df.apply(lambda row : add_time(row['Frame'],row[begin],row[end]), axis = 1) #add time column
    quantification_df= quantification_df.dropna(subset=['Time_temp'],axis='index') #selects only worms in between begin and and

    condition1_df=quantification_df[quantification_df["Condition"]==label_condition1].pivot(index="Time",columns="Position",values="Intensity_BGsub")
    condition2_df=quantification_df[quantification_df["Condition"]==label_condition2].pivot(index="Time",columns="Position",values="Intensity_BGsub")

    condition1=condition1_df.to_numpy()
    condition1_index=condition1_df.index.to_numpy()

    condition2=condition2_df.to_numpy()
    condition2_index=condition2_df.index.to_numpy()

    return condition1,condition1_index,condition2,condition2_index


[all_1, idx_all_1, all_2, idx_all_2]=filter_data(CSV_quantification,'Hatch','M4 exit') #all
[M1_1, idx_M1_1, M1_2, idx_M1_2]=filter_data(CSV_quantification,'M1 entry','M1 exit') #all
[M2_1, idx_M2_1, M2_2, idx_M2_2]=filter_data(CSV_quantification,'M2 entry','M2 exit') #all
[M3_1, idx_M3_1, M3_2, idx_M3_2]=filter_data(CSV_quantification,'M3 entry','M3 exit') #all
[M4_1, idx_M4_1, M4_2, idx_M4_2]=filter_data(CSV_quantification,'M4 entry','M4 exit') #all


if (Nc1==np.shape(M1_1)[1]==np.shape(M2_1)[1]==np.shape(M3_1)[1]==np.shape(M4_1)[1])==False:
    print("Data condition1 is not equal!")

if (Nc2==np.shape(M1_2)[1]==np.shape(M2_2)[1]==np.shape(M3_2)[1]==np.shape(M4_2)[1])==False:
    print("Data condition1 is not equal!")

#%% plot all worms together
#condition1
plt.plot(idx_all_1,all_1,color="blue",linewidth=0.05)
plt.plot(idx_M1_1,M1_1,color="red",linewidth=0.1) #M1
plt.plot(idx_M2_1,M2_1,color="red",linewidth=0.1) #M2
plt.plot(idx_M3_1,M3_1,color="red",linewidth=0.1) #M3
plt.plot(idx_M4_1,M4_1,color="red",linewidth=0.1) #M3

plt.plot(idx_all_1,np.mean(all_1,1),color="blue",linewidth=1) #average

#condition2
plt.plot(idx_all_2,all_2,color="orange",linewidth=0.05)
plt.plot(idx_M1_2,M1_2,color="black",linewidth=0.1) #M1
plt.plot(idx_M2_2,M2_2,color="black",linewidth=0.1) #M2
plt.plot(idx_M3_2,M3_2,color="black",linewidth=0.1) #M3
plt.plot(idx_M4_2,M4_2,color="black",linewidth=0.1) #M3

plt.plot(idx_all_2,np.mean(all_2,1),color="orange",linewidth=1) #average
# %% plot only certain positions
if Nc1>Nc2:
    fig, axs = plt.subplots(Nc1,2, figsize=(15,80),sharex=True)
else:
    fig, axs = plt.subplots(Nc2,2, figsize=(15,80),sharex=True)

for plotpos in range(Nc1):
    #condition1
    axs[plotpos,0].plot(idx_all_1,all_1[:,plotpos],color="blue",linewidth=1)
    axs[plotpos,0].plot(idx_M1_1,M1_1[:,plotpos],color="red",linewidth=1) #M1
    axs[plotpos,0].plot(idx_M2_1,M2_1[:,plotpos],color="red",linewidth=1) #M2
    axs[plotpos,0].plot(idx_M3_1,M3_1[:,plotpos],color="red",linewidth=1) #M3
    axs[plotpos,0].plot(idx_M4_1,M4_1[:,plotpos],color="red",linewidth=1) #M3
    axs[plotpos,0].plot(idx_all_1,np.mean(all_1,1),color="black",linestyle='dashed',linewidth=1) #average
    axs[plotpos,0].title.set_text(label_condition1+"_position "+str(goodworms["Position"].iloc[plotpos]))

for plotpos in range(Nc2):
    #condition2
    axs[plotpos,1].plot(idx_all_2,all_2[:,plotpos],color="orange",linewidth=1)
    axs[plotpos,1].plot(idx_M1_2,M1_2[:,plotpos],color="black",linewidth=1) #M1
    axs[plotpos,1].plot(idx_M2_2,M2_2[:,plotpos],color="black",linewidth=1) #M2
    axs[plotpos,1].plot(idx_M3_2,M3_2[:,plotpos],color="black",linewidth=1) #M3
    axs[plotpos,1].plot(idx_M4_2,M4_2[:,plotpos],color="black",linewidth=1) #M3
    axs[plotpos,1].plot(idx_all_2,np.mean(all_2,1),color="black",linestyle='dashed',linewidth=1) #average
    axs[plotpos,1].title.set_text(label_condition2+"_position "+str(goodworms["Position"].iloc[plotpos]))


plt.tight_layout()
plt.savefig("all worms")




#%% calculate timings molt
def timing_molt(hatch,molt):
    return (molt-hatch)*min_image_taken/60

def filter_data2(goodworms_df,quantification_df,begin,end):
    """

    Args:
        goodworms_df
        quantification_df (pandas dataframe)
        begin (str): Give the name of the column in the dataframe you have annotated the beginning
        end (str): Give the name of the column in the dataframe where you have annotated the end

    Returns:
        output_df (pandas dataframe): gives the names of the
    """


    #create arrays to plot for complete dataset
    goodworms["t_temp_entry"]=goodworms.apply(lambda row : timing_molt(row['Hatch'],row[begin]), axis = 1) #add overall time
    goodworms["t_temp_exit"]=goodworms.apply(lambda row : timing_molt(row['Hatch'],row[end]), axis = 1) #add overall time

    mean_timing=goodworms_df.groupby(['Condition']).mean()
    

    #calculating entry and exit time
    entry_c1=mean_timing["t_temp_entry"][label_condition1]
    exit_c1=mean_timing["t_temp_exit"][label_condition1]
    entry_c2=mean_timing["t_temp_entry"][label_condition2]
    exit_c2=mean_timing["t_temp_exit"][label_condition2]

    output_df=CSV_quantification[
    ((quantification_df.Condition==label_condition1) & (quantification_df.Time> entry_c1) & (quantification_df.Time< exit_c1))
    |
    ((quantification_df.Condition==label_condition2) & (quantification_df.Time> entry_c2) & (quantification_df.Time< exit_c2))
    ]


    return output_df


all=filter_data2(goodworms,CSV_quantification,"Hatch","M4 exit")

M1=filter_data2(goodworms,CSV_quantification,"M1 entry","M1 exit")
M2=filter_data2(goodworms,CSV_quantification,"M2 entry","M2 exit")
M3=filter_data2(goodworms,CSV_quantification,"M3 entry","M3 exit")
M4=filter_data2(goodworms,CSV_quantification,"M4 entry","M4 exit")






#%% 


plt.figure(1)
svm=sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=all,ci="sd")
svm=sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=M1, palette=['red','black'],legend=False,ci="sd")
svm=sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=M2, palette=['red','black'],legend=False,ci="sd")
svm=sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=M3, palette=['red','black'],legend=False,ci="sd")
svm=sns.lineplot(x='Time',y='Intensity_BGsub',hue='Condition',data=M4, palette=['red','black'],legend=False,ci="sd")

plt.legend(labels=[label_condition1+ " N="+str(Nc1),label_condition2+ " N="+str(Nc2)])


figure = svm.get_figure()    
figure.savefig('mean_and_std.png')


