#%% import packages
#plotting SWI of experiment done on 12.2.21 seg1 and 2 vs seg 1 and 2scr
#import packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#%% import csv files
#import csv
CSV_quantification = pd.read_csv("kymograph_stats.csv")
CSV_annotation= pd.read_csv("12221_valid_escape.csv")
CSV_quantification.head()
CSV_annotation.head()
min_image_taken=10 #time interval in minutes

#%% get molting points from quantification file


# %% filtering dataset for quality and Time
def select_timepoints(frame,start,end):
    if (frame>=start)and(frame<end):
        return (frame-start)*min_image_taken/60
    else:
        return None

def molt_annotation(frame,molt):
    return (frame-molt)*min_image_taken/60


CSV_quantification=CSV_quantification.merge(CSV_annotation,on=['Position'],how='inner') #merges two datasets
CSV_quantification['TimeHatch'] = CSV_quantification.apply(lambda row : select_timepoints(row['Frame'],row['hatch'],row['escape']), axis = 1) #add time column
CSV_quantification= CSV_quantification.dropna(subset=['TimeHatch'],axis='index') #selects only worms in between hatch and excape
CSV_quantification['TimeMolt1'] = CSV_quantification.apply(lambda row : molt_annotation(row['Frame'],row['M1_exit']), axis = 1) #set timing based on molt1
CSV_quantification['Intermolt1']=CSV_quantification.apply(lambda row:calculate_molt(row['hatch'],row['M1_exit']),axis=1) #add timing intermolt
Molting1=CSV_quantification[CSV_quantification.TimeMolt1==0] # make a dataframe with only the points of the intermolt.
#CSV_quantification=CSV_quantification[CSV_quantification['Quality']==1] #selects only worms where quality is 1


#%% calculating where molt is
def calculate_molt(hatch,m1):
        return (m1-hatch)*min_image_taken/60

CSV_annotation['Intermolt1']=CSV_annotation.apply(lambda row:calculate_molt(row['hatch'],row['M1_exit']),axis=1)
condition1= CSV_annotation[CSV_annotation['condition']=='seg1and2']
condition2= CSV_annotation[CSV_annotation['condition']=='seg1and2scr']
condition1_mean=condition1['Intermolt1'].mean()
condition1_std=condition1['Intermolt1'].std()
condition2_mean=condition2['Intermolt1'].mean()
condition2_std=condition2['Intermolt1'].std()

#%% filtering

filtered_CSV_quantification=CSV_quantification[
    ((CSV_quantification.TimeHatch>(condition1_mean-condition1_std))
    & (CSV_quantification.TimeHatch<(condition1_mean+condition1_std))
    & (CSV_quantification.condition=='seg1and2'))
    | ((CSV_quantification.TimeHatch>(condition2_mean-condition2_std))
    & (CSV_quantification.TimeHatch<(condition2_mean+condition2_std))
    & (CSV_quantification.condition=='seg1and2scr'))
    ]


#%% # plots whole dataset
plt.figure(1)
sns.lineplot(x='TimeMolt1',y='Intensity_BGsub',hue='condition',data=CSV_quantification).set_title('alligned on first molt')

plt.figure(2)
sns.lineplot(x='TimeHatch',y='Intensity_BGsub',hue='condition',data=CSV_quantification).set_title('alligned on hatch')
sns.lineplot(x='TimeHatch',y='Intensity_BGsub',hue='condition', data=filtered_CSV_quantification,palette=['black','red'])

plt.figure(3)
sns.lineplot(x='TimeHatch',y='Intensity_BGsub',hue='condition',data=CSV_quantification).set_title('alligned on hatch')
plt.axvspan(condition1_mean-condition1_std, condition1_mean+condition1_std,facecolor='blue',alpha=0.2)
plt.axvspan(condition2_mean-condition2_std, condition2_mean+condition2_std,facecolor='orange',alpha=0.2)
plt.axvline(condition1_mean,color='blue')
plt.axvline(condition2_mean,color='orange')


