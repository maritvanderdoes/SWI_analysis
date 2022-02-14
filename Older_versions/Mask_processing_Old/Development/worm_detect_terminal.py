# %%
# FUll code
from utils import image_lists, single_image_lists, read_image_and_metadata
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import imageio

# Creating some functions

#%%

# load parameters
# dirpath         = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210305_Compressed'
dirpath         = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210305_Compressed'
dirsel = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_HetPathQuant_Results'

channel_GFP     = 'w1Lucas-sim-488-561'
channel_mcherry = 'w2Lucas-sim-561-488'

filename = '/cleanup_Results_LJ_Lin28_Project_210305_1.csv'

motion_level = 8
motion_window = 5
ec_th = 0.9
vol_change_th = 8

# Opening the file
print('Opening the file')
preanalysis_csv = pd.read_csv(dirsel+filename)
preanalysis_csv = preanalysis_csv.sort_values(['Position','Frame'])

# Checking NAs
# print(preanalysis_csv.isnull().sum())

#%%
print('Checking the list')
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)

#%% Iterating files
print('Sorting the list')
colset = ['Frame', 'Position', 'Condition', 'Name']
list_df = pd.DataFrame(columns = colset)

for k_f in np.arange(0,len(list_mcherry)):
    if np.mod(k_f,200)==0:
        print(k_f)
    # (directory, filename) = list_mcherry[k_f].split('\\')
    filename = list_mcherry[k_f].split('/')[-1]
    (basename, rest) = filename.split('_t')
    (frame_set, rest) = rest.split('_s')
    (position_set, rest) = rest.split('_l')
    (condition_set, rest) = rest.split('p')

    list_df=list_df.append(dict(zip(colset,[int(frame_set), int(position_set), condition_set, list_mcherry[k_f]])), ignore_index=True)


#%%

print('Iterating')
# dirsel = 'D:/GitHub/SWI_analysis/output/'

# iterating_vec = np.arange(4,preanalysis_csv['Position'].max()+1)

iterating_vec = preanalysis_csv.groupby('Position').median()['mask_eccentricity'].sort_values(ascending=False).index

for position_sel in iterating_vec[0:]:
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
    plt.savefig(dirsel+'/worm_detect/'+basename+'_s'+str(position_sel)+'_Feature_analysis'+'.pdf')
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


    fig.savefig(dirsel+'/worm_detect/'+basename+'_s'+str(position_sel)+'_Sample_analysis'+'.pdf', dpi = 200)
    plt.close()

# %%
