# %%
# FUll code
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import imageio
import multiprocessing

import os
import glob

# load parameters
filespath       = '/tungstenfs/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_compression'
dirsel          = '/tungstenfs/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_summary'
filename        = 'Worm_segmentation_results_Michi.csv'
channel_BF      = 'w2smita-Brightfield-GFP'
data_format     = 'st'

motion_level = 5
motion_window = 5
ec_th = 0.9
vol_change_th = 8
entry_th = 2

n_plots_th = 6
quality_th = 0.3
n_workers = 12
dpi_res = 200

# RUNNING CODE
def run_main():
    # Opening the file
    print('Opening the file')
    preanalysis_csv = pd.read_csv(dirsel+'/'+filename)
    preanalysis_csv = preanalysis_csv.sort_values(['Position','Frame'])

    # Create the mask list
    print('File indexing')
    list_BF=sorted(glob.glob(os.path.join(filespath, "*"+channel_BF+"*.tif")))

    #%% Iterating files
    print('Sorting the list')
    p = multiprocessing.Pool(n_workers)
    results = p.map(folder_documentation,list_BF)
    list_df = pd.DataFrame(results)

    #%% Computing the files
    preanalysis_csv_list = [preanalysis_csv[preanalysis_csv['Position']==k] for k in np.unique(list_df['Position'])]
    list_df_vec          = [list_df[list_df['Position']==k] for k in np.unique(list_df['Position'])]

    # Creating the graphs
    print('Creating the graphs')
    goodworms_temp = p.starmap(image_generation,zip(preanalysis_csv_list,list_df_vec))
    p.close()
    
    # Saving files
    df = pd.DataFrame(goodworms_temp)
    df.to_csv(dirsel+'/goodworms_temp.csv', index = False)

    print('Finished')


#%% FUNCTIONS
def folder_documentation(list_file):
    filename = list_file.split('/')[-1]

    if data_format == 'st':
        (basename, rest) = filename.split('_s')
        (position_set, rest) = rest.split('_t')
        (frame_set, rest) = rest.split('.tif')
    elif data_format == 'ts':
        (basename, rest) = filename.split('_t')
        (frame_set, rest) = rest.split('_s')
        (position_set, rest) = rest.split('.tif')

    # Creating dictionary
    colset = ['Frame', 'Position', 'Basename', 'Name']
    filename_info = dict(zip(colset,[int(frame_set), int(position_set),  basename, list_file]))

    return filename_info

def image_generation(selection_csv,subset):
    names = ['hatching_estimates', 'escape_estimates', 'worm_entry_estimates']  
    position = pd.unique(subset['Position'])[0]
    basename = pd.unique(subset['Basename'])[0]
    print('Position selected: '+str(position))

    # New metrics
    selection_csv['displacement'] = \
        np.sqrt((selection_csv['mask_centroid_x'].diff())**2+\
                (selection_csv['mask_centroid_y'].diff())**2)
    selection_csv['volume_diff'] = selection_csv['mask_area'].diff()

    # Computing the metrics
    new_hatches_m = hatching_determination_motion(selection_csv, motion_level)
    new_hatches_e = hatching_determination_eccentricity(selection_csv, threshold = 2)
    hatching_estimates = np.unique(np.concatenate([new_hatches_m,new_hatches_e]))
    escape_estimates = escape_determination(selection_csv, vol_change_th)
    worm_entry_estimates = worm_entry_determination(selection_csv, threshold= entry_th)

    # Check images
    n_plots = max(min(hatching_estimates.shape[0],n_plots_th),
                min(escape_estimates.shape[0],n_plots_th),
                min(worm_entry_estimates.shape[0],n_plots_th),)

    # Plotting hatching estimates
    fig = plt.figure(figsize = (18,n_plots*5))
    gs = fig.add_gridspec(n_plots+2,9)

    axs_growth = fig.add_subplot(gs[0:2, 1:-1])
    axs_growth.semilogy(selection_csv['Frame'],selection_csv['mask_total'],'k:')
    axs_growth.semilogy(selection_csv['Frame'],selection_csv['mask_area'],'k')
    axs_growth.set_xlabel('Time')
    # Detecting where the worm hatches
    for k_est, estimates in enumerate([hatching_estimates, escape_estimates, worm_entry_estimates]):
        color = np.zeros(3)
        color[k_est] = 1
        [axs_growth.axvline(estimates[k],color = color) for k in np.arange(0,min(estimates.shape[0],n_plots_th))]
        [axs_growth.axvline(estimates[k], linestyle = ':',color = color + (1-color)*.3) for k in np.arange(min(estimates.shape[0],n_plots_th),estimates.shape[0])]

    for k_est, estimates in enumerate([hatching_estimates, escape_estimates, worm_entry_estimates]):
        
        for prob_set in range(0,min(estimates.shape[0],n_plots_th)):
            t_sel = estimates[prob_set]

            for k_t, t_set in enumerate([t_sel-1,t_sel,t_sel+1]):
                if (t_set>0)&(t_set<=subset['Frame'].max()):
                    filename_it = subset.loc[subset['Frame']==t_set,'Name'].values[0]
                    img_mcherry = np.array(imageio.mimread(filename_it))

                    axs_image = fig.add_subplot(gs[prob_set+2, k_t + 3*k_est])
                    axs_image.imshow(img_mcherry.max(axis = 0), cmap='gray')
                    axs_image.set_title('f: '+str(t_set))
                    axs_image.axis('off')

    line = plt.Line2D((.38,.38),(n_plots/(n_plots+2)*.9,.12), color="k", linewidth=3)
    fig.add_artist(line)
    line = plt.Line2D((.645,.645),(n_plots/(n_plots+2)*.9,.12), color="k", linewidth=3)
    fig.add_artist(line)


    axs_growth.set_title('Sample: '+str(position)+\
        ', NAs: '+str(selection_csv.isnull().sum()['mask_mean'])+ ', mean_ecc: '+str(np.round(selection_csv['mask_eccentricity'].mean(),2)),
        fontsize = 20)

    figure_to_save = dirsel+'/'+basename+'_s'+str(position)+'_Full_analysis'+'.pdf'
    print('Saving position: '+str(position))
    fig.savefig(figure_to_save, dpi = dpi_res)


    if hatching_estimates.shape[0] == 1: hatch_sel = hatching_estimates[0]
    else: hatch_sel = hatching_estimates[1]
    if escape_estimates.shape[0] == 1: escape_sel = escape_estimates[0]
    else: escape_sel = escape_estimates[1]
    
    # hatch_sel = hatching_estimates[np.min(1,len(hatching_estimates)-1)]
    # escape_sel =  escape_estimates[np.min(1,len(escape_estimates)-1)]

    quality = min(selection_csv['Frame'].max()-selection_csv.isnull().sum()['mask_mean'],(escape_sel-hatch_sel))>selection_csv['Frame'].max()*quality_th

    # Creating dataframe
    colset = ['Position', 'Quality', 'Hatch', 'Escape', 'Condition']
    current_res = dict(zip(colset,[int(position), 
                int(quality),
                int(hatch_sel), int(escape_sel), ' ']))

    # plt.close()
    return current_res


#%% SUPPORTING FUNCTIONS
def escape_determination(df, threshold = 8, tps = 'Frame', volume_diff = 'volume_diff', area = 'mask_area', add_end = True):
    exit_likelihood = (-(df[volume_diff]/df[area])>threshold) | \
                      (  df[volume_diff].isnull().astype(float).diff()>0)
    
    escape_estimates = df.loc[exit_likelihood,tps].values

    if add_end or (escape_estimates.size==0):
        escape_estimates = np.append(np.array([df[tps].max()-1]),escape_estimates)

    return escape_estimates

def hatching_determination_motion(df, threshold = 8, window = 3, eccentricity_th = .9, tps = 'Frame', displacement = 'displacement', eccentricity = 'mask_eccentricity', add_start = True):
    # On motion
    df["egg_likelihood_motion"] = ((df[displacement]<threshold).astype(float).rolling(window,center = True).mean().values>.5).astype(float)
    hatching_estimates_motion = df["egg_likelihood_motion"].diff()<0

    hatching_estimates = df.loc[hatching_estimates_motion,tps].values

    if add_start or (df.loc[df[tps]==1,eccentricity].values > eccentricity_th):
        hatching_estimates = np.append([2],hatching_estimates)

    return hatching_estimates

def hatching_determination_eccentricity(df, threshold = 2, window = 1, eccentricity_th = .9, tps = 'Frame', eccentricity = 'mask_eccentricity', add_start = True):
    fun = np.log(df[eccentricity]).rolling(window, center = True).mean().diff()
    hatching_estimates = df.loc[fun<-threshold*fun.std(),tps].values

    if add_start or (df.loc[df[tps]==1,eccentricity].values > eccentricity_th):
        hatching_estimates = np.append([2],hatching_estimates)

    return hatching_estimates

def worm_entry_determination(df, threshold = 1, tps = 'Frame', mask_total = 'mask_total', mask_area = 'mask_area'):
    fun = (df[mask_total]-df[mask_area])
    worm_entry_likelihood = (fun>threshold*fun.mean()).astype(float).diff()>0
    worm_entry_estimates = df.loc[worm_entry_likelihood,tps].values

    return worm_entry_estimates

# %%
# Actual run
if __name__=='__main__':
    run_main()
