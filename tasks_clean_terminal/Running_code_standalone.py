#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_final/Running_Code.py

# Loading custom libraries
from utils import image_lists, file_validation, read_image_and_metadata
from utils import adaptive_masking, make_mask_background, calculate_worm_properties
from utils import timeout, saving_log
from utils.benchmarking import tic, toc

# Import additional libraries
import multiprocessing
import numpy as np
import pandas as pd

# load parameters
dirpath         = '/tungstenfs/scratch/ggrossha/Marit/HetP_Quant/20201224_HBL1'
outputpath      = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_HetPathQuant_Results'
output_filename = 'Results_20201224_HBL1_test.csv'
debugpath       = outputpath + '/Debug'
channel_GFP     = 'w1Lucas-sim-488-561'
channel_mcherry = 'w2Lucas-sim-561-488'

# Additional parameters
n_workers = 12
data_format = 'ts'
time_out = 200 # seconds

# Scaling factors
xdim = 0.533
ydim = 0.533
zdim = 2
counter = 0

#%% MAIN FUNCTION
def main_function(imglist_mcherry,imglist_GFP):
    with timeout(time_out):
        # For debugging
        start0 = tic()

        # Determining the files to read
        files = (imglist_mcherry, imglist_GFP)

        # Reading the image and metadata
        print('MAIN: Reading image. File selected: '+files[0])
        (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

        # Creating debug log
        status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'File read')

        try:
            # Create mask
            img_binary, _  = adaptive_masking(img_mcherry, sorting=False, mm_th = 2, exp_size = 3)

            # Computing background
            background_binary= make_mask_background(img_binary)

            # max intensity projection
            img_GFP_max=np.max(img_gfp*img_binary,0)
            img_bin_max=np.max(img_binary,0)

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'Masked Image', runtime = toc(start0, False))
 
            #calculate properties
            (mean_mip, volume_mip)               = calculate_worm_properties(img_bin_max, img_GFP_max) 
            (mean_curved, volume_curved)         = calculate_worm_properties(img_binary, img_gfp) 
            (mean_background, volume_background) = calculate_worm_properties(background_binary, img_gfp) 
            
            # Storing the properties in current results for curved worm
            current_res['mean_curved']     = mean_curved  #calculate volume
            current_res['volume_curved']   = volume_curved /(xdim*ydim*zdim)
            current_res['mean_background'] = mean_background
            current_res['mean_mip']        = mean_mip
            current_res['volume_mip']      = volume_mip /(xdim*ydim*zdim)
   
            # Saving logs
            stop_alg = toc(start0)
            status = saving_log(debugpath, 'COMPLETED', current_res['Frame'], current_res['Position'], 'Completed', runtime = stop_alg)

            # global counter
            # global list_mcherry
            # counter = counter + 1
            # print(str(counter) + "out of "+ str(len(list_mcherry)/12) + "is done")

        except:
            print('MAIN: Some computations have not finished.')
            status = saving_log(debugpath, 'ERROR', current_res['Frame'], current_res['Position'], status +' (after error)', runtime = toc(start0, False))

        return current_res

# %% RUNNING CODE
# List for all channels the image files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)

# File validation
file_validation(dirpath, outputpath, output_filename, list_mcherry, list_GFP)

# Actual run
if __name__=='__main__':
    p = multiprocessing.Pool(n_workers)
    results = p.starmap(main_function,zip(list_mcherry,list_GFP))
    p.close()

    # Saving files
    df = pd.DataFrame(results)
    df.to_csv(outputpath+'/'+output_filename, index = False)

