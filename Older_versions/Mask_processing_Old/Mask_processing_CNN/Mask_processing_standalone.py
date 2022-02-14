#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_final/Running_Code.py

# Import additional libraries
import multiprocessing
import numpy as np
import pandas as pd
import os
import glob
import imageio

from skimage.measure import regionprops
from scipy.ndimage.measurements import label 

# load parameters
# dirpath         = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_segmentation'
dirpath         = '/tungstenfs/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_segmentation'
outputpath      = '/tungstenfs/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_summary'
output_filename = 'Worm_segmentation_results_Michi.csv'
data_format     = 'st'

if os.path.exists(outputpath):
    print('Output folder already existing')
else:
    print('Folder output has been created')
    os.makedirs(outputpath)

# Additional parameters
n_workers = 12
num_th = 64

# RUNNING CODE
def run_main():
    # Create the mask list
    list_mask=sorted(glob.glob(os.path.join(dirpath, "*.tif")))

    p = multiprocessing.Pool(n_workers)
    results = p.map(main_function,list_mask)
    p.close()

    # Saving files
    df = pd.DataFrame(results)
    df.to_csv(outputpath+'/'+output_filename, index = False)

#%% MAIN FUNCTION
def main_function(file_sel):
    # Reading the image and metadata        
    # print('MAIN: Reading image. File selected: '+file_sel)
    # (directory, filename) = file_sel.split('\\')
    filename = file_sel.split('/')[-1]
    if data_format == 'st':
        (basename, rest) = filename.split('_s')
        (position_set, rest) = rest.split('_t')
        (frame_set, rest) = rest.split('.tif')
    elif data_format == 'ts':
        (basename, rest) = filename.split('_t')
        (frame_set, rest) = rest.split('_s')
        (position_set, rest) = rest.split('.tif')

    # Creating dataframe
    colset = ['Frame', 'Position']
    current_res = dict(zip(colset,[int(frame_set), int(position_set)]))

    # Read image
    img_binary = np.array(imageio.imread(file_sel))>num_th
    # Calculating properties of the mask
    ccs, num_ccs = label(img_binary) #set labels in binary image
    
    if num_ccs>0:
        properties=regionprops(ccs,img_binary,['area',"mean_intensity",'centroid','eccentricity','perimeter']) #calculates the properties of the different areas
        best_region = max(properties, key=lambda region: region.area) #selects the biggest region
    
    # Saving
        current_res['mask_mean']         = best_region.mean_intensity
        current_res['mask_area']         = best_region.area
        current_res['mask_eccentricity'] = best_region.eccentricity
        current_res['mask_centroid_x']   = best_region.centroid[0]
        current_res['mask_centroid_y']   = best_region.centroid[1]
        current_res['mask_perimeter']    = best_region.perimeter
        current_res['mask_total']        = np.sum(ccs) 
        current_res['mask_regions']      = num_ccs

    return current_res

# %%
# Actual run
if __name__=='__main__':
    run_main()