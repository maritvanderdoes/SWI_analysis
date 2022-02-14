#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Loading tools
from utils import read_image_and_metadata
# Computing tools
from utils import adaptive_masking, make_mask_background, calculate_worm_properties
from utils import calculate_worm_properties_additional
from utils import create_skeleton, create_spline, head2tail_masking
from utils import straighten_image2D_dual_fast, arc_length
from skimage.measure import regionprops
from scipy.ndimage.measurements import label 

# Time out
from utils import timeout, saving_log
from utils.benchmarking import tic, toc

# Import additional libraries
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt
import os

# load parameters
from _parameters import data_format, debugging, sngpln
from _parameters import xdim, ydim, zdim   
from _parameters import outputpath, debugpath, data_format
from _parameters import time_out, steps

#%% define the main function
def main_function(list_mcherry,list_GFP):
    with timeout(time_out):
        # For debugging
        start0 = tic()

        # Determining the files to read
        files = (list_mcherry, list_GFP)

        # Reading the image and metadata
        print('MAIN: Reading image. File selected: '+files[0])
        (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

        # Creating debug log
        status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'File read')

        try:
            # Create mask
            img_binary, _  = adaptive_masking(img_mcherry, sorting=False, mm_th = 1.8, exp_size = 3)

            # dilate mask
            background_binary= make_mask_background(img_binary)

            # Maximum intensity projection
            # img_GFP_max=np.max(img_gfp*img_binary,0)
            # img_bin_max=np.max(img_binary,0)

            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'Masked Image', runtime = toc(start0, False))

            # Calculating curved worms properties 
            # (mean_mip, volume_mip)               = calculate_worm_properties(img_bin_max, img_GFP_max) 
            (mean_curved, volume_curved)         = calculate_worm_properties(img_binary, img_gfp) 
            (mean_background, volume_background) = calculate_worm_properties(background_binary, img_gfp) 
            current_res['mean_curved']     = mean_curved  #calculate volume
            current_res['volume_curved']   = volume_curved
            current_res['mean_background'] = mean_background
            # current_res['mean_mip']        = mean_mip
            # current_res['volume_mip']      = volume_mip

            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'Old props', runtime = toc(start0, False))

            # Computing all the features of the mask to determine where the worm is
            best_plane = np.sum(np.sum(img_binary,axis=2),axis = 1).argmax()
            if not(best_plane): best_plane = int(img_binary.shape[0]/2)

            props = calculate_worm_properties_additional(img_binary[best_plane,:,:], img_mcherry[best_plane,:,:]) 
            current_res['mask_mean'] = props[0]
            current_res['mask_area'] = props[1]
            current_res['mask_eccentricity'] = props[2]
            current_res['mask_centroid_x'] = props[3]
            current_res['mask_centroid_y'] = props[4]
            current_res['mask_perimeter'] = props[5]
            current_res['mask_total'] = props[6]
            current_res['mcherry_total'] = img_mcherry.sum()

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'New Props', runtime = toc(start0, False))
    
            # Saving logs
            stop_alg = toc(start0)
            status = saving_log(debugpath, 'COMPLETED', current_res['Frame'], current_res['Position'], 'Completed', runtime = stop_alg)


        except:
            print('MAIN: Some computations have not finished.')
            status = saving_log(debugpath, 'ERROR', current_res['Frame'], current_res['Position'], status +' (after error)', runtime = toc(start0, False))

            if debugging:
                # Creating the folder for the worm
                if sngpln:
                    foldername = outputpath+'/Single_Plane_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
                else:
                    foldername = outputpath+'/Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
                # Creating the individual debugging folder
                if not os.path.exists(foldername):
                    print('MAIN: Creating debugging folder.')
                    os.makedirs(foldername)

            # Saving Mask
            if debugging:
                imsave(foldername+'/A10_Mask'+'.tiff',np.float32(255*img_binary), check_contrast = False)
                imsave(foldername+'/A11_Masked_data'+'.tiff',np.float32(img_binary*img_gfp), check_contrast = False)
                imsave(foldername+'/A12_Data'+'.tiff',np.float32(img_mcherry), check_contrast = False)


        return current_res
