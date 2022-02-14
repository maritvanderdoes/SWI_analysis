
#%%
# To run the code type on cmd/terminal:
# python tasks_straightening_single/Running_Code.py

# Setting the directories
#from skimage.util.dtype import img_as_float
import _setup

# Loading tools
from utils import image_lists, read_image_and_metadata
# Computing tools
from utils import adaptive_masking, calculate_worm_properties, crop_image
# Benchmarking
from utils.benchmarking import tic, toc
# straightening
from utils import create_spline, head2tail_masking
from utils import straighten_image2D_dual_fast

from utils import arc_length
#from utils.deprecated import create_skeleton
from utils.core_utils import create_skeleton

# Time out
from utils import timeout, raise_timeout, saving_log

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Saving the mask
from skimage.io import imsave
from skimage.measure import regionprops, regionprops_table
from scipy.ndimage.measurements import label 
import skimage.morphology as skimorph

# load parameters
from _parameters import dirpath, outputpath, debugpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity, debugging, sngpln, straightenedname
from _parameters import time_out, steps
from _parameters import xdim, ydim, zdim    

#%% define the main function
def main_function(list_mcherry,list_GFP):
# with timeout(time_out):
    # For debugging
    start0 = tic()

    # Determining the files to read
    files = (list_mcherry, list_GFP)

    # Reading the image and metadata
    print('MAIN: Reading image. File selected: '+files[0])
    (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

    if np.mod(int(current_res['Frame']),10) == 0:
        # Creating debug log
        status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'File read')

        try:
            # Create mask
            # print('MAIN: Adaptive Masking')
            img_binary, _  = adaptive_masking(img_mcherry, mm_th = 2.5, th_sel = 0.3, exp_size = 3)

            # Saving Mask
            if debugging:
                imsave(straightenedname+'/A10_Mask'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',np.float32(255*img_binary), check_contrast = False)
                # Saving Masked Data
                imsave(straightenedname+'/A11_Masked_data'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',np.float32(img_binary*img_gfp), check_contrast = False)
            
            ccs, num_ccs = label(img_binary)
            properties=regionprops(ccs,img_gfp,['area','mean_intensity','min_intensity'])

            areas = [o.area for o in properties]
            selected = ccs == areas.index(max(areas))+1

            krn = skimorph.ball(3)
            surroundings = skimorph.binary_dilation(selected, krn)
            surroundings = surroundings*(np.invert(selected))

            current_res['mean_curved'] = np.mean(img_gfp[selected])
            current_res['median_curved'] = np.median(img_gfp[selected])
            current_res['min_curved'] = np.min(img_gfp[selected])
            current_res['max_curved'] = np.max(img_gfp[selected])
            current_res['mean_bg'] = np.mean(img_gfp[np.invert(selected)])
            current_res['median_bg'] = np.median(img_gfp[np.invert(selected)])
            current_res['min_bg'] = np.min(img_gfp[np.invert(selected)])
            current_res['max_bg'] = np.max(img_gfp[np.invert(selected)])
            current_res['mean_surr'] = np.mean(img_gfp[surroundings])
            current_res['median_surr'] = np.median(img_gfp[surroundings])
            current_res['min_surr'] = np.min(img_gfp[surroundings])
            current_res['max_surr'] = np.max(img_gfp[surroundings])

            if sngpln:
                print('Select single plane')
                area_cropped = np.sum(np.sum(selected,axis=2),axis = 1)
                best_plane = area_cropped.argmax()
                # Selecting the frame
                cropped_binary = selected[best_plane,:,:]
                cropped_surr = surroundings[best_plane,:,:]
                cropped_image = img_gfp[best_plane,:,:]
                # Reshaping
                cropped_binary = cropped_binary[None,:,:]
                cropped_surr = cropped_surr[None,:,:]
                cropped_image = cropped_image[None,:,:]
                

            # Calculating curved worms properties 
            current_res['mean_curved_sp'] = np.mean(cropped_image[cropped_binary])
            current_res['median_curved_sp'] = np.median(cropped_image[cropped_binary])
            current_res['min_curved_sp'] = np.min(cropped_image[cropped_binary])
            current_res['max_curved_sp'] = np.max(cropped_image[cropped_binary])
            current_res['mean_bg_sp'] = np.mean(cropped_image[np.invert(cropped_binary)])
            current_res['median_bg_sp'] = np.median(cropped_image[np.invert(cropped_binary)])
            current_res['min_bg_sp'] = np.min(cropped_image[np.invert(cropped_binary)])
            current_res['max_bg_sp'] = np.max(cropped_image[np.invert(cropped_binary)])
            current_res['mean_surr_sp'] = np.mean(cropped_image[cropped_surr])
            current_res['median_surr_sp'] = np.median(cropped_image[cropped_surr])
            current_res['min_surr_sp'] = np.min(cropped_image[cropped_surr])
            current_res['max_surr_sp'] = np.max(cropped_image[cropped_surr])
            # if mean_curved != mean_curved: return current_res
            
            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Cropped Image (Comp)', runtime = toc(start0, False))
            if steps == 'Cropping': return current_res


            # Saving logs
            stop_alg = toc(start0)
            status = saving_log(debugpath, 'COMPLETED', current_res['Frame'], current_res['Position'], 'Completed', runtime = stop_alg)


        except:
            print('MAIN: Some computations have not finished.')
            status = saving_log(debugpath, 'ERROR', current_res['Frame'], current_res['Position'], status +' (after error)', runtime = toc(start0, False))


    return current_res

# %%
