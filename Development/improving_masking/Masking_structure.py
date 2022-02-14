#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
# To run the code type on cmd/terminal:
# python tasks_clean/Running_Code.py

import _setup
# Loading tools
from utils import image_lists, single_image_lists, read_image_and_metadata
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

# load parameters
from _parameters import dirpath, outputpath, debugpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity, debugging, sngpln
from _parameters import time_out, steps
from _parameters import xdim, ydim, zdim    

# % Loading new data
# Initialisation
results = []
current_res_df = []
image = 0

# Selecting sample
s_sel = 1
t_sel = 100

# list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel_GFP, s_sel = s_sel, t_sel = t_sel, channel3 = None, data_format = data_format)

#%%

# For debugging
start0 = tic()

files = (list_mcherry[0], list_GFP[0])

# Reading the image and metadata
print('MAIN: Reading image. File selected: '+files[0])
(img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

mm_th = 4;
th_sel = 0.1;
krn_size = 3;
exp_size = 3;


# Create mask
# print('MAIN: Adaptive Masking')
img_binary, _  = adaptive_masking(img_mcherry,  mm_th = mm_th, th_sel = th_sel, krn_size = krn_size,
    krn_type = 'Disk', exp_size = exp_size, fill_holes = True, z_threshold = 0.7, 
    sorting = False, verbose = False)

# Direct calculations
current_res['total_direct'] = np.sum(img_binary*img_gfp)  #calculate volume
current_res['volume_direct'] = np.sum(img_binary)/(xdim*ydim*zdim)

# Calculating curved worms properties 
(mean_curved, volume_curved) = calculate_worm_properties(img_binary, img_gfp)
current_res['mean_curved'] = mean_curved  #calculate volume
current_res['volume_curved'] = volume_curved/(xdim*ydim*zdim)# if mean_curved != mean_curved: return current_res

plt.imshow(np.sum(img_binary,axis = 0))

area_cropped = np.sum(np.sum(img_binary,axis=2),axis = 1)
best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
# Selecting the frame
cropped_binary = img_binary[best_plane,:,:]
cropped_image = img_gfp[best_plane,:,:]
        
plt.figure()
plt.imshow(cropped_binary)


#%%

# For debugging
start0 = tic()

files = (list_mcherry[0], list_GFP[0])

# Reading the image and metadata
print('MAIN: Reading image. Fi/le selected: '+files[0])
(img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

# Parameters to change
mm_th_vec = np.linspace(1.5,5,11)
th_sel_vec = np.linspace(0.1,0.9,9)
krn_size_vec = np.linspace(3,3,1)
exp_size_vec = np.linspace(3,3,1)

calculations = np.zeros([11,9,5,5,6])

for (mm_th_k,mm_th) in enumerate(mm_th_vec):
    print('mm_th_k='+str(mm_th))
    for (th_sel_k,th_sel) in enumerate(th_sel_vec):
        print('th_sel_k='+str(th_sel))
        for (krn_size_k,krn_size) in enumerate(krn_size_vec):    
            for (exp_size_k,exp_size) in enumerate(exp_size_vec):
                # Create mask
                # print('MAIN: Adaptive Masking')
                img_binary, _  = adaptive_masking(img_mcherry,  mm_th = mm_th, th_sel = th_sel, krn_size = krn_size,
                    krn_type = 'Disk', exp_size = exp_size, fill_holes = True, z_threshold = 0.7, 
                    sorting = False, verbose = False)

                # Direct calculations
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,0] = np.sum(img_binary*img_gfp)  #calculate volume
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,1] = np.sum(img_binary)/(xdim*ydim*zdim)

                # Calculating curved worms properties 
                (mean_curved, volume_curved) = calculate_worm_properties(img_binary, img_gfp)
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,2] = mean_curved  #calculate volume
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,3] = volume_curved/(xdim*ydim*zdim)
                # if mean_curved != mean_curved: return current_res

                area_cropped = np.sum(np.sum(img_binary,axis=2),axis = 1)
                best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
                # Selecting the frame
                cropped_binary = img_binary[best_plane,:,:]
                cropped_image = img_gfp[best_plane,:,:]

                (mean_curved, volume_curved) = calculate_worm_properties(cropped_binary, cropped_image)
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,4] = mean_curved  #calculate volume
                calculations[mm_th_k,th_sel_k,krn_size_k,exp_size_k,5] = volume_curved/(xdim*ydim*zdim)
                        
                # plt.imshow(cropped_binary)


_ = toc(start0)
# %%
import seaborn as sns
sns.heatmap(np.log10(calculations[:,:,1,1,1]),xticklabels=th_sel_vec.round(2),
            yticklabels=mm_th_vec.round(2))

plt.figure()
sns.heatmap(np.log10(calculations[:,:,1,1,2]),xticklabels=th_sel_vec.round(2),
            yticklabels=mm_th_vec.round(2))

plt.figure()
plt.plot(mm_th_vec.round(2),np.log10(calculations[:,3,1,1,1]))

plt.figure()
plt.plot(th_sel_vec.round(2),(calculations[5,:,1,1,1]))
# %%
