#%% importing
import pandas as pd
from marit_functions import image_lists_mcherry_GFP
from marit_functions import read_image
from marit_functions import img_thresholding
from marit_functions import select_zslides
from marit_functions import calculate_worm_properties
from marit_functions import get_meta_info_temp

from lucas_functions import generating_ball
from lucas_functions import adaptive_masking
from lucas_functions import dataset_comparison
from lucas_functions import masking_summary

#%% additional imports for coding
import matplotlib.pyplot as plt
import skimage.filters as skf
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import glob
import os
import re
from skimage import io
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops
import skimage

print(skimage.__version__)


#%% load parameters
dirpath = 'C:/Users/moraluca/Desktop/Lin28_test'
#dirpath = 'C:\Users\moraluca\Desktop\Lin28_test'
# directory = 'W:/tungsten/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210116/'

outputpath = 'C:/Users/moraluca/Desktop/Lin28_test/Output'
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#save retults
results = []
image=0

#%% to run

# list for all channels the stk files in folder
list_mcherry, list_GFP = image_lists_mcherry_GFP(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%%
#open mcherry and segment on signal
for i,(file1, file2) in enumerate(zip(list_mcherry, list_GFP)):
    print (i)
    if (i==image):
        #imageprint(file1)

        name = 'C:/Users/moraluca/Desktop/Lin28_test/Lin28_Project_210116_t1_s1_lHW2668p1_1_w2Lucas-sim-561-488.stk'
        img_mcherry2 = io.imread(name)
        print(img_mcherry2.shape)


        img_mcherry = read_image(file1)
        img_gfp = read_image(file2)

        # 0. Loading random data of same dimensions
        # Dimensions
        nz, nx, ny = (28, 120, 120)
        # Offsets
        sz, sx, sy = (0, .1, -.1)
        # Threshold level
        th_level = .5
        # Background level
        b_level = 10
        # Noise level
        n_level = 1e-1
        # Generating synthetic dataset
        mask_data, test_data, X_ball, Y_ball, Z_ball = generating_ball(nx, ny, nz, sx, sy, sz, th_level, n_level, b_level)

        # Parameters
        # Min/Max threshold
        mm_th = 1.05
        # Threshold_selection
        th_sel = .4

        # Running the function
        THPH, SORTED, ADPT, PRETH = adaptive_masking(test_data, mm_th, th_sel)

        # 4. Refining

        # Comparison of synthetic datasets
        dataset_comparison(mask_data, test_data, THPH)

        # Presenting outputs
        masking_summary(PRETH, ADPT, mm_th)

        # 4. Additional thresholding
        # krn = strel('disk',5);
        # for pl = 1:28
        #     THPX(:,:,pl) = imdilate(imerode(THPX(:,:,pl),krn,'same'),krn,'same');
        #     THPX(:,:,pl) = imfill(THPX(:,:,pl),'holes');
        # end

        # #preprocessing and thresholding to make image binary
        # img_binary = img_thresholding(img_mcherry)

        # #select slides based on mean signal
        # img_signal, img_binary = select_zslides(img_gfp, img_binary)

        # #calculates properties of the segmented worm
        # binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
        # img_binary, img_signal)

        # #create overlay of binary image with GFP image
        # img_overlay = img_signal * img_binary

        # #add properties in current results
        # current_res = get_meta_info_temp(file2)  #get metadata
        # current_res['volume'] = area  #calculate volume
        # current_res['mean_intensity'] = mean_intensity
        # current_res['min_intensity'] = min_intensity
        # current_res[final_intensity'] = mean_intensity - min_intensity  #calculate intensity

        # #save in resultarray
        # results.append(current_res)

        break
 # %%
