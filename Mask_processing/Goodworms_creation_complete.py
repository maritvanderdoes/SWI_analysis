#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_final/Running_Code.py

# Import libraries
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import imageio

from skimage.measure import regionprops
from scipy.ndimage.measurements import label 

# Select the segmentation mode (mcherry or cnn)
segment_mode    = 'mcherry'

# load parameters
# dirpath         = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/woelmich/20211126_nhr-23GFP_nhr25i/20211126_segmentation'
segmentpath     = '/tungstenfs/scratch/ggrossha/woelmich/20211203_blmp10mu_nhr25GFP/20211203_segmentation'
filespath       = '/tungstenfs/scratch/ggrossha/woelmich/20211203_blmp10mu_nhr25GFP/20211203_compressed'
dirsel          = '/tungstenfs/scratch/ggrossha/woelmich/20211203_blmp10mu_nhr25GFP/20211203_summary_2'
filename        = '20211203_Worm_segmentation_results_Michi.csv'
channel_BF      = 'w2smita-Brightfield-GFP'
channel_GFP     = 'w2chiara-488-561sim'
channel_mcherry = 'w1chiara-561-488sim'
data_format     = 'st'
file_format     = 'tif'

motion_level = 4 #5
motion_window = 5
ec_th = 0.9
vol_change_th = 6 #8
entry_th = 2

(n_plots_th, quality_th)     = (8,0.3)
(n_workers, dpi_res, num_th) = (12,200,64)

# print(dirsel)

# Checking for directory existence
if os.path.exists(dirsel):
    print('Output folder already existing.', end = ' ')

    if os.path.exists(dirsel+'/'+filename):
        print('Statistics file already existing. Running worm detection.')
        run_mask_processing = False
    else:
        print('Statistics file does not exists. Running full code.')
        run_mask_processing = True

else:
    print('Folder output has been created. Running full code.')
    run_mask_processing = True
    os.makedirs(dirsel)