#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
To ensure the code is as compatible between both users, Marit and 
Lucas, this function takes care of setting that for both of them.

To ensure clarity, the code also set ups the different parameters
that do not require to be changed constantly.
"""

import os

#%% Finding directories
# Default working directories
Lucas_path = 'C:/Users/moraluca/Desktop/Lin28_test'
Lucas_data_format = 'ts'
Marit_path = '/Users/Marit/Documents/work/test_images/HBL1gfp_worm6'
Marit_data_format = 'ts'

# Setting the directoy
if os.path.exists(Lucas_path):
    dirpath = Lucas_path
    data_format = Lucas_data_format
    print('Directory set for Lucas, as', end=" ")
elif os.path.exists(Marit_path):
    dirpath = Marit_path
    data_format = Marit_data_format
    print('Directory set for Marit, as', end=" ")
else:
    print('No available paths have been found. No data_format has selected')

if os.path.exists(dirpath):
    print(dirpath)
    print('Selected data_format: '+data_format)
   
# Create output folder
outputpath = os.path.abspath(os.path.join(dirpath,'Output'))
if os.path.exists(outputpath):
    print('Output folder already existing')
else:
    print('Folder output has been created')
    os.makedirs(os.path.abspath(os.path.join(dirpath,'Output')))

#%% Defining channels
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#%% Numerical Parameters
# Min/Max threshold
#mm_th = 1.8
# Threshold_selection
#th_sel = .4
# Kernel
#krn_size = 5
#krn_type = 'Disk'
# Expansion kernel
#exp_size = 5

# Further refinement of the mask in z-planes
#z_threshold = 0.2