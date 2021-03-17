#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
To load paralell folders it is required to append the folder utils to
the path. This is done by this function.

Normally the path is removed after reloading the kernel, so it should
not be necessary to manually remove it. In any case, I have put the 
if statement to ensure there is no additional instances of the path.

To manually remove the path, check below.
"""

import os

#%% Finding directories
# Default working directories
Lucas_path = 'C:/Users/moraluca/Desktop/Lin28_test'
Marit_path = '/Users/Marit/Documents/work/test_images/HBL1gfp_worm6'

# Setting the directoy
if os.path.exists(Lucas_path):
    dirpath = Lucas_path
    print('Directory set for Lucas, as', end=" ")
elif os.path.exists(Marit_path):
    dirpath = Marit_path
    print('Directory set for Marit, as', end=" ")
else:
    print('No available paths have been found')

if os.path.exists(dirpath):
    print(dirpath)
   
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

# Numerical Parameters
# Min/Max threshold
mm_th = 1.8
# Threshold_selection
th_sel = .5
# Kernel
krn_size = 5
krn_type = 'Disk'