#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
# To run the code type on cmd/terminal:
# python tasks_clean/Running_Code.py

# Setting the directories
import _setup

# Loading custom libraries
from utils import image_lists, single_image_lists
from utils.benchmarking import tic, toc

# Import additional libraries
import numpy as np
import pandas as pd

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format

# Import main functions
from _main_function import main_function

# %% MAIN CODE
# Initialisation
results = []
current_res_df = []
image = 0

# Selecting sample
s_sel = 35
t_sel = 114

# list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel_GFP, s_sel = s_sel, t_sel = None, channel2 = None, channel3 = None, data_format = data_format)

#%% Running for a single sample
current_res = main_function(list_mcherry[0],list_GFP[0])

area_cropped = np.sum(np.sum(cropped_binary,axis=2),axis = 1)
            best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
            # Selecting the frame
            cropped_binary = cropped_binary[best_plane,:,:]
            cropped_image = cropped_image[best_plane,:,:]
            # Reshaping
            cropped_binary = cropped_binary[None,:,:]
            cropped_image = cropped_image[None,:,:]
            

# Saving results
df = pd.DataFrame(current_res)
df.to_csv(outputpath+'/Results_t'+str(t_sel)+'_s'+str(s_sel)+'.csv', index = False)
