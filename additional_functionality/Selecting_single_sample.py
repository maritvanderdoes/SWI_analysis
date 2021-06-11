#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_clean/Running_Code.py

# Setting the directories
import _setup

# Loading custom libraries
from utils import image_lists, read_image_and_metadata
from utils.benchmarking import tic, toc

# Import additional libraries
import numpy as np
import pandas as pd
import os
import re
import glob

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import parallel, n_workers, batch_size

# Import main functions
from _main_function import main_function

# %% MAIN CODE
# Initialisation
results = []
current_res_df = []
image = 0

# list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)


#%%
ksel = 1
files = (list_mcherry[ksel], list_GFP[ksel])
(img_mcherry, img_gfp), current_res = read_image_and_metadata(files)

print(current_res)
#%%
pictureinfo = re.split('_t(\d+)_s(\d+)\..+', files[0])

#pictureinfo = re.split('_t(\d+)_s(\d+)\..+', list_mcherry[ksel])
s_info = 0
t_info = 1

print(pictureinfo[0])

#%%
s_sel = 1
t_sel = None

if t_sel == None:
    list_1=sorted(glob.glob(os.path.join(dirpath, "*"+"_s"+str(s_sel)+"*"+channel_mcherry+"*")))
    
else:
    list_1=sorted(glob.glob(os.path.join(dirpath, "*"+"_t"+str(t_sel)+"_s"+str(s_sel)+"*"+channel_mcherry+"*")))


print(list_1)

#%%
t_sel = 1
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel1 = channel_mcherry, s_sel= s_sel, t_sel = None, channel2 = channel_GFP)

print(list_mcherry)
# %%
