#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_final/Running_Code.py

# Loading custom libraries
from utils import image_lists, file_validation, single_image_lists

# Import additional libraries
import pandas as pd
import multiprocessing

from _parameters import dirpath, outputpath, output_filename, data_format
from _parameters import channel_GFP, channel_mcherry
from _parameters import n_workers

# Import main functions
from _main_function import main_function

# %% MAIN CODE
# Selecting sample
s_sel = 35
t_sel = 114

# list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel_GFP, s_sel = s_sel, t_sel = None, channel2 = None, channel3 = None, data_format = data_format)

# File validation
file_validation(dirpath, outputpath, output_filename, list_mcherry, list_GFP)

# Actual run
if __name__=='__main__':
    p = multiprocessing.Pool(n_workers)
    results = p.starmap(main_function,zip(list_mcherry,list_GFP))
    p.close()

    # Saving files
    df = pd.DataFrame(results)
    df.to_csv(outputpath+'/'+output_filename, index = False)

