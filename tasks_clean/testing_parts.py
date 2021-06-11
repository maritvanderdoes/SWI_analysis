#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_clean/testing_parts.py

# Setting the directories
import _setup

# Loading custom libraries
from utils import image_lists
from utils.benchmarking import tic, toc

# Import additional libraries
import numpy as np
import pandas as pd
import multiprocessing

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
ksel = 0
current_res = main_function(list_mcherry[ksel],list_GFP[ksel])