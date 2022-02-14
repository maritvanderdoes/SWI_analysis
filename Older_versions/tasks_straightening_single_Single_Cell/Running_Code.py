#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_straightening_single/Running_Code.py

# Setting the directories
from numpy.core.fromnumeric import sort
import _setup

# Loading custom libraries
from utils import image_lists, single_image_lists
from utils.benchmarking import tic, toc, saving_log
from utils import image_lists, read_image_and_metadata, adaptive_masking

# Import additional libraries
import numpy as np
import pandas as pd
import multiprocessing
import os
import matplotlib.pyplot as plt

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import parallel, n_workers, batch_size, data_format
from _parameters import debugpath, debugging

# Import main functions
from _main_function import main_function

# %%
dirpath = 'D:/SWI_Test'
outputpath = 'D:/SWI_Test/Output'
channel_GFP =   'w1Lucas-confmCherry'
channel_mcherry= 'w1Lucas-confmCherry'
data_format = 'st'

# %% MAIN CODE
# Initialisation
results = []
current_res_df = []
image = 0

# Selecting sample
s_sel = 1

# list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
# (list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel_GFP, s_sel = s_sel, t_sel = None, channel3 = None, data_format = data_format)

# %%

# Sample sel 
k = 1

# For debugging
start0 = tic()

# Determining the files to read
files = (list_mcherry[k], list_GFP[k])
print(files)

(img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)


img_binary, _  = adaptive_masking(img_mcherry, mm_th = 1.1, exp_size = 3)

#%
plt.figure()
plt.imshow(np.sum(img_binary[:,:,:],axis = 0))
plt.figure()
plt.imshow(np.max(img_mcherry[:,:,:],axis = 0))

#%% Checking different thresholds
n_pt = 21
mm_th_vec = np.linspace(1.06, 1.8, n_pt)
vals_sum = np.zeros(n_pt)

for k_th, mm_th in enumerate(mm_th_vec):
    img_binary, _  = adaptive_masking(img_mcherry, mm_th = mm_th, exp_size = 3)
    vals_sum[k_th] = np.sum(img_binary)

# %%
plt.semilogy(mm_th_vec, vals_sum)

#%%

# Main code
if parallel:        
    print('MAIN: Running in Parallel mode. Make sure you are running this in the terminal.')
    n_samples = len(list_mcherry)
    n_batches = int(np.ceil(n_samples/batch_size))
    start_parallel = tic()

    for batch_sel in np.arange(0,n_batches):
        print('MAIN: Running batch '+ str(batch_sel) + ' out of '+ str(n_batches-1))
        #batch_pos = np.arange(batch_sel*batch_size,(batch_sel+1)*batch_size)
        if batch_sel == n_batches-1:
            batch_mcherry = list_mcherry[batch_sel*batch_size:]
            batch_GFP = list_GFP[batch_sel*batch_size:]
        else:
            batch_mcherry = list_mcherry[batch_sel*batch_size:(batch_sel+1)*batch_size]
            batch_GFP = list_GFP[batch_sel*batch_size:(batch_sel+1)*batch_size]
        
        start_batch = tic()
        if __name__=='__main__':
            p = multiprocessing.Pool(n_workers)
            results = p.starmap(main_function,zip(batch_mcherry,batch_GFP))
            p.close()

        print('Batch finished. ', end = ' ')
        stop = toc(start_parallel)
        
        # Saving batch
        df = pd.DataFrame(results)
        df.to_csv(outputpath+'/Results_batch_straightened_'+str(batch_sel)+'_'+str(n_batches)+'.csv', index = False)

    print('Parallel code finished. ', end = ' ')
    stop = toc(start_parallel)

else:
    print('MAIN: Running in Sequential mode.')
    for k,files in enumerate(zip(list_mcherry, list_GFP)):
        print('Sample selected: '+str(k))
        if k == image :
            current_res = main_function(files[0],files[1])
            results.append(current_res)
            break
        if (not image) and (k != image):
            current_res = main_function(files[0],files[1])
            results.append(current_res)

    # Saving files
    df = pd.DataFrame(results)
    df.to_csv(outputpath+'/Results.csv', index = False)



# %%
