#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
# To run the code type on cmd/terminal:
# python tasks_clean/Running_Code.py

# Setting the directories
# import _setup

# Loading custom libraries
# from utils import image_lists, single_image_lists
# from utils.benchmarking import tic, toc

# Import additional libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
plt.rcParams["axes.grid"] = False
plt.rcParams["axes.facecolor"] = 'black'
from skimage.io import imsave, imread

# load parameters
# from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
# from _parameters import data_format

# Import main functions
# from _main_function import main_function

#%%
s_sel = 1
straightenedname = 'C:/Users/moraluca/Desktop/Lin28_test/Output/Straightened'

CURVED = []
STRAIGHT = []
for t_sel in range(0,240):
    try:
        cropped_binary = imread(straightenedname+'/A20_Cropped_mask'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        cropped_image = imread(straightenedname+'/A21_Cropped_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        area_cropped = np.sum(np.sum(cropped_binary,axis=2),axis = 1)
        best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
        # Selecting the frame
        cropped_binary = cropped_binary[best_plane,:,:]
        cropped_image = cropped_image[best_plane,:,:]

        cropped_image[cropped_image == 0] = np.nan

        straight = imread(straightenedname+'/A51_Straightened_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        straight[straight == 0] = np.nan

        CURVED = [CURVED,cropped_image]
        STRAIGHT = [STRAIGHT,straight]
    except:
        print(t_sel,end = ',')

#%%
cmap = matplotlib.cm.gray
cmap.set_bad('black',1.)
# %%
vmin = 100
vmax = 150
k_sel = 2
fig, axs = plt.subplots(1,2,figsize=(15,15))   
axs[0].imshow(CURVED[k_sel].astype(float),cmap = 'gray', vmin = vmin, vmax=vmax)
axs[0].grid(False)

axs[1].imshow(STRAIGHT[k_sel].astype(float),cmap = 'gray', vmin = vmin, vmax=vmax)
axs[1].grid(False)

# %%
imsave(straightenedname+'/Worm_comp'+'_s'+str(s_sel)+'.tiff',np.float16(STRAIGHT), check_contrast = False)

# %%


#%%
s_sel = 1
straightenedname = 'C:/Users/moraluca/Desktop/Lin28_test/Output/Straightened'

CURVED = []
STRAIGHT = []

fig, axs = plt.subplots(1,2,figsize=(15,15))   

for t_sel in range(0,240):
    try:
        cropped_binary = imread(straightenedname+'/A20_Cropped_mask'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        cropped_image = imread(straightenedname+'/A21_Cropped_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        area_cropped = np.sum(np.sum(cropped_binary,axis=2),axis = 1)
        best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
        # Selecting the frame
        cropped_binary = cropped_binary[best_plane,:,:]
        cropped_image = cropped_image[best_plane,:,:]

        cropped_image[cropped_image == 0] = np.nan

        straight = imread(straightenedname+'/A51_Straightened_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        straight[straight == 0] = np.nan

        axs[0].imshow(cropped_image.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax)
        axs[0].grid(False)

        axs[1].imshow(straight.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax)
        axs[1].grid(False)

        frames.append([fig, cmap='gray' ,animated=True)])
    except:
        print(t_sel,end = ',')

#%%
s_sel = 58
straightenedname = 'C:/Users/moraluca/Desktop/Lin28_test/Output/Straightened'
straightenedname = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210116/Output/Straightened'

CURVED = []
STRAIGHT = []

img = []
frames = []
fig = plt.figure(figsize=(15,15))   

for t_sel in range(0,240):
    try:
        cropped_binary = imread(straightenedname+'/A20_Cropped_mask'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        cropped_image = imread(straightenedname+'/A21_Cropped_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        area_cropped = np.sum(np.sum(cropped_binary,axis=2),axis = 1)
        best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
        # Selecting the frame
        cropped_binary = cropped_binary[best_plane,:,:]
        cropped_image = cropped_image[best_plane,:,:]

        cropped_image[cropped_image == 0] = np.nan

        straight = imread(straightenedname+'/A51_Straightened_data'+'_t'+str(t_sel)+'_s'+str(s_sel)+'.tiff')
        straight[straight == 0] = np.nan

        # asd = plt.imshow(cropped_image.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax ,animated=True)
        # asd = plt.imshow(cropped_image.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax ,animated=True)
        asd = plt.imshow(straight.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax ,animated=True)
        #axs[0].grid(False)

        #axs[1].imshow(straight.astype(float),cmap = 'gray', vmin = vmin, vmax=vmax)
        #axs[1].grid(False)

        frames.append([asd])
    except:
        print(t_sel,end = ',')

print('\n Sorted') 



ani = animation.ArtistAnimation(fig, frames, interval=50, blit=False,
                                repeat_delay=1000)
ani.save(straightenedname+'\im.tiff')
# plt.show()
# %%
