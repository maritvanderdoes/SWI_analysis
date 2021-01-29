#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:03:25 2020

@author: Marit
"""

import re
from skimage import io
import matplotlib.pyplot as plt
import numpy as np

#import own functions
from marit_functions import *


#folder with images
dirpath= "/Users/Marit/Documents/work/swi"

#channels
channel_GFP='1Lucas-sim-488-561'
channel_mcherry='2Lucas-sim-561-488'

#showzslides
zslide=15
showworm=1
showimg=[2,3,4,5]
counter=0


# list for all channels the stk files in folder
list_mcherry, list_GFP=image_lists_mcherry_GFP(dirpath, channel_mcherry,channel_GFP)

fig, ax = plt.subplots(3,np.size(showimg), figsize=(np.size(showimg)*3,9))
        

for (file1,file2) in zip(list_mcherry,list_GFP):
    pictureinfo=re.split('_s(\d+)_t(\d+)\..+', file1)
    
    if  (int(pictureinfo[1])==showworm and int(pictureinfo[2]) in showimg):
        print(file1)
    
        #open GFP and mCherry image
        img_mcherry=io.imread(file1)
        img_gfp=io.imread(file2)

        #preprocessing and thresholding to make image binary
        img_binary=img_thresholding(img_mcherry)
        
        #select slides based on mean signal
        img_binary, idx =select_zslides(img_mcherry,img_binary)
        
        
        # plot pictures certain z slides
        fig.suptitle(file1)
        ax[0,counter].imshow(img_mcherry[zslide,:,:],cmap='gray')
        ax[0,counter].title.set_text("timepoint ="+ str(showimg[counter]))
        ax[0,counter].set_yticks([])
        ax[0,counter].set_xticks([])
        ax[0,0].set_ylabel("mCherry")
        
          
        ax[1,counter].imshow(img_binary[zslide,:,:],cmap="gray")
        ax[1,counter].set_yticks([])
        ax[1,counter].set_xticks([])
        ax[1,0].set_ylabel("segmentation")
        
        ax[2,counter].imshow(img_gfp[zslide,:,:],cmap="gray")
        ax[2,counter].set_yticks([])
        ax[2,counter].set_xticks([])
        ax[2,0].set_ylabel("GFP")  
        counter=counter+1
        
