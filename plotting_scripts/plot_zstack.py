#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 19:56:49 2020

@author: Marit
"""

import os
import re

#import numpy as np
from skimage import io

#import own functions
from marit_functions import *



#folder with images
dirpath= "/Users/Marit/Documents/work/swi"

#channels
channel_mcherry='2Lucas-sim-561-488'
channel_GFP='1Lucas-sim-488-561'


#showworm
showworm=1
showtime=10

# list for all channels the stk files in folder
list_mcherry, list_GFP=image_lists_mcherry_GFP(dirpath, channel_mcherry,channel_GFP)


for (file1,file2) in zip(list_mcherry,list_GFP):
    pictureinfo=re.split('_s(\d+)_t(\d+)\..+', file1) #make to function get metadata
    
    if  (int(pictureinfo[1])==showworm and int(pictureinfo[2])==showtime):
        print(file1)
        
        #open GFP and mCherry image
        img_mcherry=io.imread(file1)
        img_gfp=io.imread(file2)

        #preprocessing and thresholding to make image binary
        img_binary=img_thresholding(img_mcherry)
        
        #select slides based on mean signal
        img_binary, idx =select_zslides(img_mcherry,img_binary)
        
        # plot pictures certain z slides
        zslide=[idx[0]-1,idx[0]+1, int((idx[0]+idx[1])/2), idx[1]-1, idx[1]+1]
        #zslide=[0,15,-2]
        plotzslides(file1,zslide,img_mcherry,img_binary,img_gfp)
        
        
        



  