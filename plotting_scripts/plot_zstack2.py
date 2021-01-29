#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 19:56:49 2020

@author: Marit
"""

from skimage import io
import re

#import numpy as np
from skimage import io

#import own functions
from marit_functions import *



#folder with images
dirpath_Segm= '/Users/Marit/Documents/SWI_oldsettings/segm'
dirpath_GFP= '/Users/Marit/Documents/SWI_oldsettings/img'

channel_BF='w2Marit-BF-488-Cam1' #BF
channel_GFP='w1Marit-488-BF-Cam0' #GFP


#showworm
showworm=2
showtime=10

# list for all channels the stk files in folder
list_BF, list_GFP=image_lists_BF_GFP(dirpath_Segm, channel_BF, dirpath_GFP, channel_GFP)


for (file1,file2) in zip(list_BF,list_GFP):
    pictureinfo=re.split('_s(\d+)_t(\d+)\..+', file1) #make to function get metadata
    
    if  (int(pictureinfo[1])==showworm and int(pictureinfo[2])==showtime):
        print(file1)
        print('yes')
        
        #open GFP and mCherry image
        img_BF = read_image(file1)
        img_GFP = read_image(file2)

        #select slides based on mean signal
        img_GFP =select_zslides2(img_GFP)
        
        # plot pictures certain z slides
        zslide=[1,5,10,15,20]
        #zslide=[0,15,-2]
        plotzslides(file1,zslide,img_GFP,img_GFP,img_GFP)
        
        
        



  