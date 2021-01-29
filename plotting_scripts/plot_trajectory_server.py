#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 21:53:18 2020

@author: Marit
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 19:56:49 2020

@author: Marit
"""

import os
import re
from skimage import io
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd

#import own functions
from marit_functions import img_thresholding
from marit_functions import select_zslides
from marit_functions import image_lists_mcherry_GFP_BF

#folder with images
dirpath= "/Users/Marit/Desktop/worm6"

#channels
channel_GFP='1_w1Lucas-sim-488-561'
channel_mcherry='1_w2Lucas-sim-561-488'
channel_BF= '2_Lucas-Brightfield-488-561'

# list for all channels the stk files in folder
list_mcherry, list_GFP, list_BF=image_lists_mcherry_GFP_BF(dirpath, channel_mcherry,channel_GFP,channel_BF)

zslide=10
showworm=1
showimg=[1,10,11,12]

counter=0

intensity=np.array([])
volume=np.array([])
time=np. array([])

fig1 = plt.figure(figsize=(np.size(showimg)*3,14),constrained_layout=True)
gs= fig1.add_gridspec(5,np.size(showimg))



for (file1,file2,file3) in zip(list_mcherry,list_GFP,list_BF):
    pictureinfo=re.split('_t(\d+)_s(\d+)', file1)
    
    #take only worm showworm and calculate 
    if  (int(pictureinfo[1])==showworm):
        print(file1)
        
        #open GFP and mCherry image
        img_mcherry=io.imread(file1)
        img_gfp=io.imread(file2)
        img_bf=io.imread(file3)
                

        #preprocessing and thresholding to make image binary
        img_binary=img_thresholding(img_mcherry)
        
        #select slides based on mean signal
        img_binary, idx =select_zslides(img_mcherry,img_binary)
        
        #create overlay of binary image with GFP image
        img_overlay= img_gfp*img_binary

        #time
        time=np.append(time,int(pictureinfo[2]))
        
        #calculate volume
        volume=np.append(volume, np.sum(img_binary))
        
        #calculate GFPintenstiy
        intensity=np.append(intensity, np.sum(img_overlay))
        
        if (int(pictureinfo[2]) in showimg):
            fig1.suptitle(pictureinfo[0])
            
            axMcherry=fig1.add_subplot(gs[0,counter])
            axMcherry.imshow(img_bf[zslide,:,:],cmap='gray')
            axMcherry.title.set_text("timepoint ="+ str(showimg[counter]))
            axMcherry.set_yticks([])
            axMcherry.set_xticks([])
            if counter==0:
                axMcherry.set_ylabel("Brightfield")
        
            axSegmentation=fig1.add_subplot(gs[1,counter])
            axSegmentation.imshow(img_binary[zslide,:,:],cmap="gray")
            axSegmentation.set_yticks([])
            axSegmentation.set_xticks([])
            if counter==0:
                axSegmentation.set_ylabel("segmentation")
    
            axGFP=fig1.add_subplot(gs[2,counter])
            axGFP.imshow(img_gfp[zslide,:,:],cmap="gray")
            axGFP.set_yticks([])
            axGFP.set_xticks([])
            if counter==0:
                axGFP.set_ylabel("GFP")  
            
            counter=counter+1
            print(counter)
            if counter==4:
                break
            
       

#sorting
indices=time.argsort()
sorted_time = time[indices[:]]
sorted_volume = volume[indices[:]]
sorted_intensity= intensity[indices[:]]

axVolume=fig1.add_subplot(gs[3,:])   
axVolume.plot(sorted_time,sorted_volume,".")
axVolume.set_xlabel("time (h)")
axVolume.set_ylabel("Volume (pixels)")   
axVolume.title.set_text("volume of wrom over time")    

        
axInt=fig1.add_subplot(gs[4,:]) 
axInt.plot(sorted_time,sorted_intensity/sorted_volume,".")
axInt.set_xlabel("time (h)")
axInt.set_ylabel("intensity per volume") 
axInt.title.set_text("avarage GFP intensity of worm over time")      
            
fig1.savefig(pictureinfo[0]+".png")
df = pd.DataFrame({"time" : sorted_time, "volume" : sorted_volume,"intensity": sorted_intensity/sorted_volume})
df.to_csv(pictureinfo[0]+".csv", index=False)
#save csv with quantification from hatching