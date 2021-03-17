#%% importing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from skimage import data
import matplotlib.pyplot as plt
from skimage.util import invert
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#import cv2
from fil_finder import FilFinder2D
import astropy.units as u
import timeit

from utils import image_lists_mcherry_GFP
from utils import read_image
from utils import img_thresholding
from utils import select_zslides
from utils import calculate_worm_properties
from utils import get_meta_info_temp
from skimage.morphology import skeletonize, skeletonize_3d

from utils import tic
from utils import toc

#%% load parameters
dirpath = 'C:/Users/moraluca/Desktop/Lin28_test'
outputpath = 'C:/Users/moraluca/Desktop/Lin28_test\Output'
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#save retults
results = []
image=1

#%% to run
# list for all channels the stk files in folder
list_mcherry, list_GFP = image_lists_mcherry_GFP(dirpath, channel_mcherry, channel_GFP)

#open mcherry and segment on signal
for i,(file1, file2) in enumerate(zip(list_mcherry, list_GFP)):
    if (i==0):
        print(file1)

        img_mcherry = read_image(file1)
        img_gfp = read_image(file2)

        #preprocessing and thresholding to make image binary
        img_binary = img_thresholding(img_mcherry)
        image=np.max(img_binary,0)
  
        #find and plot skeleton
        start = tic()
        skeleton = skeletonize(image)
        skeleton3d = skeletonize_3d(img_binary)

        toc(start)
        
        fig, axes = plt.subplots(1, 3, figsize=(8, 4), sharex=True, sharey=True)
        ax = axes.ravel()
        
        ax[0].imshow(image, cmap=plt.cm.gray, interpolation='nearest')
        ax[0].set_title('original')
        ax[0].axis('off')
        
        ax[1].imshow(skeleton, cmap=plt.cm.gray, interpolation='nearest')
        ax[1].set_title('skeletonize')
        ax[1].axis('off')
        
        ax[2].imshow(np.max(skeleton3d,0), cmap=plt.cm.gray, interpolation='nearest')
        ax[2].set_title('skeletonize_3d')
        ax[2].axis('off')
        
        fig.tight_layout()
        plt.show()
        
        # find longest distance
        start = tic()

        fil = FilFinder2D(image, distance=250 * u.pc, mask=image)
        fil.preprocess_image(flatten_percent=85)
        fil.create_mask(border_masking=True, verbose=False,
        use_existing_mask=True)
        fil.medskel(verbose=False)
        fil.analyze_skeletons(branch_thresh=40* u.pix, skel_thresh=10 * u.pix, prune_criteria='length')

        plt.figure(2)
        plt.imshow(image, cmap='gray')
        plt.contour(fil.skeleton_longpath, colors='r')
        plt.axis('off')
        plt.show()
 

        toc(start)
        


        #select slides based on mean signal
        img_signal, img_binary = select_zslides(img_gfp, img_binary)
    

        #calculates properties of the segmented worm
        binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
        img_binary, img_signal)

        #create overlay of binary image with GFP image
        img_overlay = img_signal * img_binary

        #add properties in current results
        current_res = get_meta_info_temp(file2)  #get metadata
        current_res['volume'] = area  #calculate volume
        current_res['mean_intensity'] = mean_intensity
        current_res['min_intensity'] = min_intensity
        current_res['final_intensity'] = mean_intensity - min_intensity  #calculate intensity

        #save in resultarray
        results.append(current_res)


# %%

# %%
potato = fil.skeleton_longpath

np.sum(potato)
# %%
