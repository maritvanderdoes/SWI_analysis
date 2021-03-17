#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Importing required libraries
import skimage.filters as skf
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import glob
import os
import re
from skimage import io
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops
import scipy.ndimage.morphology as scimorph
import skimage.morphology as skimorph

#-----------------------------------------------------------------------------
# Loading Datasets

def get_meta_info(path, data_format = 'st'):
    if data_format == 'st':
        pictureinfo = re.split('_s(\d+)_t(\d+)\..+', path)
        s_info = 0
        t_info = 1
    if data_format == 'ts':
        pictureinfo = re.split('t(\d+)_s(\d+)_', path)
        s_info = 1
        t_info = 0

    return {'Position': pictureinfo[s_info], 'Frame': pictureinfo[t_info]}


def get_meta_info_temp(path):
    pictureinfo = re.split('t(\d+)_s(\d+)_', path)
    return {'Position': pictureinfo[2], 'Frame': pictureinfo[1]}

    
def read_image(path):
    if np.ndim(path) == 0:
        img = io.imread(path)
        img_output = img.squeeze()
    else:
        img_output = []
        for k in enumerate(path):
            img = io.imread(path[k[0]])
            img_output.append(img.squeeze())

    return img_output

#-----------------------------------------------------------------------------
# read files from folder
def image_lists_BF_GFP(dir1, channel1, dir2):
    """
    lists GFP images in dir1/channel1 and find corresponding GFP images in dir2

    Parameters
    ----------
    dir1 : directory where the GFP images are located 
    channel1 : name of channel1 images(GFP) 
    dir2 : directory where segmented images are located

    Returns
    -------
    list_1 : list of files in GFP 
    list_2 : list of files in BF

    """
    #creating list1 with GFP images           
    list_1=sorted(glob.glob(os.path.join(dir1, "*"+channel1+"*")))
    
    #checks if the channel pattern is present in the directory
    if not list_1:
        raise FileNotFoundError(
         'No files in image_folder match the given image_file_pattern={}'
         .format(channel1))
                
    #creating list2 with segmented images           
    list_2=[]
    for name in list_1:
        match=re.split('(\S*_)(\S*)(_s\S*)', os.path.basename(name)) #space=[0],exp_name=[1],channelname=[2],wrom&tp=[3]
        #finds the matched image pairs for exp_name and worm&tp in dir2
        matching_BF=glob.glob(os.path.join(dir2, match[1]+"*"+ match[3])) 
        
        #checks if the segmented image pair is present
        if not matching_BF:
            raise FileNotFoundError(
                'the segmented pair was not found for ={}'
                .format(name))
        
        #if image pair is present, first entry of the mathing file [0] is added to list2
        list_2.append(matching_BF[0])
    return list_1, list_2

def image_lists_mcherry_GFP(directory, channel1,channel2):
    """
    lists mcherry images in dir1/channel1 and corresponding GFP images in dir2/channel2

    Parameters
    ----------
    dir1 : directory where brightfield images are located
    channel1 : name of channel1 images(BF)
    dir2 : directory where GFP images are located
    channel2 : name of channel2 images(GFP)

    Returns
    -------
    list_1 : list of files in BF
    list_2 : list of files in GFP

    """
    
    list_1=sorted(glob.glob(os.path.join(directory, "*"+channel1+"*")))
    list_2=[name.replace(channel1, channel2) for name in list_1]     
    return list_1, list_2

def image_lists_mcherry_GFP_BF(directory, channel1,channel2,channel3):
    """
    lists mcherry images in dir1/channel1 and corresponding GFP and BF images in channel2 and 3 

    Parameters
    ----------
    dir1 : directory where all images are located
    channel1 : name of channel1 images(mcherry)
    channel2 : name of channel2 images(GFP)
    channel3 : name of channel2 images(BF)

    Returns
    -------
    list_1 : list of files in mcherry
    list_2 : list of files in GFP
    list_3 : list of files in BF

    """
    
    list_1=sorted(glob.glob(os.path.join(directory, "*"+channel1+"*")))
    list_2=[name.replace(channel1, channel2) for name in list_1]     
    list_3= [name.replace(channel1, channel3) for name in list_1]  
    return list_1, list_2, list_3

#----------------------------------------------------------------------------
# selecting good chambers and  slides in focus

def select_zslides(signal_image,binary_image):
    """
    select z slides that are in focus, where the mean intensity of every plane is above the avarage over all planes
    
    Parameters
    ----------
    signal_image : signal that determines where the worm is in focus
    binary_image : binary image
    
    Returns
    -------
    selected_slides_signal : signal image with only sz-lides that are in focus
    selected_slides_binary : binary image with only z-slides that are in focus

    """
    
    mean_img = np.mean(signal_image, 
                       axis=tuple(ax for ax in range(1, signal_image.ndim))) #calculates mean of every zslide
    threshold = (np.max(mean_img)+np.min(mean_img))/2 #calulates the mean of mean, and determines the treshold

    start, stop = find_best_interval(mean_img, threshold) #finds the best interval where intensity is above avarage strashold
    
    selected_slides_signal = signal_image[start:stop] # removes the slides that are out of focus
    selected_slides_binary = binary_image[start:stop] 
    
    # plt.figure()
    # plt.plot(mean_img)
    # plt.hlines(threshold,start,stop,'r')
    # plt.plot([start,start],[0,threshold],'--r')
    # plt.plot([stop,stop],[0,threshold],'--r')
    # plt.title("mean intensity of zlides per plane")
    # plt.xlabel("z-slide")
    # plt.ylabel("mean intensity")
            
    return selected_slides_signal, selected_slides_binary

def find_best_interval(values, threshold):
    '''returns the indices of the interval over the given values that maximizes

        \sum_{i \in [start:stop)} sign(values[i] - threshold)

    Please note that the intervals stop is always exclusive. You can
    therefore use values[start:stop] to access all values of the
    interval.

    '''
    class Interval:
        def __init__(self, start, stop, score=0.0):
            self.start = start
            self.stop = stop
            self.score = score

        def reset(self, idx):
            self.start = idx
            self.stop = idx
            self.score = 0.0

    assert values.ndim == 1

    current = Interval(0, 0)
    best = Interval(0, 0)

    for idx, val in enumerate(np.sign(values - threshold)):
        current.score += val
        if current.score < 0.0:
            current.reset(idx + 1)

        elif current.score > best.score:
            best = Interval(current.start, idx + 1, current.score)

    return best.start, best.stop

#----------------------------------------------------------------------------
# make image binary

def img_thresholding(input_image):
    """
    img_thresholding is a function that applies gausian filter,
    creates a threshold based on a treshold method (otsu),
    and plots the histogram of intensities with the calculated threshold

    Parameters
    ----------
    input_image : 3D image of the worm
        should be a image on which the signal needs to be thresholded.
        
    Returns
    -------
    binary_image : 3D binary image of the worm
        is a mask of the image in binary.

    """
    
    img_gaus=gaussian_filter(input_image, sigma=5)
    
    treshold = skf.threshold_otsu(img_gaus)
    binary_image= img_gaus > treshold
    
    histogram, bin_edges = np.histogram(img_gaus, bins=256)
        
    # plt.figure()
    # plt.plot(bin_edges[0:-1], histogram) 
    # plt.axvline(treshold,color='r')
    # plt.title("Grayscale Histogram")
    # plt.xlabel("grayscale value")
    # plt.ylabel("pixels")
    # plt.yscale('log')
    # plt.text(treshold+20,50, 'threshold = '+ str(treshold))

    return binary_image