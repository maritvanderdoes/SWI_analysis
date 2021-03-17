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