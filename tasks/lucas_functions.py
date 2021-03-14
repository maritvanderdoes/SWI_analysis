#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
#import cv2


#-----------------------------------------------------------------------------
# Generating synthetic datasets

# Generating a simple ball
def generating_ball(nx, ny, nz, sx, sy, sz, th_level, n_level, b_level):
    x = np.linspace(-1, 1, nx) + sx
    y = np.linspace(-1, 1, ny) + sy
    z = np.linspace(-1, 1, nz) + sz

    Z = np.zeros([nz,nx,ny])
    X = np.zeros([nz,nx,ny])
    Y = np.zeros([nz,nx,ny])

    for z_p in range(0,nz):
        Xt, Yt = np.meshgrid(x, y)
        X[z_p,:,:] = Xt
        Y[z_p,:,:] = Yt
        Z[z_p,:,:] = z[z_p]

    mask_data = (np.sqrt(Z**2+X**2+Y**2)<th_level) 
    test_data = mask_data + n_level*np.random.randn(nz,nx,ny) + b_level

    return mask_data, test_data, X, Y, Z

#-----------------------------------------------------------------------------
# Validating synthetic datasets
def dataset_comparison(ground_truth, noisy_input, output_image):
        # Plotting results
        fig, axs = plt.subplots(2, 3, sharex = True)
        # Ground truth
        axs[0,0].imshow(noisy_input[int(np.round(noisy_input.shape[0]/2)),:,:])   
        axs[0,0].set_title('Noisy input') 
        axs[1,0].imshow(noisy_input[:,int(np.round(noisy_input.shape[2]/2)),:])   

        # Noisy input
        axs[0,1].imshow(ground_truth[int(np.round(ground_truth.shape[0]/2)),:,:])      
        axs[0,1].set_title('Ground truth') 
        axs[1,1].imshow(ground_truth[:,int(np.round(ground_truth.shape[1]/2)),:])  

        # Masked
        axs[0,2].imshow(output_image[int(np.round(output_image.shape[0]/2)),:,:])      
        axs[0,2].set_title('Masked Data') 
        axs[1,2].imshow(output_image[:,int(np.round(output_image.shape[2]/2)),:])  

        # Summary
        print('Non-matching pixels = '+str(np.sum(np.abs(THPX-mask_data))))

#-----------------------------------------------------------------------------
# Masking

# Adaptive masking
def adaptive_masking(input_image, mm_th, th_sel):
    """
    adaptive_masking is a function that takes a 3D image (z,x,y) and it proceeds
    to mask it using an adaptive threshold for each pixel.

    This function should replace img_thresholding and select_zslides functions 
    from the Marit package.

    Parameters
    ----------
    input_image : 3D image of the worm
        should be a image on which the signal needs to be thresholded.

    mm_th: the absolute threshold

    th_sel: the adaptive threshold
        
    Returns
    -------
    THPX : 3D binary image of the worm
        is a mask of the image in binary.

    """
    # Retrieving dimensions
    datdim = input_image.shape

    # Storing matrices
    SORTED = np.zeros([datdim[0],datdim[1]*datdim[2],3])
    IMGSEL = np.zeros([datdim[1]*datdim[2],3])
    THPX = np.zeros([datdim[0],datdim[1],datdim[2]])  
    ADPT = np.zeros([datdim[1]*datdim[2]])

    # Generating matrix coordinates
    x = np.arange(0, datdim[1])
    y = np.arange(0, datdim[2])
    X, Y = np.meshgrid(x, y)

    # 1. Sorting values
    for z_plane in np.arange(0,datdim[0]):
        TEMPSEL = input_image[z_plane,:,:]
        IMGSEL[:,0] = np.ndarray.flatten(TEMPSEL)
        IMGSEL[:,1] = np.ndarray.flatten(X)
        IMGSEL[:,2] = np.ndarray.flatten(Y)

        # a = np.array([[1,4,4], [3,1,1], [1,5,2], [2,1,0]]) for test
        # srt = a[a[:,0].argsort()] for test
        SORTED[z_plane,:,:] = IMGSEL[IMGSEL[:,0].argsort()[::-1]]

    # 2. Computing the thresholds
    MAXSORTED = np.max(SORTED[:,:,0], axis=0)
    MINSORTED = np.min(SORTED[:,:,0], axis=0)
    PRETH = MAXSORTED/MINSORTED
    MTH = PRETH > mm_th
    TH = PRETH*MTH

    # 3. Thresholding
    for px in np.arange(0,np.sum(MTH)):
        for z_plane in np.arange(0,datdim[0]):
            #ADPT[px] = th_sel*MINSORTED[px]*TH[px]
            ADPT[px] = MINSORTED[px]*(1+(TH[px]-1)*th_sel)
            THPX[z_plane,int(SORTED[z_plane,px,2]),int(SORTED[z_plane,px,1])] = \
                SORTED[z_plane,px,0]>ADPT[z_plane]

    # 4. Refining

    return THPX, SORTED, ADPT, PRETH

# Summary masking
def masking_summary(PRETH, ADPT, mm_th):
    # Ordered pixels
    fig, axs = plt.subplots(2, 1, sharex = True)   
    axs[0].plot(PRETH,'r')
    axs[0].plot([0, PRETH.shape[0]],[mm_th, mm_th],'k:')
    axs[0].set_xlabel('Pixel index')
    axs[0].set_ylabel('Intensity range')

    # Treshold
    axs[1].plot(ADPT,'r')
    axs[1].set_xlabel('Pixel index')
    axs[1].set_ylabel('Threshold value')

    # Examples
    fig, axs = plt.subplots(1,3)
    axs[0].plot(SORTED[:,0,0])
    axs[0].plot([0, SORTED.shape[0]],[ADPT[0], ADPT[0]])
    axs[0].set_xlabel('Z plane')
    axs[0].set_title('Highest pixel')

    midpx = np.sum(PRETH>mm_th)-1
    axs[1].plot(SORTED[:,midpx,0])
    axs[1].plot([0, SORTED.shape[0]],[ADPT[midpx], ADPT[midpx]])
    axs[1].set_xlabel('Z plane')
    axs[1].set_title('Last pixel')

    axs[2].plot(SORTED[:,-1,0])
    axs[2].plot([0, SORTED.shape[0]],[ADPT[-1], ADPT[-1]])
    axs[2].set_xlabel('Z plane')
    axs[2].set_title('Lowest pixel')
    #axs[2].set_ylim([np.min(SORTED[:,-1,0]), np.max(SORTED[:,-1,0])])
