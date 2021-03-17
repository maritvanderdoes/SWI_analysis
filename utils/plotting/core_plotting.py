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

#----------------------------------------------------------------------------   
# plotting functions

def plotzslides(title,zslides,img1,img2,img3):
    """
    plotzslides creates a figure that shows image 1 in row 1, image 2 in row 2
    and image 3 in row 3. It shows in the column the z slides that selected, and
    give the title as the name of the figure

    Parameters
    ----------
    title : string
        DESCRIPTION.
    zslides : array
        DESCRIPTION.
    img1 : 3D mCherry image
        DESCRIPTION.
    img2 : 3D binary image
        DESCRIPTION.
    img3 : 3D GFP image
        DESCRIPTION.

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots(3,np.size(zslides), figsize=(np.size(zslides)*3,9))
    fig.suptitle(title)
        
    minimum_img1=np.min(img1)
    maximum_img1=np.max(img1)
           
    for counter,z in enumerate(zslides):
        ax[0,counter].imshow(img1[z,:,:],cmap='gray',vmin=minimum_img1, vmax=maximum_img1)
        ax[0,counter].title.set_text("zslide="+str(z))
        ax[0,counter].set_yticks([])
        ax[0,counter].set_xticks([])
        ax[0,0].set_ylabel("mCherry")
                  
        ax[1,counter].imshow(img2[z,:,:],cmap="gray")
        ax[1,counter].set_yticks([])
        ax[1,counter].set_xticks([])
        ax[1,0].set_ylabel("segmentation")
            
        ax[2,counter].imshow(img3[z,:,:],cmap="gray")
        ax[2,counter].set_yticks([])
        ax[2,counter].set_xticks([])
        ax[2,0].set_ylabel("GFP")  
           
    plt.tight_layout()

# Summary masking
def masking_summary(PRETH, SORTED, ADPT, mm_th):
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

    # SORTED distribution
    plt.figure()
    plt.imshow(SORTED[:,:,0],aspect='auto')
    plt.xlabel('Pixel index')
    plt.ylabel('Z-plane')

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

    
#-----------------------------------------------------------------------------
# Validating synthetic datasets
def dataset_comparison(ground_truth, noisy_input, output_image):
        # Plotting results
        fig, axs = plt.subplots(2, 3, sharex = True)
        # Ground truth
        axs[0,0].imshow(noisy_input[int(np.round(noisy_input.shape[0]/2)),:,:])   
        axs[0,0].set_title('Noisy input') 
        axs[1,0].imshow(noisy_input[:,int(np.round(noisy_input.shape[2]/2)),:],aspect = 'auto')   

        # Noisy input
        axs[0,1].imshow(ground_truth[int(np.round(ground_truth.shape[0]/2)),:,:])      
        axs[0,1].set_title('Ground truth') 
        axs[1,1].imshow(ground_truth[:,int(np.round(ground_truth.shape[1]/2)),:],aspect = 'auto')  

        # Masked
        axs[0,2].imshow(output_image[int(np.round(output_image.shape[0]/2)),:,:])      
        axs[0,2].set_title('Masked Data') 
        axs[1,2].imshow(output_image[:,int(np.round(output_image.shape[2]/2)),:],aspect = 'auto')  
        
        # Comparison
        fig, axs = plt.subplots(2, 1, sharex = True)
        comparison = np.abs(output_image-ground_truth)
        axs[0].imshow(np.sum(comparison[:,:,:],axis = 0))      
        axs[0].set_title('Projected Difference') 
        axs[1].imshow(np.sum(comparison[:,:,:],axis = 2),aspect = 'auto')      

        # Summary
        print('Non-matching pixels = '+str(np.sum(np.abs(output_image-ground_truth))))
