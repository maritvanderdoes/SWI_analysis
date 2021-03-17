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

def get_meta_info(path):
    pictureinfo=re.split('_s(\d+)_t(\d+)\..+', path)
    return {'Position': pictureinfo[1], 'Frame': pictureinfo[2]}

def get_meta_info_temp(path):
    pictureinfo=re.split('t(\d+)_s(\d+)_', path)
    return {'Position': pictureinfo[2], 'Frame': pictureinfo[1]}

def read_image(path):
    img = io.imread(path)
    return img.squeeze()

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

#-----------------------------------------------------------------------------
# Calculating worm properties
def calculate_worm_properties(img_binary,img_signal):
    '''
    worm_proprties_output is a function that  calculates different properties 
    (area, mean_intensity, min_intensity) of the signal images, based on the 
    biggest segmented area in the binary image

    Parameters
    ----------
    img_binary : binary image
    img_signal : signal image

    Returns
    -------
    binary_image : image that only contains the best area
    area         : area of the best area
    mean_intensity: mean intensity of the signal image in the best area
    min_intensity: minimum intensity of the signal image in the best area

    '''
    
    #select biggest area
    ccs, num_ccs = label(img_binary) #set labels in binary image
    properties=regionprops(ccs,img_signal,['area','mean_intensity','min_intensity']) #calculates the properties of the different areas
    best_region = max(properties, key=lambda region: region.area) #selects the biggest region
    
    binary_image= (ccs == best_region.label).astype(np.uint8)
    min_intensity=best_region.min_intensity
    mean_intensity=best_region.mean_intensity
    area=best_region.area
    
    return binary_image, area, mean_intensity, min_intensity

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
            #adpt = th_sel*MINSORTED[z_plane]*TH[z_plane]
            adpt = MINSORTED[px]*(1+(TH[px]-1)*th_sel)
            ADPT[px] = adpt
            THPX[z_plane,int(SORTED[z_plane,px,2]),int(SORTED[z_plane,px,1])] = \
                SORTED[z_plane,px,0]>adpt

    return THPX, SORTED, ADPT, PRETH

# Mask post-processing
def mask_postprocessing(input_mask, krn):
    for z_plane in np.arange(0,input_mask.shape[0]):
        # Erosion
        input_mask[0,:,:] = skimorph.binary_erosion(input_mask[0,:,:], krn)
        # Dilation
        input_mask[0,:,:] = skimorph.binary_dilation(input_mask[0,:,:], krn)
        # Filling holes
        input_mask[0,:,:] = scimorph.binary_fill_holes(input_mask[0,:,:])

    return input_mask