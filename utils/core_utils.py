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

def read_image_and_metadata(path, data_format = 'ts'):
    '''
    Reads a STK or TIFF file for a path (include filename) and return
    the (3D) image and the values for position and frame.
    
    This function allows for more than on path name to be read. This 
    allows for all channels to be read simultaneously. Following Marit's 
    convention for ordering. For instance:
    (mCherry_img, GFP_img, BF_imag), metadata_out = 
        read_image_and_metadata(mCherry_path, GFP_path, BF_path)

    Filenames are assume to be of the form:
        name_sx_tx_otherinfo.extension
    Other filenames would not yield correct results.

    Version note: This function replaces read_image, get_meta_info 
    and get_meta_info_temp from the Marit module.

    Parameters
    ----------
    path : Path of the file(s) to read.
    data_format: Declaring the filename structure. By default it is 
        assumed for time to preceed position ('st' format). 
        If position preceeds time, then define data_format = 'ts'.

    Returns
    -------
    img_output: (3D) array(s) of the selected image(s). For more than
        one image, each array is contained in a tuple.
    metadata_out: Dictionary containing the position (s) and the 
        frame (t).

    '''
    # Image
    if np.ndim(path) == 0:
        img = io.imread(path)
        img_output = img.squeeze()
        
        # Metadata
        if data_format == 'st':
            pictureinfo = re.split('_s(\d+)_t(\d+)\..+', path)
            s_info = 0
            t_info = 1
        if data_format == 'ts':
            pictureinfo = re.split('t(\d+)_s(\d+)_', path)
            s_info = 1
            t_info = 0

    else:
        img_output = []
        for k in enumerate(path):
            img = io.imread(path[k[0]])
            img_output.append(img.squeeze())
        
        # Metadata
        if data_format == 'st':
            pictureinfo = re.split('_s(\d+)_t(\d+)\..+', path[1])
            s_info = 1
            t_info = 2
        if data_format == 'ts':
            pictureinfo = re.split('t(\d+)_s(\d+)_', path[1])
            s_info = 2
            t_info = 1

    metadata_out = {'Position': pictureinfo[s_info], 'Frame': pictureinfo[t_info]}

    return img_output, metadata_out

#-----------------------------------------------------------------------------
# read files from folder
def image_lists(directory, channel1, channel2 = None, channel3 = None):
    '''
    List images for different channels of the same image. Images are
    assumed to be in the same directory. Marit's convention tend to 
    order channels as mCherry, GFP and BF (see parenthesis).

    It can take a single channel and up to three.

    Version note: This function replaces image_lists_BF_GFP,
    image_lists_mcherry_GFP and image_lists_mcherry_GFP_BF from  the 
    Marit module.

    Parameters
    ----------
    dir1 : directory where all images are located
    channel1 : name of channel1 images (mcherry)
    channel2 : name of channel2 images (GFP)
    channel3 : name of channel3 images (BF)

    Returns
    -------
    list_1 : list of files with the channel1 (mCherry)
    list_2 : list of files with the channel2 (GFP)
    list_3 : list of files with the channel3 (BF)

    '''
    list_set = []
    list_1=sorted(glob.glob(os.path.join(directory, "*"+channel1+"*")))
    list_set.append(list_1)

    if channel2 is not None:
        list_2=[name.replace(channel1, channel2) for name in list_1] 
        list_set.append(list_2)

    if channel3 is not None:
        list_3= [name.replace(channel1, channel3) for name in list_1] 
        list_set.append(list_3)

    return list_set

#-----------------------------------------------------------------------------
# Calculating worm properties
def calculate_worm_properties(img_binary, img_signal):
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
def adaptive_masking(input_image, mm_th = 3, th_sel = 0.3, krn_size = 2, krn_type = 'Disk', exp_size = 1, z_threshold = 0.7, sorting = False, verbose = False):
    """
    adaptive_masking is a function that takes a 3D image (z,x,y) and it 
    proceeds to mask it using an adaptive threshold for each pixel.

    After masking, the image is then processed to:
        - Remove noisy pixels by eroding and then dilating.
        - Filling holes.
        - Bridging gaps by dilating and then eroding.
        - Removing z planes that do present areas below a threshold.

    Version note: This function replaces img_thresholding and 
    select_zslides functions from the Marit module.

    Parameters
    ----------
    input_image : 3D image of the worm
        should be a image on which the signal needs to be thresholded.

    mm_th: absolute threshold to determine background for all the z 
        planes. By default it is 1.8. However, for certain synthetic data
        it is better to define it as 1.05.

    th_sel: adaptive threshold to determine which z planes present signal.
        It has a default value of 0.4.

    krn_type: type of kernel. By default a disk kernel is used.

    krn_size: kernel size for the removal of noisy pixels. Default value 
        of 1 implies this step is not performed.

    exp_size: kernel size for bridging disconnected parts and smoothing 
        the shape. Default value of 1 implies this step is not performed.

    z_threshold:  threshold for the removal of spurious regions in certain
        z planes. Default value of 0.2 eliminates z planes with abnormal
        areas below size z_threshold*max(area_zplane).
    
    verbose: Output progression through steps and time taken. By default
        it is off.
        
    Returns
    -------
    output_mask: (3D) matrix representing the binary mask using the given
        input_image

    sorted_values: sorted intensity values for the input_image in descending
        order. It has dimensions of [zplanes, number_of_pixels,3]. 
        The third dimension contains:
        - [0] Intensity value
        - [1] X coordinate
        - [2] Y coordinate

    pixel_threshold: threshold value for each pixel in sorted_values. It
        has dimensions of [number_of_pixels] (ordered).

    pixel_range: ratio between the highest intensity value and the lowest
        intensity value for each pixel along the z axis. It has dimensions
        of [number_of_pixels] (ordered).

    area_zplane: area of the mask for each z plane.

    """
    # Debugging and benchmarking
    if verbose:
        from utils.benchmarking import tic,toc
        start0 = tic()
        print('Adaptive masking. Verbose mode')

    # Retrieving dimensions
    datdim = input_image.shape

    # Storing matrices
    sorted_values = np.zeros([datdim[0],datdim[1]*datdim[2],3])
    IMGSEL = np.zeros([datdim[1]*datdim[2],3])
    output_mask = np.zeros([datdim[0],datdim[1],datdim[2]])  
    pixel_threshold = np.zeros([datdim[1]*datdim[2]])

    # Generating matrix coordinates
    x = np.arange(0, datdim[1])
    y = np.arange(0, datdim[2])
    X, Y = np.meshgrid(x, y)

    # 1. Sorting values
    if sorting:
        for z_plane in np.arange(0,datdim[0]):
            TEMPSEL = input_image[z_plane,:,:]
            IMGSEL[:,0] = np.ndarray.flatten(TEMPSEL)
            IMGSEL[:,1] = np.ndarray.flatten(X)
            IMGSEL[:,2] = np.ndarray.flatten(Y)

            # a = np.array([[1,4,4], [3,1,1], [1,5,2], [2,1,0]]) for test
            # srt = a[a[:,0].argsort()] for test
            sorted_values[z_plane,:,:] = IMGSEL[IMGSEL[:,0].argsort()[::-1]]

    else:
        # Reshape
        sorted_values[:,:,0] = input_image.reshape(datdim[0], datdim[1]*datdim[2])
        # Introduce the X and Y
        ones_mat = np.ones([datdim[0],datdim[1],datdim[2]])
        coord_mat = X[None,:,:]*ones_mat
        sorted_values[:,:,1] = coord_mat.reshape(datdim[0], datdim[1]*datdim[2])
        coord_mat = Y[None,:,:]*ones_mat
        sorted_values[:,:,2] = coord_mat.reshape(datdim[0], datdim[1]*datdim[2])

    # Debugging and benchmarking
        if verbose:
            if sorting:
                print('Values sorted.', end = " ")
            else:
                print('Values not sorted. Matrix reshaped.', end = " ")
            stop = toc(start0)
            start = tic()

    # 2. Computing the thresholds
    MAXSORTED = np.max(sorted_values[:,:,0], axis=0)
    MINSORTED = np.min(sorted_values[:,:,0], axis=0)
    pixel_range = MAXSORTED/MINSORTED
    MTH = pixel_range > mm_th
    TH = pixel_range*MTH

    # 3. Thresholding
    px_vals = np.arange(0,datdim[1]*datdim[2])
    for px in px_vals[MTH]:
        for z_plane in np.arange(0,datdim[0]):
            #adpt = th_sel*MINSORTED[z_plane]*TH[z_plane]
            adpt = MINSORTED[px]*(1+(TH[px]-1)*th_sel)
            pixel_threshold[px] = adpt
            output_mask[z_plane,int(sorted_values[z_plane,px,2]),int(sorted_values[z_plane,px,1])] = \
                sorted_values[z_plane,px,0]>adpt

    # Debugging and benchmarking
    if verbose:
        print('Pixels thresholded.', end = " ")
        stop = toc(start)
        start = tic()

    # 4. Removing holes
    output_mask = _mask_postprocessing(output_mask, krn_size = krn_size, krn_type = krn_type, exp_size = exp_size)

    # Debugging and benchmarking
    if verbose:
        print('XY postprocessing.', end = " ")
        stop = toc(start)
        start = tic()

    # 5. Removing
    output_mask, area_zplane = _mask_refinement(output_mask, z_threshold = z_threshold)

    # Debugging and benchmarking
    if verbose:
        print('Z postprocessing.', end = " ")
        stop = toc(start)

        print('Adaptive masking finished.', end = " ")
        stop = toc(start0)

    return output_mask, sorted_values, pixel_threshold, pixel_range, area_zplane

# Supporting functions for masking
def _mask_postprocessing(input_mask, krn_size = 1, krn_type = 'Disk', exp_size = 1):
    '''
    Processes a binary mask to:
        - Remove noisy pixels by eroding and then dilating.
        - Filling holes.
        - Bridging gaps by dilating and then eroding.

    Parameters
    ----------
    input_mask : 3D binary mask of the worm.

    krn_type: type of kernel. By default a disk kernel is used.

    krn_size: kernel size for the removal of noisy pixels. Default value 
        of 1 implies this step is not performed.

    exp_size: kernel size for bridging disconnected parts and smoothing 
        the shape. Default value of 1 implies this step is not performed.
        
    Returns
    -------
    input_maks: Processed mask.

    '''
    if krn_size>1 :
        # Kernel selection
        if krn_type == 'Disk':
            krn = skimorph.disk(krn_size)
        elif krn_type == 'Square':
            krn = skimorph.square(krn_size)

        for z_plane in np.arange(0,input_mask.shape[0]):
            # Erosion
            input_mask[z_plane,:,:] = skimorph.binary_erosion(input_mask[z_plane,:,:], krn)
            # Dilation
            input_mask[z_plane,:,:] = skimorph.binary_dilation(input_mask[z_plane,:,:], krn)
            # Filling holes
            input_mask[z_plane,:,:] = scimorph.binary_fill_holes(input_mask[z_plane,:,:])

            # Expanding the mask
            
    if exp_size>1:
        # Kernel selection
        if krn_type == 'Disk':
            krn = skimorph.disk(exp_size)
        elif krn_type == 'Square':
            krn = skimorph.square(exp_size)

        for z_plane in np.arange(0,input_mask.shape[0]):
            # Dilation
                input_mask[z_plane,:,:] = skimorph.binary_dilation(input_mask[z_plane,:,:], krn)
            # Erosion
                input_mask[z_plane,:,:] = skimorph.binary_erosion(input_mask[z_plane,:,:], krn)

    return input_mask

def _mask_refinement(input_mask, z_threshold = 0.2):
    """
    Processes a binary mask to:
        - Removing z planes that do present areas below a threshold.

    Parameters
    ----------
    input_mask : 3D binary mask of the worm.

    z_threshold:  threshold for the removal of spurious regions in certain
        z planes. Default value of 0.2 eliminates z planes with abnormal
        areas below size z_threshold*max(area_zplane).
        
    Returns
    -------
    output_mask: Processed (3D) mask.

    area_zplane: area of the mask for each z plane.

    """
    # Obtaining the number of pixels
    area_zplane = np.sum(np.sum(input_mask,axis=2),axis = 1)

    # Find the maxima
    max_pixels = np.max(area_zplane)

    # Find which planes to keep
    maintain_values = area_zplane > max_pixels*z_threshold

    # Setting values to zero
    output_mask = input_mask*maintain_values[:,None,None]

    # Obtaining the number of pixels
    area_zplane = np.sum(np.sum(output_mask,axis=2),axis = 1)

    return output_mask, area_zplane