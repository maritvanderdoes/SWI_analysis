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
from scipy import ndimage




#-----------------------------------------------------------------------------
#read image and get information from image

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


def calculatesometing(a,b):
    return a+b
## warping images
# Copyright 2007 Zachary Pincus
# This file is part of CellTool.
#
# CellTool is free software; you can redistribute it and/or modify
# it under the terms of version 2 of the GNU General Public License as
# published by the Free Software Foundation.



def warp_images(from_points, to_points, images, output_region, interpolation_order = 1, approximate_grid=2):
    """Define a thin-plate-spline warping transform that warps from the from_points
    to the to_points, and then warp the given images by that transform. This
    transform is described in the paper: "Principal Warps: Thin-Plate Splines and
    the Decomposition of Deformations" by F.L. Bookstein.
    Parameters:
        - from_points and to_points: Nx2 arrays containing N 2D landmark points.
        - images: list of images to warp with the given warp transform.
        - output_region: the (xmin, ymin, xmax, ymax) region of the output
                image that should be produced. (Note: The region is inclusive, i.e.
                xmin <= x <= xmax)
        - interpolation_order: if 1, then use linear interpolation; if 0 then use
                nearest-neighbor.
        - approximate_grid: defining the warping transform is slow. If approximate_grid
                is greater than 1, then the transform is defined on a grid 'approximate_grid'
                times smaller than the output image region, and then the transform is
                bilinearly interpolated to the larger region. This is fairly accurate
                for values up to 10 or so.
    """
    transform = _make_inverse_warp(from_points, to_points, output_region, approximate_grid)
    return [ndimage.map_coordinates(np.asarray(image), transform, order=interpolation_order) for image in images]

def _make_inverse_warp(from_points, to_points, output_region, approximate_grid):
    x_min, y_min, x_max, y_max = output_region
    if approximate_grid is None: approximate_grid = 1
    x_steps = (x_max - x_min) // approximate_grid
    y_steps = (y_max - y_min) // approximate_grid
    x, y = np.mgrid[x_min:x_max:x_steps*1j, y_min:y_max:y_steps*1j]

    # make the reverse transform warping from the to_points to the from_points, because we
    # do image interpolation in this reverse fashion
    transform = _make_warp(to_points, from_points, x, y)

    if approximate_grid != 1:
        # linearly interpolate the zoomed transform grid
        new_x, new_y = np.mgrid[x_min:x_max+1, y_min:y_max+1]
        x_fracs, x_indices = np.modf((x_steps-1)*(new_x-x_min)/float(x_max-x_min))
        y_fracs, y_indices = np.modf((y_steps-1)*(new_y-y_min)/float(y_max-y_min))
        x_indices = x_indices.astype(int)
        y_indices = y_indices.astype(int)
        x1 = 1 - x_fracs
        y1 = 1 - y_fracs
        ix1 = (x_indices+1).clip(0, x_steps-1)
        iy1 = (y_indices+1).clip(0, y_steps-1)
        t00 = transform[0][(x_indices, y_indices)]
        t01 = transform[0][(x_indices, iy1)]
        t10 = transform[0][(ix1, y_indices)]
        t11 = transform[0][(ix1, iy1)]
        transform_x = t00*x1*y1 + t01*x1*y_fracs + t10*x_fracs*y1 + t11*x_fracs*y_fracs
        t00 = transform[1][(x_indices, y_indices)]
        t01 = transform[1][(x_indices, iy1)]
        t10 = transform[1][(ix1, y_indices)]
        t11 = transform[1][(ix1, iy1)]
        transform_y = t00*x1*y1 + t01*x1*y_fracs + t10*x_fracs*y1 + t11*x_fracs*y_fracs
        transform = [transform_x, transform_y]
    return transform

_small = 1e-100
def _U(x):
    return (x**2) * np.where(x<_small, 0, np.log(x))

def _interpoint_distances(points):
    xd = np.subtract.outer(points[:,0], points[:,0])
    yd = np.subtract.outer(points[:,1], points[:,1])
    return np.sqrt(xd**2 + yd**2)

def _make_L_matrix(points):
    n = len(points)
    K = _U(_interpoint_distances(points))
    P = np.ones((n, 3))
    P[:,1:] = points
    O = np.zeros((3, 3))
    L = np.asarray(np.bmat([[K, P],[P.transpose(), O]]))
    return L

def _calculate_f(coeffs, points, x, y):
    w = coeffs[:-3]
    a1, ax, ay = coeffs[-3:]
    # The following uses too much RAM:
    # distances = _U(numpy.sqrt((points[:,0]-x[...,numpy.newaxis])**2 + (points[:,1]-y[...,numpy.newaxis])**2))
    # summation = (w * distances).sum(axis=-1)
    summation = np.zeros(x.shape)
    for wi, Pi in zip(w, points):
     summation += wi * _U(np.sqrt((x-Pi[0])**2 + (y-Pi[1])**2))
    return a1 + ax*x + ay*y + summation

def _make_warp(from_points, to_points, x_vals, y_vals):
    from_points, to_points = np.asarray(from_points), np.asarray(to_points)
    err = np.seterr(divide='ignore')
    L = _make_L_matrix(from_points)
    V = np.resize(to_points, (len(to_points)+3, 2))
    V[-3:, :] = 0
    coeffs = np.dot(np.linalg.pinv(L), V)
    x_warp = _calculate_f(coeffs[:,0], from_points, x_vals, y_vals)
    y_warp = _calculate_f(coeffs[:,1], from_points, x_vals, y_vals)
    np.seterr(**err)
    return [x_warp, y_warp]
    
    

    