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
from fil_finder import FilFinder2D
import astropy.units as u
from utils.benchmarking import tic,toc


#-----------------------------------------------------------------------------
# Loading Datasets

def get_meta_info(path, data_format = 'st'):
    if data_format == 'st':
        pictureinfo = re.split('_s(\d+)_t(\d+)\..+', path)
        s_info = 1
        t_info = 2
    if data_format == 'ts':
        pictureinfo = re.split('t(\d+)_s(\d+)_', path)
        s_info = 2
        t_info = 1

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

def calculate_worm_properties(img_binary, img_signal):
    '''
    worm_proprties_output is a function that  calculates different properties 
    (volume, mean_intensity, min_intensity) of the signal images, based on the 
    biggest segmented area in the binary image

    Parameters
    ----------
    img_binary : binary image
    img_signal : signal image

    Returns
    -------
    binary_image : image that only contains the best area
    volume         : area of the best area
    mean_intensity: mean intensity of the signal image in the best area
    min_intensity: minimum intensity of the signal image in the best area

    '''
    
    #select biggest area
    ccs, num_ccs = label(img_binary) #set labels in binary image
    properties=regionprops(ccs,img_signal,['area','mean_intensity','min_intensity']) #calculates the properties of the different areas
    best_region = max(properties, key=lambda region: region.area) #selects the biggest region
    
    #outputs
    cropped_binary=best_region.image
    cropped_image=best_region.intensity_image

    min_intensity=best_region.min_intensity
    mean_intensity=best_region.mean_intensity
    volume=best_region.area

    #possible outputs
    centroid = best_region.centroid
    binary_image= (ccs == best_region.label).astype(np.uint8)
    metrics = (cropped_image, cropped_binary, volume, mean_intensity, min_intensity)
    
    return  metrics

def create_skeleton(image, mask, verbose=False):
    """ function uses the package https://fil-finder.readthedocs.io/en/latest/ 
    to create skeleton and prune the she skeleton based on the colored image. 
    Images can be in 3D, but it creates a skeleton based on the maximum intensity projection of the 3D images.
    So the function generates a skeleton in X and Y coordinates. 
    For extra accuracy, the mask is also loaded. 
    

    Args:
        image ([3D image]): [Image with signal used for segmentation]
        mask ([3D image]): [generated mask for the image]
        verbose (bool, optional): [When true, skeletion will be plotted in the image ]. Defaults to False.

    Returns:
        X: 
        Y:
    
    """
    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Computing the skeleton. Verbose mode.', end = " ")

    zslide_to_focus = np.int(mask.shape[0]/2)

    image2D = image[zslide_to_focus,:,:]
    mask2D = mask[zslide_to_focus,:,:]
 

    fil = FilFinder2D(image2D, mask=mask2D) #creates a class, add beamwith?!
    # idea, put the fill finder on one plane where the worm is in focus 
    fil.preprocess_image()
    fil.create_mask(use_existing_mask=True)
    fil.medskel(verbose=False)
    fil.analyze_skeletons(branch_thresh=40* u.pix, skel_thresh=10 * u.pix, prune_criteria='length')

    skeletonized_image=fil.skeleton_longpath
    [X,Y]=np.where(skeletonized_image==1)

    if verbose:
        stop = toc(start)

    return X, Y

def head2tail_masking(X, Y, dx, dy, img_binary, cut_th=0.2, verbose=False):

    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Cutting the worm at '+str(cut_th*100)+'%. Verbose mode.', end = " ")

    # # check for sizes
    # if np.ndim(img_binary)==2:
    #     shapez = 1
    #     shapex = img_binary.shape[0]
    #     shapey = img_binary.shape[1]
    # elif np.ndim(img_binary)==3:
    #     shapez = img_binary.shape[0]
    #     shapex = img_binary.shape[1]
    #     shapey = img_binary.shape[2]

    shapex = img_binary.shape[-2]
    shapey = img_binary.shape[-1]
    
    # if verbose:
    #     print("Number of z_planes: "+number_z)

    # binary_to_test
    # Define the points
    grid_y, grid_x = np.meshgrid(np.arange(0,shapey),np.arange(0,shapex))

    # Find the first line
    points2mask = np.linspace(0,1,np.shape(X)[0])

    # Find the upperline
    upper = np.sum(points2mask<cut_th)
    lower = np.sum(points2mask<(1-cut_th))

    # Reference point
    ref_up = X[upper]*dx[upper]+Y[upper]*dy[upper]
    ref_low = X[lower]*dx[lower]+Y[lower]*dy[lower]

    # define the upper points to be zero
    fun_u = (grid_y*dy[upper]+grid_x*dx[upper])
    fun_l = (grid_y*dy[lower]+grid_x*dx[lower])
    binary_upper = fun_u>ref_up
    binary_lower = fun_l<ref_low

    #new image
    binary_grid = binary_upper*binary_lower

    # Final product
    binary_new = binary_grid*img_binary
    # if np.ndim(img_binary)==2:
    #     binary_new = binary_grid*img_binary
    # elif np.ndim(img_binary)==3:
    #     binary_new = binary_grid[None,:,:]*img_binary


    # plt.figure()
    # plt.contour((fun_u)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')
    # plt.figure()
    # plt.imshow((mat_upper)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')

    # plt.figure()
    # plt.contour((fun_l)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')
    # plt.figure()
    # plt.imshow((mat_lower)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')


    #adapt if it is 2D!!!
    if verbose:
        stop = toc(start)
        
    return binary_new

def _arc_length(x, y):
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)
    
    return arc