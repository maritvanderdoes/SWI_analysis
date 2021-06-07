#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Importing required libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import glob
import os
import re
import networkx as nx
from fil_finder import FilFinder2D
import astropy.units as u


from skimage import io
import skimage.filters as skf
import skimage.morphology as skimorph
from skimage.morphology import skeletonize

from skimage.measure import regionprops

from scipy.ndimage.measurements import label 
import scipy.ndimage.morphology as scimorph

from scipy import interpolate

from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import KNeighborsRegressor

from utils.benchmarking import tic,toc



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
def crop_image(img_binary, img_signal):
    """crop worm is a function that gets a cropped image of the biggest
    segmented region in the image.

    Args:
        img_binary (3D array): [description]
        img_signal (3D array): [description]

    Returns:
        cropped_image (3D arary): [description]
        cropped_binary (3D arary): [description]
    """

    
    ccs, num_ccs = label(img_binary) #set labels in binary image
    properties=regionprops(ccs,img_signal,['area']) #calculates the properties of the different areas

    if (len(properties)==0):
        cropped_binary=img_binary
        cropped_image=img_signal
    else:
        best_region = max(properties, key=lambda region: region.area) #selects the biggest region
        cropped_binary=best_region.image
        cropped_image=best_region.intensity_image

    return cropped_binary,cropped_image

def calculate_worm_properties(img_binary, img_signal, verbose = False):
    """calculate worm properties is a function that calculate the area of the segmented area
    in the binary image, and calculate the mean intensity of this segmented area in the image signal.
    This mean intensity is substracted by the minimum intensity as background substraction

    Args:
        img_binary (3D array): [description]
        img_signal (3D array): [description]

    Returns:
        mean_intensity: [description]
        volume        : 
    """
    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Calculating worm properties. Verbose mode.', end = " ")

    ccs, num_ccs = label(img_binary) #set labels in binary image

    if (num_ccs==1):
        properties=regionprops(ccs,img_signal,['area','mean_intensity','min_intensity'])
 
        min_intensity=properties[0].min_intensity
        mean_intensity=properties[0].mean_intensity-properties[0].min_intensity
        volume=properties[0].area

        # Debugging and benchmarking
        if verbose:
            stop = toc(start)

        return (mean_intensity, volume)
    
    else:
        # Debugging and benchmarking
        if verbose:
            stop = toc(start)

        print('Multiple areas are selected, segmentation not good.')
        return (np.nan, np.nan)

#-----------------------------------------------------------------------------
# Masking

# Adaptive masking
def adaptive_masking(input_image, mm_th = 3, th_sel = 0.3, krn_size = 2,
    krn_type = 'Disk', exp_size = 3, fill_holes = True, z_threshold = 0.7, 
    sorting = False, verbose = False):
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

    fill_holes: Are holes filled?

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
        start0 = tic()
        print('Adaptive masking. Verbose mode')

    # Retrieving dimensions
    datdim = input_image.shape

    # Storing matrices
    sorted_values = np.zeros([datdim[0],datdim[1]*datdim[2],3])
    IMGSEL = np.zeros([datdim[1]*datdim[2],3])
    output_mask = np.zeros([datdim[0],datdim[1],datdim[2]])  
    pixel_threshold = np.zeros([datdim[1]*datdim[2]])

    # 1. Sorting values
    if sorting:
        # Generating matrix coordinates
        x = np.arange(0, datdim[1])
        y = np.arange(0, datdim[2])
        X, Y = np.meshgrid(x, y)

        for z_plane in np.arange(0,datdim[0]):
            TEMPSEL = input_image[z_plane,:,:]
            IMGSEL[:,0] = np.ndarray.flatten(TEMPSEL)
            IMGSEL[:,1] = np.ndarray.flatten(X)
            IMGSEL[:,2] = np.ndarray.flatten(Y)

            # Sorting values
            sorted_values[z_plane,:,:] = IMGSEL[IMGSEL[:,0].argsort()[::-1]]
            # a = np.array([[1,4,4], [3,1,1], [1,5,2], [2,1,0]]) for test
            # srt = a[a[:,0].argsort()] for test

    else:
        # Reshape
        sorted_values[:,:,0] = input_image.reshape(datdim[0], datdim[1]*datdim[2])

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
    pixel_threshold = MINSORTED*(1+(TH-1)*th_sel)

    if sorting:
        px_vals = np.arange(0,datdim[1]*datdim[2])
        for px in px_vals[MTH]:
            for z_plane in np.arange(0,datdim[0]):
                #adpt = th_sel*MINSORTED[z_plane]*TH[z_plane]
                #adpt = MINSORTED[px]*(1+(TH[px]-1)*th_sel)
                #pixel_threshold[px] = adpt
                output_mask[z_plane,int(sorted_values[z_plane,px,2]),int(sorted_values[z_plane,px,1])] = \
                    sorted_values[z_plane,px,0] > pixel_threshold[px] #formerly adpt

    else:
        output_mask = (sorted_values[:,:,0]>pixel_threshold)*MTH
        output_mask = output_mask.reshape(datdim)

    # Debugging and benchmarking
    if verbose:
        print('Pixels thresholded.', end = " ")
        stop = toc(start)
        start = tic()

    # 4. Removing holes
    output_mask = _mask_postprocessing(output_mask, krn_size = krn_size, krn_type = krn_type, exp_size = exp_size, fill_holes = fill_holes)

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

    return output_mask, (sorted_values, pixel_threshold, pixel_range, area_zplane)

# Supporting functions for masking
def _mask_postprocessing(input_mask, krn_size = 1, krn_type = 'Disk', exp_size = 1, fill_holes = True):
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

    fill_holes: Are holes filled?
        
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

    # Filling holes
    if fill_holes:
        for z_plane in np.arange(0,input_mask.shape[0]):
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

#------------------------------------------------------------------------------
# straightening


def create_skeleton(mask, verbose=False):
    """[summary]

    Args:
        mask (3D image): 3D binary image/mask of the worm
        verbose (bool, optional): When verbose is true, it plots the skeleton in red over the masked image. Defaults to False.

    Returns:
        X: unsorted x coordinates of skeleton
        Y: unsorted y coordinates of skeleton
    """
    zslide_to_focus=np.int(mask.shape[0]/2)

    image=mask[zslide_to_focus,:,:]
    skeletonized_image=skeletonize(image)
    [X,Y]=np.where(skeletonized_image==1)
    
    if (verbose==True):
        plt.figure()
        plt.imshow(image, cmap='gray')
        plt.contour(skeletonized_image, colors='r')
        plt.axis('off')
        plt.title('check skeleton')
        plt.show()
    
    return X,Y

def create_spline(X,Y, sampling_space=1 , s=100, k=3, n_neighbors = 2, verbose=False):
    """Creates a spline through the X and Y coordinates by a number of splinepoints.
    This gives as output the rescaled x and y, and the derivatives of the line.
    Idea explained from 
    https://stackoverflow.com/questions/37742358/sorting-points-to-form-a-continuous-line

    Args:
        X ([1D array]): xcoordinates through with spline need to be splitted
        Y ([1D array]): ycoordinates through with spline need to be splitted
        verbose (bool, optional): [description]. Defaults to False.

    Returns:
        x_new: 
        y_new:
        dx:
        dy:
    """
    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Creating spline. Verbose mode.', end = " ")

    #variables
    nsplinepoints=int(X.shape[0]/sampling_space)
    inputarray = np.c_[X, Y]

    # runs a nearest neighbors algorithm on the coordinate array
    clf = NearestNeighbors(n_neighbors=n_neighbors).fit(inputarray)
    G = clf.kneighbors_graph()
    T = nx.from_scipy_sparse_matrix(G)

    # sorts coordinates according to their nearest neighbors order
    order = list(nx.dfs_preorder_nodes(T, 0))
    xx = X[order]
    yy = Y[order]

    # Loops over all points in the coordinate array as origin, determining which results in the shortest path
    paths = [list(nx.dfs_preorder_nodes(T, i)) for i in range(len(inputarray))]
    mindist = np.inf
    minidx = 0

    for i in range(len(inputarray)):
        p = paths[i]           # order of nodes
        ordered = inputarray[p]    # ordered nodes
        # find cost of that order by the sum of euclidean distances between points (i) and (i+1)
        cost = (((ordered[:-1] - ordered[1:])**2).sum(1)).sum()
        #print(cost)
        if cost < mindist:
            mindist = cost
            minidx = i

    opt_order = paths[minidx]

    xxx = X[opt_order]
    yyy = Y[opt_order]

    # plt.figure()
    # plt.plot(xxx,yyy)
    # plt.show()

    # fits a spline to the ordered coordinates
    tckp, u = interpolate.splprep([xxx, yyy], s = s, k = k) #s = smoothing is 100, nest=-1
    x_new, y_new = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp)
    dx, dy = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp,der=1)

    #calculating normalized derivative points
    magnitude=np.sqrt(dx**2+dy**2)
    dx=dx/magnitude
    dy=dy/magnitude

    if verbose:
        stop = toc(start)

    return x_new, y_new, dx, dy

def straighten_image3D(image_to_be_straighten,X,Y,dx,dy,sampling_space=1, verbose=False):
    """
    streighten the image based on the spline and the derivatives of the spline.
    Goes over the z stack, and applies this 2D spline on the whole image.

    Args:
        image_to_be_straighten ([type]): [description]
        X ([type]): [description]
        Y ([type]): [description]
        dx ([type]): [description]
        dy ([type]): [description]
        verbose (bool, optional): [description]. Defaults to False.
    
    Returns:
        straightened_image: 3D straightened image

    """
    width_worm=100  #find optimal worm_with
    
    shapex=image_to_be_straighten.shape[2]
    shapey=image_to_be_straighten.shape[1]
    shapez=image_to_be_straighten.shape[0]

    #create new coordinate system
    # new coord system= old coordinate system + (-dy,dx)*new coordinate system
    n = np.arange(-width_worm,width_worm,sampling_space)
    xcoord = X[:,None] - n[None,:]*dy[:,None]
    ycoord = Y[:,None] + n[None,:]*dx[:,None]
    
    grid_x, grid_y = np.meshgrid(np.arange(0,shapex),np.arange(0,shapey))
    coord_old = np.array([grid_x.reshape(shapex*shapey),grid_y.reshape(shapex*shapey)]).T

    #create empty image
    straightened_image=np.empty([shapez,xcoord.shape[0],xcoord.shape[1]])

    #go over each slide and fill the straightened image
    for zslide in range(shapez):
        print(zslide)
        intensity_old = image_to_be_straighten[zslide,:,:].reshape(shapex*shapey)
        intensity_new = interpolate.griddata(coord_old, intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)))
        straightened_image[zslide,:,:]= np.squeeze(intensity_new).reshape(xcoord.shape)

    if (verbose==True):
        plt.figure()
        plt.imshow(np.max(image_to_be_straighten,0).T)
        plt.plot(X,Y)
        plt.plot(xcoord.T,ycoord.T)
        plt.title("original image")
        plt.show()

        plt.figure()
        plt.imshow(np.max(straightened_image,0).T)
        plt.title("straightened image")
        plt.show()

    return straightened_image

def straighten_image2D(image_to_be_straighten,X,Y,dx,dy,sampling_space=1,width_worm=150, verbose=False):
    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Straightening the worm. Verbose mode.')

    # create new coordinate system
    # new coord system= old coordinate system + (-dy,dx)*orthogonal coordinate
    n = np.arange(-width_worm,width_worm,sampling_space)
    xcoord = X[:,None] - n[None,:]*dy[:,None]
    ycoord = Y[:,None] + n[None,:]*dx[:,None]

    shapex=image_to_be_straighten.shape[1]
    shapey=image_to_be_straighten.shape[0]
    
    grid_x, grid_y = np.meshgrid(np.arange(0,shapex),np.arange(0,shapey))
    coord_old = np.array([grid_x.reshape(shapex*shapey),grid_y.reshape(shapex*shapey)]).T

    #new image
    intensity_old = image_to_be_straighten.reshape(shapex*shapey)
    intensity_new = interpolate.griddata(coord_old, intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)))
    straightened_image= np.squeeze(intensity_new).reshape(xcoord.shape)
    straightened_image=np.nan_to_num(straightened_image)

    if (verbose==True):
        stop = toc(start)    

    return straightened_image, (xcoord, ycoord)

def straighten_image2D_dual(images_to_be_straightened,X,Y,dx,dy,sampling_space=1,width_worm=150, verbose=False):
    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Straightening the worm. Verbose mode.')


    straightened_images = []

    # create new coordinate system
    # new coord system= old coordinate system + (-dy,dx)*orthogonal coordinate
    n = np.arange(-width_worm,width_worm,sampling_space)
    xcoord = X[:,None] - n[None,:]*dy[:,None]
    ycoord = Y[:,None] + n[None,:]*dx[:,None]

    shapex=np.shape(images_to_be_straightened)[2]
    shapey=np.shape(images_to_be_straightened)[1]
    
    grid_x, grid_y = np.meshgrid(np.arange(0,shapex),np.arange(0,shapey))
    coord_old = np.array([grid_x.reshape(shapex*shapey),grid_y.reshape(shapex*shapey)]).T

    if verbose:
        print('Creating the grid.', end = " ")
        stop = toc(start)
        start = tic()

    for k, image_to_be_straighten in enumerate(images_to_be_straightened):
        #new image
        intensity_old = image_to_be_straighten.reshape(shapex*shapey)
        intensity_new = interpolate.griddata(coord_old, intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)))
        straightened_image= np.squeeze(intensity_new).reshape(xcoord.shape)
        straightened_image=np.nan_to_num(straightened_image)

        # Stacking
        straightened_images.append(straightened_image)

        if verbose:
            print('Straightening a worm.', end = " ")
            stop = toc(start)
            start = tic()

    return straightened_images, (xcoord, ycoord)

def head2tail_masking(X, Y, dx, dy, img_binary, cut_th=0.2, verbose=False):

    # Debugging and benchmarking
    if verbose:
        start = tic()
        print('Cutting the worm at '+str(cut_th*100)+'%. Verbose mode.')

    shapex = img_binary.shape[-2]
    shapey = img_binary.shape[-1]

    # Define the points
    grid_y, grid_x = np.meshgrid(np.arange(0,shapey),np.arange(0,shapex))

    # Find the first line
    points2mask = np.linspace(0,1,np.shape(X)[0])

    # Creating storing matrices
    npoints = 21
    adj_mat = np.zeros([npoints, shapex, shapey])
    k_mat = np.zeros([shapex, shapey])

    # Computing distances
    for k, cut_ths in enumerate(np.linspace(0,1,npoints)):
        point_ref = np.sum(points2mask<cut_ths)
        adj_mat[k,:,:] = np.sqrt( (X[point_ref]-grid_x)**2+(Y[point_ref]-grid_y)**2)

    # Computing the region of thresholds
    val_mat = np.where(adj_mat == np.min(adj_mat,axis = 0))
    k_mat[val_mat[1],val_mat[2]] = val_mat[0]/(npoints-1)

    
    if verbose:
        print('Computing the domains.', end = " ")
        stop = toc(start)
        start = tic()

    # Find the upperline
    upper = np.sum(points2mask<=cut_th)
    lower = np.sum(points2mask<=(1-cut_th))

    # Reference point
    ref_up = X[upper]*dx[upper]+Y[upper]*dy[upper]
    ref_low = X[lower]*dx[lower]+Y[lower]*dy[lower]

    # define the upper points to be zero
    fun_u = (grid_y*dy[upper]+grid_x*dx[upper])
    fun_l = (grid_y*dy[lower]+grid_x*dx[lower])
    binary_upper = fun_u<ref_up
    binary_lower = fun_l>ref_low

    # Computing the new thresholded areas
    final_upper = (1-((binary_upper).astype(int)*(k_mat<=cut_th).astype(int)))
    final_lower = (1-((binary_lower).astype(int)*(k_mat>=(1-cut_th)).astype(int)))

    #new image
    binary_grid = final_upper*final_lower

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
        print('Computing the derivative domains.', end = " ")
        stop = toc(start)
        
    return binary_new

def arc_length(x, y):
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)
    
    return arc