#%% importing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from skimage import data
import matplotlib.pyplot as plt

from sklearn.neighbors import NearestNeighbors
from skimage import io
from scipy import interpolate
import matplotlib.pyplot as plt
import networkx as nx

from marit_functions import image_lists_mcherry_GFP
from marit_functions import read_image
from marit_functions import img_thresholding
from marit_functions import select_zslides
from marit_functions import calculate_worm_properties
from marit_functions import get_meta_info_temp
from marit_functions import calculatesometing
from marit_functions import warp_images

from scipy.interpolate import griddata
from skimage.morphology import skeletonize, skeletonize_3d
from sklearn.neighbors import KNeighborsRegressor


def create_skeleton(twoDimage,twoDmask, verbose=False):
  #creates a skeleton from an image and a mask in 2D, and returs x and y coordinates for the skeletion. Coordinates are not ordered
  from fil_finder import FilFinder2D
  import astropy.units as u

  fil = FilFinder2D(twoDimage, mask=twoDmask) #creates a class, add beamwith?!
  # idea, put the fill finder on one plane where the worm is in focus 
  fil.preprocess_image()
  fil.create_mask(use_existing_mask=True)
  fil.medskel(verbose=False)
  fil.analyze_skeletons(branch_thresh=40* u.pix, skel_thresh=10 * u.pix, prune_criteria='length')

  skeletonized_image=fil.skeleton_longpath
  [input_x,input_y]=np.where(skeletonized_image==1)

  if (verbose==True):
    #plot
    plt.figure()
    plt.imshow(twoDimage, cmap='gray')
    plt.contour(skeletonized_image, colors='r')
    plt.axis('off')
    plt.show()

  
  return input_x, input_y

def arc_length(x, y):
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)

    return arc

def create_spline(X,Y, verbose=False):
  # creates a spline through X and Y coordiantes by a number of splinepoints. Gives as output the rescaled x and y, and the derivatives'"""
  #idea explained on
  #https://stackoverflow.com/questions/37742358/sorting-points-to-form-a-continuous-line
  
  #variables
  s=1
  nsplinepoints=100 # optimize number of splinepoints

  
  
  inputarray = np.c_[X, Y]

  # runs a nearest neighbors algorithm on the coordinate array
  clf = NearestNeighbors(2).fit(inputarray)
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

  # fits a spline to the ordered coordinates
  tckp, u = interpolate.splprep([xxx, yyy], s=s, k=3, nest=-1)
  xpointsnew, ypointsnew = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp)
  dx, dy = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp,der=1)

  #calculating normalized derivative points
  magnitude=np.sqrt(dx**2+dy**2)
  dx=dx/magnitude
  dy=dy/magnitude

  if (verbose==True):
    plt.plot(xpointsnew,ypointsnew)
    plt.title('fittet spline')

  return xpointsnew,ypointsnew, dx,dy

def straighten_image(image_to_be_straighten,X,Y,dx,dy, verbose=False):
        #image to be straightened can be 3D. variables X,Y, dx, dy  come from create spline . Returned  straightened image"""
        width_worm=40  #find optimal worm_with
        #create new coordinate system
        # new coord system= old coordinate system + (-dy,dx)*new coordinate system
        n = np.arange(-width_worm,width_worm,1)
        xcoord = X[:,None] - n[None,:]*dy[:,None]
        ycoord = Y[:,None] + n[None,:]*dx[:,None]
        
        grid_x, grid_y = np.meshgrid(np.arange(0,image_to_be_straighten.shape[1]),np.arange(0,image_to_be_straighten.shape[2]))
        coord_old = np.array([grid_x.reshape(image_to_be_straighten.shape[1]*image_to_be_straighten.shape[2]),grid_y.reshape(image_to_be_straighten.shape[1]*image_to_be_straighten.shape[2])]).T

        #create empty image
        straightened_image=np.empty([image_to_be_straighten.shape[0],xcoord.shape[0],xcoord.shape[1]])

        #go over each slide and fill the straightened image
        for zslide in range(image_to_be_straighten.shape[0]):
          intensity_old = image_to_be_straighten[zslide,:,:].reshape(image_to_be_straighten.shape[1]*image_to_be_straighten.shape[2])
          intensity_new = griddata(coord_old, intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)))
          straightened_image[zslide,:,:]= np.squeeze(intensity_new).reshape(xcoord.shape)

        if (verbose==True):
          plt.figure()
          plt.imshow(np.max(image_to_be_straighten,0).T)
          plt.plot(X,Y)
          plt.plot(xcoord.T,ycoord.T)
          plt.title("original image")
          plt.show()

          plt.figure()
          plt.imshow(np.max(straightened_image,0))
          plt.titel("straightened image")

        return straightened_image

  # plt.figure(figsize=(10,10))
  # plt.imshow(img.T)
  # plt.plot(xpointsnew, ypointsnew, 'r-')
  # plt.plot(xcoord, ycoord, 'b.')
  # plt.show()

  #plt.imshow(intensity_new)
  

#%% load parameters
dirpath = '/Users/Marit/Documents/work/test_images'
outputpath = 'Users/Marit/Documents'
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#save retults
results = []
image=1

# list for all channels the stk files in folder
list_mcherry, list_GFP = image_lists_mcherry_GFP(dirpath, channel_mcherry, channel_GFP)

#open mcherry and segment on signal
for i,(file1, file2) in enumerate(zip(list_mcherry, list_GFP)):
    if (i==0):
        print(file1)

        img_mcherry = read_image(file1)
        img_gfp = read_image(file2)

        #select slides based on mean signal
        img_signal, img_binary = select_zslides(img_gfp, img_binary)
    
        #preprocessing and thresholding to make image binary
        img_binary = img_thresholding(img_mcherry)

        #creating skeleton
        twoDimage=np.max(img_mcherry,0)
        twoDmask=np.max(img_binary,0)

        Xinput,Yinput=create_skeleton(twoDimage,twoDmask)
        X,Y,dx,dy=create_spline(Xinput,Yinput, verbose=False)
        imgtoshow=straighten_image(img_binary,X,Y,dx,dy)


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
        #results.append(current_res)
