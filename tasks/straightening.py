#%% importing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from skimage import data
import matplotlib.pyplot as plt

#from skimage.util import invert
#import cv2

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

from skimage.morphology import skeletonize, skeletonize_3d

def create_skeleton(twoDimage,twoDmask):
  from fil_finder import FilFinder2D
  import astropy.units as u

  fil = FilFinder2D(twoDimage, mask=twoDmask) #creates a class, add beamwith?!
  # idea, put the fill finder on one plane where the worm is in focus 
  fil.preprocess_image()
  fil.create_mask(use_existing_mask=True)
  fil.medskel(verbose=False)
  fil.analyze_skeletons(branch_thresh=40* u.pix, skel_thresh=10 * u.pix, prune_criteria='length')
  plt.figure()
  plt.imshow(np.max(img_mcherry,0), cmap='gray')
  plt.contour(fil.skeleton_longpath, colors='r')
  plt.axis('off')
  plt.show()

  skeletonized_image=fil.skeleton_longpath
  [input_x,input_y]=np.where(skeletonized_image==1)
  outputarray = np.c_[input_x, input_y]
  return input_x, input_y




#%% load parameters
dirpath = '/Users/Marit/Documents/work/test_images/HBL1gfp_worm6'
outputpath = 'Users/Marit/Documents'
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#save retults
results = []
image=1

#%% to run
# list for all channels the stk files in folder
list_mcherry, list_GFP = image_lists_mcherry_GFP(dirpath, channel_mcherry, channel_GFP)

#open mcherry and segment on signal
for i,(file1, file2) in enumerate(zip(list_mcherry, list_GFP)):
    if (i==15):
        print(file1)

        img_mcherry = read_image(file1)
        img_gfp = read_image(file2)

        #preprocessing and thresholding to make image binary
        img_binary = img_thresholding(img_mcherry)

        #creating skeleton
        twoDimage=np.max(img_mcherry,0)
        twoDmask=np.max(img_binary,0)
        X,Y=create_skeleton(twoDimage,twoDmask)

       # %%
        tckp, u = interpolate.splprep([X, Y], s=3, k=2, nest=-1)
        xpointsnew, ypointsnew = interpolate.splev(np.linspace(0,1,270), tckp)

        plt.figure(3)
        plt.plot(xpointsnew, ypointsnew, 'r-')
        plt.show()


       # %%
      #runs a nearest neighbors algorithm on the coordinate array or spline interpolation to order the coordinates

      

        from sklearn.neighbors import KNeighborsRegressor
        neigh = KNeighborsRegressor(n_neighbors=2)
        G = neigh.kneighbors_graph()
        #neigh.fit(X,Y)
        #print(neigh.predict([[1.5]]))
      
        
  
        # %%

        clf = NearestNeighbors(2).fit(inputarray)
        
       
        T = nx.from_scipy_sparse_matrix(G)

        # sorts coordinates according to their nearest neighbors order
        order = list(nx.dfs_preorder_nodes(T, 0))
        xx = input_x[order]
        yy = input_y[order]

        # Loops over all points in the coordinate array as origin, determining which results in the shortest path
        paths = [list(nx.dfs_preorder_nodes(T, i)) for i in range(len(inputarray))]

        mindist = np.inf
        minidx = 0

        for i in range(len(inputarray)):
            p = paths[i]           # order of nodes
            ordered = inputarray[p]    # ordered nodes
            # find cost of that order by the sum of euclidean distances between points (i) and (i+1)
            cost = (((ordered[:-1] - ordered[1:])**2).sum(1)).sum()
            if cost < mindist:
                mindist = cost
                minidx = i

        opt_order = paths[minidx]

        xxx = input_x[opt_order]
        yyy = input_y[opt_order]

        # fits a spline to the ordered coordinates
        tckp, u = interpolate.splprep([xxx, yyy], s=3, k=2, nest=-1)
        xpointsnew, ypointsnew = interpolate.splev(np.linspace(0,1,270), tckp)

        # prints spline variables
        print(tckp)

        # plots the spline
        plt.figure(3)
        plt.plot(xpointsnew, ypointsnew, 'r-')
        plt.show()




   
          # %% 
          # rest of the code


        #select slides based on mean signal
        img_signal, img_binary = select_zslides(img_gfp, img_binary)
    

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


# %%
input_x = [-2.5, 0.0, 2.5, 5.0, 7.5, 3.0, -1.0]
input_y = [0.7, -6, 5, 6.5, 0.0, 5.0, -2.0]

x, y, yaw, k, travel = spline.calc_2d_spline_interpolation(input_x, input_y, num=200)

plt.subplots(1)
plt.plot(input_x, input_y, "xb", label="input")
plt.plot(x, y, "-r", label="spline")
plt.grid(True)
plt.axis("equal")
plt.xlabel("x[m]")
plt.ylabel("y[m]")
plt.legend()
# %%


