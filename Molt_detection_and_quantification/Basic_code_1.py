# %%
# FUll code
from sklearn import metrics
from utils import image_lists, single_image_lists, read_image_and_metadata
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import imageio

from skimage.measure import regionprops
import skimage.morphology as morph 
from scipy.ndimage.measurements import label 


#%%
segmentpath     = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/woelmich/20211203_blmp10mu_nhr25GFP/20211203_segmentation'
channel_BF      = 'w2smita-Brightfield-GFP'
data_format     = 'st'
file_format     = 'tif'

position_sel = 5 # 1-26; 3,5
num_th = 64

print('Running the mask processing')
# Index the segmentation list
# list_mask=sorted(glob.glob(os.path.join(segmentpath, "*.tif")))
list_mask=sorted(glob.glob(os.path.join(segmentpath,'*_s'+str(position_sel)+ '_*.tif')))
len(list_mask)

#%% Timepoint ordering
timepoint_vec = np.zeros(len(list_mask))
for k_f in np.arange(0,len(list_mask)):
    timepoint_tmp = list_mask[k_f].split('_t')[-1]
    timepoint_vec[k_f] = int(timepoint_tmp.split('.tif')[0])

new_list_mask = [list_mask[k] for k in timepoint_vec.argsort()]

# %% Opening image
pos_sel = 210 # 150, 100
file_sel = new_list_mask[pos_sel]
img_binary = np.array(imageio.imread(file_sel))>num_th

plt.figure(figsize = (12,12))
plt.imshow(img_binary)

# %%
ccs, num_ccs = label(img_binary) #set labels in binary image
properties=regionprops(ccs,img_binary,['area',"mean_intensity",'centroid','eccentricity','perimeter']) #calculates the properties of the different areas
best_region = max(properties, key=lambda region: region.area) #selects the biggest region

[y0,x0] = best_region.centroid
axis_length = 100
x1 = best_region.centroid[1]+axis_length*np.sin(best_region.orientation)
x2 = best_region.centroid[1]-axis_length*np.cos(best_region.orientation)
y1 = best_region.centroid[0]+axis_length*np.cos(best_region.orientation)
y2 = best_region.centroid[0]+axis_length*np.sin(best_region.orientation)

plt.figure(figsize = (12,12))
plt.imshow(img_binary)
plt.plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
plt.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
plt.plot(best_region.centroid[1],best_region.centroid[0],'c.', markersize=15)

minr, minc, maxr, maxc = best_region.bbox
bx = (minc, maxc, maxc, minc, minc)
by = (minr, minr, maxr, maxr, minr)
plt.plot(bx, by, '-b', linewidth=2.5)

best_region.eccentricity

# %%

skeleton = morph.skeletonize(img_binary)
worm_length = np.sum(skeleton)
plt.figure(figsize = (12,12))
plt.imshow(skeleton)
plt.plot(best_region.centroid[1],best_region.centroid[0],'c.', markersize=15)
plt.plot(bx, by, '-b', linewidth=2.5)

# morph.binary_dilation(skeleton, morph.disk(skeleton.sum()/2))
# %%
[X,Y] = np.meshgrid(np.arange(0,img_binary.shape[1]),np.arange(0,img_binary.shape[0]))
distance_from_centroid = np.sqrt((Y-y0)**2+(X-x0)**2)
plt.figure(figsize = (12,12))
plt.imshow(distance_from_centroid)
plt.plot(best_region.centroid[1],best_region.centroid[0],'c.', markersize=15)

centre_mat = (distance_from_centroid.max()-distance_from_centroid)*skeleton
central_point = np.argwhere(centre_mat == centre_mat.max())[0]

plt.figure(figsize = (12,12))
plt.imshow(centre_mat)
plt.plot(best_region.centroid[1],best_region.centroid[0],'c.', markersize=15)
plt.plot(central_point[1],central_point[0],'rx', markersize=10, linewidth = 5)


# %%
plt.figure(figsize = (12,12))
plt.imshow(img_binary)
plt.plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
plt.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
plt.plot(best_region.centroid[1],best_region.centroid[0],'c.', markersize=15)

minr, minc, maxr, maxc = best_region.bbox
bx = (minc, maxc, maxc, minc, minc)
by = (minr, minr, maxr, maxr, minr)
plt.plot(bx, by, '-b', linewidth=2.5)
plt.plot(central_point[1],central_point[0],'rx', markersize=10)

# %%
projected_length = np.sqrt(((skeleton.sum(axis = 0)>0).sum())**2 + ((skeleton.sum(axis = 1)>0).sum())**2 )
straightness = projected_length/worm_length

# %% Final collection of metrics
# Fixed data
[X,Y] = np.meshgrid(np.arange(0,img_binary.shape[1]),np.arange(0,img_binary.shape[0]))

metrics_vector = np.zeros([len(new_list_mask), 12])
for pos_sel in np.arange(0,len(new_list_mask)):
    print(pos_sel)
    # Loading image
    file_sel = new_list_mask[pos_sel]
    img_binary = np.array(imageio.imread(file_sel))>num_th

    # Computing best region
    ccs, num_ccs = label(img_binary) #set labels in binary image
    if num_ccs > 0:
        properties=regionprops(ccs,img_binary,['area',"mean_intensity",'centroid','eccentricity','perimeter']) #calculates the properties of the different areas
        best_region = max(properties, key=lambda region: region.area) #selects the biggest region
        [y0,x0] = best_region.centroid

        # Computing line metrics
        skeleton = morph.skeletonize(img_binary)
        projected_length = np.sqrt(((skeleton.sum(axis = 0)>0).sum())**2 + ((skeleton.sum(axis = 1)>0).sum())**2 )

        # Distances
        distance_from_centroid = np.sqrt((Y-y0)**2+(X-x0)**2)
        centre_mat = (distance_from_centroid.max()-distance_from_centroid)*skeleton

        # Metrics to measure
        total_area = img_binary.sum()
        selected_area = best_region.area
        eccentricty = best_region.eccentricity
        worm_length = np.sum(skeleton)
        worm_width = worm_length/selected_area
        [yc,xc] = np.argwhere(centre_mat == centre_mat.max())[0]
        straightness = projected_length/worm_length

        metrics_vector[pos_sel, 0:10] = [total_area, selected_area, y0,x0, eccentricty,\
            worm_length,worm_width,yc,xc,straightness]

#%% extra metrics
metrics_vector[1:,10] =  (np.diff(metrics_vector[:,2:4],axis = 0)**2).sum(axis=1)
metrics_vector[1:,11] = (np.diff(metrics_vector[:,7:9],axis = 0)**2).sum(axis=1)

coord_sel = [2,3]
coord_sel = [7,8]
kern_size = 5
y_avg = np.convolve(metrics_vector[:,coord_sel[0]], np.ones(kern_size)/kern_size, 'same')
x_avg = np.convolve(metrics_vector[:,coord_sel[1]], np.ones(kern_size)/kern_size, 'same')

avg_dist = np.sqrt(np.diff(y_avg)**2+np.diff(x_avg)**2)

#%%
kern_size = 11
plt.plot(np.convolve(avg_dist,np.ones(kern_size)/kern_size))

#%%
plt.plot((metrics_vector[:,9]))

# %%
plt.plot(np.log2(metrics_vector[:,6]))
# %%
kern_size = 11
plt.plot(np.convolve(metrics_vector[:,9],np.ones(kern_size)/kern_size))
plt.grid()

# %%
kern_size = 5
plt.plot(np.convolve(metrics_vector[:,9],np.ones(kern_size)/kern_size))
plt.grid()
# %%
kern_size = 1
plt.plot(np.convolve(metrics_vector[:,9],np.ones(kern_size)/kern_size))
plt.plot(np.convolve(metrics_vector[:,10],np.ones(kern_size)/kern_size))
plt.grid()
# %%
plt.plot(metrics_vector[:,9],metrics_vector[:,10],'x')
plt.grid()
# %%
kern_size = 3

plt.plot(np.convolve(np.log2(metrics_vector[:,1]),np.ones(kern_size)/kern_size))
plt.grid()
# %%
kern_size = 11

plt.plot(np.convolve(np.log2(metrics_vector[:,10]),np.ones(kern_size)/kern_size))
plt.plot(np.convolve(np.log2(metrics_vector[:,11]),np.ones(kern_size)/kern_size))
plt.grid()
# %%
kern_size = 11

plt.plot(np.convolve(np.log2(metrics_vector[:,1]),np.ones(kern_size)/kern_size))
plt.plot(np.convolve(np.log2(metrics_vector[:,10]),np.ones(kern_size)/kern_size))
plt.plot(np.convolve(np.log2(metrics_vector[:,11]),np.ones(kern_size)/kern_size))
plt.grid()
# %%
plt.plot(np.convolve(np.log2(avg_dist),np.ones(kern_size)/kern_size))
plt.plot(np.convolve(np.log2(metrics_vector[:,10]),np.ones(kern_size)/kern_size))
plt.plot(np.convolve(np.log2(metrics_vector[:,11]),np.ones(kern_size)/kern_size))
plt.grid()

# %%
kern_size = 11

plt.plot(np.convolve(np.log2(metrics_vector[:,1]),np.ones(kern_size)/kern_size,'same'))
plt.plot(np.convolve(np.log2(avg_dist),np.ones(kern_size)/kern_size,'same'))
plt.grid()
# %%
