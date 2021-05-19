#%% Setting up
# Loading utils package
import _setup
# Loading tools
from utils import image_lists, read_image_and_metadata
# Computing tools
from utils import adaptive_masking, calculate_worm_properties, crop_image
# Plotting tools
from utils.plotting import dataset_comparison, masking_summary, plotzslides
# Benchmarking
from utils import downscaling
# straightening
from utils import create_spline,straighten_image3D, straighten_image2D, head2tail_masking
from utils import straighten_image2D_dual


from utils import arc_length
from utils.deprecated import create_skeleton
from utils.benchmarking import tic,toc
#from utils.core_utils import create_skeleton

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import os
import skimage.morphology as skimorph


# Saving the mask
from skimage.io import imsave

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity, debugging
from _parameters import xdim, ydim, zdim


#totest
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops

# Coding parameters
# 'Masking', 'Cropping', 'Skeletonisation', 'Straightening', 'Cutting_Straightened'
steps = 'Complete'

# Other coding parameters
dwnscl = False
svngtm = False
sngpln = True

# Setting the parallel mode
parallel = False
n_workers = 6

# Parameters in change
sorting = False
mm_th = 1.8 #2.5
th_sel = 0.3
krn_size = 2
fill_holes = True
exp_size = 3 # a 19 seem to be able to bridge, but slows down the 
             # code considreably. A faster and better implementation
             # is to reduce the z_threshold.
z_threshold = 0.6

# Sorting mm_th = 1.8, th_sel = 0.3 and z_threshold = 0.3
# Not Sorting mm_th = 2.5, th_sel = 0.3 and z_threshold = 0.6

#save retults
results = []
current_res_df = []
image = 4

# list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
#print(list_mcherry)

#%%
if parallel:
    print('MAIN: Running in Parallel mode. Make sure you are running this in the terminal.')
    if __name__=='__main__':
        p = multiprocessing.Pool(n_workers)
        results = p.starmap(main_function,zip(list_mcherry,list_GFP))
        p.close()
else:
    print('MAIN: Running in Sequential mode.')
    for k,files in enumerate(zip(list_mcherry, list_GFP)):
        print('Sample selected: '+str(k))
        if k == image :
            # For debugging
            start0 = tic()

            # Reading the image and metadata
            print('MAIN: Reading image. File selected :'+files[0])
            (img_mcherry, img_gfp), current_res = \
                read_image_and_metadata(files, data_format = data_format)

            if debugging:
                # Creating the folder for the worm
                if sngpln:
                    foldername = outputpath+'\Single_Plane_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
                else:
                    foldername = outputpath+'\Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
                # Creating the individual debugging folder
                if not os.path.exists(foldername):
                    print('MAIN: Creating debugging folder.')
                    os.makedirs(foldername)

            status = 'File read'

            # Downscaling (for testing purposes)
            if dwnscl:
                (img_mcherry, img_gfp) = downscaling((img_mcherry, img_gfp), verbose = verbosity)


            # Create mask
            print('MAIN: Adaptive Masking')
            img_binary, (sorted_values, pixel_threshold, pixel_range, area_zplane)  = \
                adaptive_masking(img_mcherry, mm_th = mm_th, th_sel = th_sel, krn_size = krn_size, 
                exp_size = exp_size, fill_holes = fill_holes, z_threshold = z_threshold, 
                sorting = sorting, verbose = verbosity)

            # Creating overlay of binary image with GFP image
            img_overlay = img_gfp * img_binary
            
            status = 'Masked Image'

            
            # krn = skimorph.disk(3)
            # img_binary_temp = img_binary
            # for z_plane in np.arange(0,img_binary.shape[0]):
            #     img_binary_temp[z_plane,:,:] = skimorph.binary_dilation(img_binary[z_plane,:,:], krn)
            #     img_binary_temp[z_plane,:,:] = skimorph.binary_erosion(img_binary[z_plane,:,:], krn)

            # Cropping image for further processing
            print('MAIN: Cropping data.')
            cropped_binary, cropped_image = \
                crop_image(img_binary, img_gfp, verbose = verbosity)

            status = 'Cropped Image'

            # Find the best plane (to study other metrics)
            if sngpln:
                area_cropped = np.sum(np.sum(cropped_binary,axis=2),axis = 1)
                best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
                # Selecting the frame
                cropped_binary = cropped_binary[best_plane,:,:]
                cropped_image = cropped_image[best_plane,:,:]
                # Reshaping
                cropped_binary = cropped_binary[None,:,:]
                cropped_image = cropped_image[None,:,:]

            # Calculating curved worms properties 
            (mean_curved, volume_curved) = \
                calculate_worm_properties(cropped_binary, cropped_image, verbose = verbosity)
            
            # Storing the properties in current results
            current_res['mean_curved'] = mean_curved  #calculate volume
            current_res['volume_curved'] = volume_curved/(xdim*ydim*zdim)
            
            status = 'Cropped Image (Comp)'

            # Skelotinisation
            print('MAIN: Skeletonisation and Spline.')
            Xinput, Yinput = create_skeleton(cropped_image, cropped_binary, verbose = verbosity)
            status = 'Skeleton created'
            
            # Spline fitting
            X, Y, dx, dy = create_spline(Xinput, Yinput, verbose = verbosity)
            status = 'Spline created'

            # Cutting off head and tail
            print('MAIN: Head to tail.')
            cropped_binary_ht = head2tail_masking(X,Y,dx,dy,cropped_binary,cut_th=0.2, verbose = verbosity)
            status = 'Cut worm'
            
            # Calculating properties
            (mean_curved_ht, volume_curved_ht) = calculate_worm_properties(cropped_binary_ht, cropped_image, verbose = verbosity)

            # Storing the properties in current results
            current_res['mean_curved_ht'] = mean_curved_ht
            current_res['volume_curved_ht'] = volume_curved_ht/(xdim*ydim*zdim)

            # Maximum intensity projection
            max_binary = np.max(cropped_binary,0)
            max_image = np.max(cropped_image,0)

            status = 'Cut worm'

            # Calculating the width of the worm as 1/10th approximately of the worm length
            length = arc_length(X,Y)
            #adp_width = int([np.min(150,length/10)])
            adp_width = int(length/10)

            # Straightening the worm
            print('MAIN: Straightening.')
            (straightened_image, straightened_binary), (xcoord, ycoord) = \
                straighten_image2D_dual((max_image, max_binary), X, Y, dx, dy, width_worm = adp_width,verbose = verbosity)
            status = 'Straightened worm'

            # Calculating intensity straightened worm
            (mean_straightened, volume_straightened) = calculate_worm_properties(straightened_binary,straightened_image, verbose = verbosity)
            
            # Storing the properties in current results
            current_res['mean_straightened'] = mean_straightened
            current_res['area_straightened'] = volume_straightened/(xdim*ydim)

            status = 'Straightened worm (Comp)'

            #cutting off head and tail, calculating properties
            print('MAIN: Cutting Straigthened.')
            cutoff=int(straightened_image.shape[0]*0.2)
            straightened_binary_ht = straightened_binary[cutoff:-cutoff,:]
            straightened_image_ht = straightened_image[cutoff:-cutoff,:]
            (mean_straightened_ht, volume_traightened_ht)=calculate_worm_properties(straightened_binary_ht, straightened_image_ht, verbose = verbosity)
            
            # Storing the properties in current results
            current_res['mean_straightened_ht'] = mean_straightened_ht
            current_res['area_straightened_ht'] = volume_traightened_ht

            status = 'Cut straightened'
        
            #length of worm
            print('MAIN: Calculating length.')
            length = arc_length(X,Y)/xdim
            bens_volume=np.pi*current_res['area_straightened']**2/(4*length)

            current_res['bens_volume'] = bens_volume

            status = 'Completed'

            # Debugging and benchmarking
            start = tic()
            print('MAIN: Saving outputs. Verbose mode.', end = " ")

            # Computing better metrics
            vmin_val = np.min(cropped_image[cropped_binary])*.95
            img_dim = xcoord.shape
            img_dim2 = max_image.shape
            dwncoord = np.linspace(0,img_dim[0]-1,20).astype(int)

            # Saving status
            #with open(foldername+'\A00_Status.txt', 'w') as f:
            #    f.write(status + ', with time of '+ str(np.round(stop_alg,2)) + ' (seconds)')

            # Presenting masking outputs
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane,
                mm_th = mm_th, scale = 'log')

            # Saving Mask
            #imsave(foldername+'\A10_Mask'+'.tiff',np.float32(255*img_binary), check_contrast = False)
            # Saving Masked Data
            #imsave(foldername+'\A11_Masked_data'+'.tiff',np.float32(img_overlay), check_contrast = False)
            
            # Saving Cropped Mask
            #imsave(foldername+'\A20_Cropped_mask'+'.tiff',np.float32(cropped_binary), check_contrast = False)
            # Saving Cropped Data
            #imsave(foldername+'\A21_Cropped_data'+'.tiff',np.float32(cropped_image), check_contrast = False)
            
            # Saving the skeleton
            fig = plt.figure()
            if img_dim[0] > img_dim[1]:
                plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:], cmap='gray')
                plt.plot(Yinput,Xinput,'r.')
            else:
                plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:].T, cmap='gray')
                plt.plot(Xinput,Yinput,'r.')
            plt.axis('off')
            plt.title('Old Skeleton')
            #plt.savefig(foldername+'\A30_Pre_Skeleton.png')
            #plt.close(fig)

            # Saving the new skeleton
            fig = plt.figure()
            if img_dim[0] > img_dim[1]:
                plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:], cmap='gray')
                plt.plot(Y, X, 'r.')
            else:
                plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:].T, cmap='gray')
                plt.plot(X, Y, 'r.')
            plt.axis('off')
            plt.title('New Skeleton')
            #plt.savefig(foldername+'\A31_Post_Skeleton.png')
            #plt.close(fig)

            # Saving the cut mask
            fig = plt.figure()
            if img_dim[0] > img_dim[1]:
                plt.imshow(np.max(cropped_binary_ht,0),cmap='gray')
            else:
                plt.imshow(np.max(cropped_binary_ht,0).T,cmap='gray')
            plt.title('Reduced length by '+str(0.2*100)+"%")
            plt.axis('off')
            #plt.savefig(foldername+'\A40_Reduced_worm.png')
            #plt.close(fig)

            fig, axs = plt.subplots(1,2)   
            if img_dim[0] < img_dim[1]:
                axs[0].imshow(max_binary.T, cmap='gray')
                axs[0].plot(xcoord[dwncoord,:].T,ycoord[dwncoord,:].T)
                axs[0].plot(X,Y,'r')
            else:
                axs[0].imshow(max_binary, cmap='gray')
                axs[0].plot(ycoord[dwncoord,:].T,xcoord[dwncoord,:].T)
                axs[0].plot(Y,X,'r')
            axs[0].set_title("Original mask with sampling")
            axs[1].imshow(straightened_binary, cmap='gray')
            axs[1].set_title("Straightened mask")
            #fig.savefig(foldername+'\A50_Straightened_mask.png')
            #plt.close(fig)

            fig, axs = plt.subplots(1,2)   
            if img_dim[0] < img_dim[1]:
                axs[0].imshow(max_image.T, cmap='gray')
                axs[0].plot(xcoord[dwncoord,:].T,ycoord[dwncoord,:].T)
                axs[0].plot(X,Y,'r')
            else:
                axs[0].imshow(max_image, cmap='gray')
                axs[0].plot(ycoord[dwncoord,:].T,xcoord[dwncoord,:].T)
                axs[0].plot(Y,X,'r')
            axs[0].set_title("Original image with sampling")
            axs[1].imshow(straightened_image, cmap='gray', vmin = vmin_val)
            axs[1].set_title("Straightened image")
            #fig.savefig(foldername+'\A51_Straightened_image.png')
            #plt.close(fig)

            stop = toc(start)


#%% Saving results
df = pd.DataFrame(results)
df.to_csv(outputpath+'/Results_parallel.csv', index = False)

#%% To run the code
# Additional notes
# import sklearn gave an error. To solve type pip install -U scikit-learn
#
#
# python tasks/Development_with_data_Lucas_parallel.py



#%%
asd = np.zeros([10,10])
asd2 = np.ones([10,10])

mats = (asd,asd2)

print(np.shape(mats)[1])

for k, mat in enumerate(mats):
    print(mat)
    print('potat')


#%% Testing parts of the code
# Skelotinisation
print('MAIN: Skeletonisation and Spline.')
Xinput, Yinput = create_skeleton(cropped_image, cropped_binary, verbose = verbosity)
status = 'Skeleton created'

# Spline fitting
X, Y, dx, dy = create_spline(Xinput, Yinput, verbose = verbosity)
status = 'Spline created'

#%%
plt.plot(Xinput,Yinput)
# %%
plt.plot(X,Y)
# %% Dissected creat_spline
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from scipy import interpolate


sampling_space=1
s=100
k=3
n_neighbors = 2
radius = 1

#variables
nsplinepoints=int(Xinput.shape[0]/sampling_space)
inputarray = np.c_[Xinput, Yinput]
#inputarray = np.c_[np.flipud(Xinput), np.flipud(Yinput)]

# runs a nearest neighbors algorithm on the coordinate array
clf = NearestNeighbors(n_neighbors=n_neighbors, radius = radius).fit(inputarray)
G = clf.kneighbors_graph()
T = nx.from_scipy_sparse_matrix(G)

# sorts coordinates according to their nearest neighbors order
order = list(nx.dfs_preorder_nodes(T, 0))
xx = Xinput[order]
yy = Yinput[order]

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

xxx = Xinput[opt_order]
yyy = Yinput[opt_order]

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

# %%
plt.plot(Xinput,Yinput,'r')
plt.plot(xx,yy)
#plt.xlim([10, 20])
#plt.ylim([160, 200])
# %%
plt.spy(G, markersize = 1)
# %%
nx.draw_circular(T)


# %% Dissected 
cropped_binary_ht_test = head2tail_masking(X,Y,dx,dy,cropped_binary,cut_th=0.5, verbose = verbosity)

# %% Dissected

cut_th = .1
# # check for sizes
# if np.ndim(img_binary)==2:
#     shapez = 1
#     shapex = img_binary.shape[0]
#     shapey = img_binary.shape[1]
# elif np.ndim(img_binary)==3:
#     shapez = img_binary.shape[0]
#     shapex = img_binary.shape[1]
#     shapey = img_binary.shape[2]

shapex = cropped_binary.shape[-2]
shapey = cropped_binary.shape[-1]

# if verbose:
#     print("Number of z_planes: "+number_z)

# binary_to_test
# Define the points
grid_y, grid_x = np.meshgrid(np.arange(0,shapey),np.arange(0,shapex))

# Find the first line
points2mask = np.linspace(0,1,np.shape(X)[0])

# Find the upperline
upper = np.sum(points2mask<=cut_th)
lower = np.sum(points2mask<=(1-cut_th))

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
cropped_binary_ht_test = binary_grid*cropped_binary
if np.ndim(cropped_binary_ht_test)==2:
    cropped_binary_ht_test = binary_grid*cropped_binary_ht_test
elif np.ndim(img_binary)==3:
    cropped_binary_ht_test = binary_grid[None,:,:]*cropped_binary_ht_test


plt.figure()
plt.imshow(cropped_binary.squeeze())
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')

plt.figure()
plt.imshow((binary_upper))
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')
plt.figure()
plt.imshow((binary_lower))
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')


plt.figure()
plt.contour((fun_u)*cropped_binary.squeeze())
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')
plt.figure()
plt.imshow((binary_upper)*cropped_binary.squeeze())
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')

plt.figure()
plt.contour((fun_l)*cropped_binary.squeeze())
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')
plt.figure()
plt.imshow((binary_lower)*cropped_binary.squeeze())
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')

# %%
fig = plt.figure()
if img_dim[0] > img_dim[1]:
    plt.imshow(np.max(cropped_binary_ht_test,0),cmap='gray')
else:
    plt.imshow(np.max(cropped_binary_ht_test,0).T,cmap='gray')
plt.title('Reduced length by '+str(0.2*100)+"%")
plt.axis('off')
#plt.savefig(foldername+'\A40_Reduced_worm.png')
#plt.close(fig)
# %% Alternative stages (not great)
npoints = 11
adj_mat = np.zeros([npoints, shapex, shapey])
com_mat = np.zeros([npoints, shapex, shapey])
for k, cut_th in enumerate(np.linspace(0,1,npoints)):
    upper = np.sum(points2mask<cut_th)
    ref_up = X[upper]*dx[upper]+Y[upper]*dy[upper]
    fun_u = (grid_y*dy[upper]+grid_x*dx[upper])
    adj_mat[k,:,:] = (fun_u-ref_up)**2
    com_mat[k,:,:] = k

val_mat = np.where(adj_mat == np.min(adj_mat,axis = 0))

# %% Alternative stages (not great)
npoints = 21
adj_mat = np.zeros([npoints, shapex, shapey])
k_mat = np.zeros([shapex, shapey])

grid_y, grid_x = np.meshgrid(np.arange(0,shapey),np.arange(0,shapex))

for k, cut_ths in enumerate(np.linspace(0,1,npoints)):
    point_ref = np.sum(points2mask<cut_ths)
    adj_mat[k,:,:] = np.sqrt( (X[point_ref]-grid_x)**2+(Y[point_ref]-grid_y)**2)

val_mat = np.where(adj_mat == np.min(adj_mat,axis = 0))
k_mat[val_mat[1],val_mat[2]] = val_mat[0]/(npoints-1)

upper = np.sum(points2mask<=cut_th)
ref_up = X[upper]*dx[upper]+Y[upper]*dy[upper]
fun_u = (grid_y*dy[upper]+grid_x*dx[upper])
binary_upper = fun_u<ref_up

final_upper = (1-((binary_upper).astype(int)*(k_mat<=cut_th).astype(int)))#*cropped_binary.squeeze()

lower = np.sum(points2mask<=(1-cut_th))
ref_low = X[lower]*dx[lower]+Y[lower]*dy[lower]
fun_l = (grid_y*dy[lower]+grid_x*dx[lower])
binary_lower = fun_l>ref_low

final_lower = (1-((binary_lower).astype(int)*(k_mat>=(1-cut_th)).astype(int)))#*cropped_binary.squeeze()

final_matrix = (final_upper)*(final_lower)*cropped_binary.squeeze()

#k_mat = np.reshape(val_mat,(shapex,shapey))
# %%
plt.imshow(binary_upper)
# %%
plt.imshow(k_mat*cropped_binary.squeeze()>=.2)
# %%
plt.imshow(adj_mat[10,:,:])
# %%
plt.figure()
#plt.imshow((binary_lower).astype(int)-(k_mat>0).astype(int))
#plt.imshow((binary_upper).astype(int)+(k_mat<=.2).astype(int))
plt.imshow(final_matrix)
plt.plot(Y[upper],X[upper],'ro')
plt.plot(Y[lower],X[lower],'rx')
# %%
