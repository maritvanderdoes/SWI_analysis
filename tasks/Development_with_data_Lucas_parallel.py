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
from utils.benchmarking import tic, toc
# straightening
from utils import create_spline,straighten_image3D, straighten_image2D, head2tail_masking
from utils import straighten_image2D_dual

from utils import arc_length
from utils.deprecated import create_skeleton
#from utils.core_utils import create_skeleton

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import os

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
# Other names will just run the code completely
# steps = 'Cropping'
steps = 'Complete'

# Other coding parameters
dwnscl = False # Downscaling the images
svngtm = False # Saving a temporal file
sngpln = True  # Checking a single plane (similar to ben's code)

# Setting the parallel mode
parallel = True
n_workers = 7

# Parameters in change
sorting = False
mm_th = 2.2 #2.5
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
image = 0

# list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
#print(list_mcherry)

#%% define the main function
def main_function(list_mcherry,list_GFP):
    # For debugging
    start0 = tic()

    # Determining the files to read
    files = (list_mcherry, list_GFP)

    # Reading the image and metadata
    print('MAIN: Reading image. File selected :'+files[0])
    (img_mcherry, img_gfp), current_res = \
        read_image_and_metadata(files, data_format = data_format)

    if debugging:
        # Creating the folder for the worm
        if sngpln:
            foldername = outputpath+'/Single_Plane_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
        else:
            foldername = outputpath+'/Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
        # Creating the individual debugging folder
        if not os.path.exists(foldername):
            print('MAIN: Creating debugging folder.')
            os.makedirs(foldername)

    status = 'File read'

    # Downscaling (for testing purposes)
    if dwnscl:
        (img_mcherry, img_gfp) = downscaling((img_mcherry, img_gfp), verbose = verbosity)

    try:
        # This loop does nothing, just to allow to break the function.
        # A better practice is to return every single time.
        for it in range(0,1):
            # Create mask
            print('MAIN: Adaptive Masking')
            img_binary, (sorted_values, pixel_threshold, pixel_range, area_zplane)  = \
                adaptive_masking(img_mcherry, mm_th = mm_th, th_sel = th_sel, krn_size = krn_size, 
                exp_size = exp_size, fill_holes = fill_holes, z_threshold = z_threshold, 
                sorting = sorting, verbose = verbosity)

            # Creating overlay of binary image with GFP image
            img_overlay = img_gfp * img_binary
            
            status = 'Masked Image'
            if steps == 'Masking':
                break

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
            if steps == 'Cropping':
                break

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
            if steps == 'Skeletonisation':
                break

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
            if steps == 'Straightening':
                break

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
            if steps == 'Cutting_Straightened':
                break
        
        #length of worm
        print('MAIN: Calculating length.')
        length = length/xdim
        bens_volume=np.pi*current_res['area_straightened']**2/(4*length)

        current_res['bens_volume'] = bens_volume

        status = 'Completed'

    except:
        print('MAIN: Some computations have not finished.')
        current_res['mean_curved'] = np.nan  #calculate volume
        current_res['volume_curved'] = np.nan 
        current_res['mean_curved_ht'] = np.nan 
        current_res['volume_curved_ht'] = np.nan 
        current_res['mean_straightened'] = np.nan 
        current_res['area_straightened'] = np.nan 
        current_res['mean_straightened_ht'] = np.nan 
        current_res['area_straightened_ht'] = np.nan 
        current_res['bens_volume'] = np.nan 

    stop_alg = toc(start0)

    # Saving temporal files
    if debugging:
        try:
            # Debugging and benchmarking
            start = tic()
            print('MAIN: Saving outputs. Verbose mode.', end = " ")

            # Computing better metrics
            vmin_val = np.min(cropped_image[cropped_binary])*.95
            img_dim = xcoord.shape
            img_dim2 = max_image.shape
            dwncoord = np.linspace(0,img_dim[0]-1,20).astype(int)

            # Saving status
            with open(foldername+'/A00_Status.txt', 'w') as f:
                f.write(status + ', with time of '+ str(np.round(stop_alg,2)) + ' (seconds)')

            # Presenting masking outputs
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane,
                mm_th = mm_th, scale = 'log', foldername = foldername)

            # Saving Mask
            imsave(foldername+'/A10_Mask'+'.tiff',np.float32(255*img_binary), check_contrast = False)
            # Saving Masked Data
            imsave(foldername+'/A11_Masked_data'+'.tiff',np.float32(img_overlay), check_contrast = False)
            
            # Saving Cropped Mask
            imsave(foldername+'/A20_Cropped_mask'+'.tiff',np.float32(cropped_binary), check_contrast = False)
            # Saving Cropped Data
            imsave(foldername+'/A21_Cropped_data'+'.tiff',np.float32(cropped_image), check_contrast = False)
            
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
            plt.savefig(foldername+'/A30_Pre_Skeleton.png')
            plt.close(fig)

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
            plt.savefig(foldername+'/A31_Post_Skeleton.png')
            plt.close(fig)

            # Saving the cut mask
            fig = plt.figure()
            if img_dim[0] > img_dim[1]:
                plt.imshow(np.max(cropped_binary_ht,0),cmap='gray')
            else:
                plt.imshow(np.max(cropped_binary_ht,0).T,cmap='gray')
            plt.title('Reduced length by '+str(0.2*100)+"%")
            plt.axis('off')
            plt.savefig(foldername+'/A40_Reduced_worm.png')
            plt.close(fig)

            fig, axs = plt.subplots(1,2)   
            if img_dim2[0] < img_dim2[1]:
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
            fig.savefig(foldername+'/A50_Straightened_mask.png')
            plt.close(fig)

            fig, axs = plt.subplots(1,2)   
            if img_dim2[0] < img_dim2[1]:
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
            fig.savefig(foldername+'/A51_Straightened_image.png')
            plt.close(fig)

            fig, axs = plt.subplots(2,1)   
            axs[0].plot(np.sum(straightened_image,axis = 1))
            axs[0].set_title("Projection over length")
            axs[0].set_ylabel('Intensity')
            axs[1].plot(np.sum(straightened_image,axis = 1)/np.sum(straightened_binary, axis = 1))
            axs[1].set_title("Mean projection over length")
            axs[1].set_xlabel('Length')
            axs[1].set_ylabel('Mean Intensity')
            fig.savefig(foldername+'/A52_Projection_over_length.png')
            plt.close(fig)

            stop = toc(start)

        except:
            print('MAIN: Some graphs have not been saved.')

            # Saving status
            with open(foldername+'/A00_Status.txt', 'a') as f:
                f.write('\n')
                f.write('There are some issues saving.')

        
    if svngtm:
        print('MAIN: Saving temp file.')
        current_res_df.append(current_res)
        df = pd.DataFrame(current_res_df)
        df.to_csv(outputpath+'/Results_temp_s'+(current_res['Position'])+'_t'+(current_res['Frame'])+'.csv', index=False)

    return current_res


#%%
if parallel:
    print('MAIN: Running in Parallel mode. Make sure you are running this in the terminal.')
    start_parallel = tic()
    if __name__=='__main__':
        p = multiprocessing.Pool(n_workers)
        results = p.starmap(main_function,zip(list_mcherry,list_GFP))
        p.close()
    print('Parallel code finished. ', end = ' ')
    stop = toc(start_parallel)

else:
    print('MAIN: Running in Sequential mode.')
    for k,files in enumerate(zip(list_mcherry, list_GFP)):
        print('Sample selected: '+str(k))
        if k == image :
            current_res = main_function(files[0],files[1])
            results.append(current_res)
            break
        if (not image) and (k != image):
            current_res = main_function(files[0],files[1])
            results.append(current_res)


#%% Saving results
df = pd.DataFrame(results)
if sngpln:
    df.to_csv(outputpath+'/Results_parallel_single_plane.csv', index = False)
else:
    df.to_csv(outputpath+'/Results_parallel.csv', index = False)

#%% To run the code
# Additional notes
# import sklearn gave an error. To solve type pip install -U scikit-learn
#
#
# python tasks/Development_with_data_Lucas_parallel.py

