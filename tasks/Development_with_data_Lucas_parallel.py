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
from utils.benchmarking import tic, toc
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
steps = 'Complete'

# Other coding parameters
dwnscl = False
svngtm = False

# Setting the parallel mode
parallel = True
n_workers = 6

# Parameters in change
sorting = False
mm_th = 3 #2.5
th_sel = 0.3
krn_size = 2
fill_holes = True
exp_size = 1 # a 19 seem to be able to bridge, but slows down the 
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
    # Determining the files to read
    files = (list_mcherry, list_GFP)

    # Reading the image and metadata
    print('MAIN: Reading image. File selected :'+files[0])
    (img_mcherry, img_gfp), current_res = \
        read_image_and_metadata(files, data_format = data_format)

    if debugging:
        # Creating the folder for the worm
        foldername = outputpath+'\Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']
        # Creating the individual debugging folder
        if not os.path.exists(foldername):
            print('MAIN: Creating debugging folder.')
            os.makedirs(foldername)

    status = 'File read'

    # Downscaling (for testing purposes)
    if dwnscl:
        (img_mcherry, img_gfp) = downscaling((img_mcherry, img_gfp), verbose = verbosity)

    try:
        for it in range(0,1):
            # Create mask
            print('MAIN: Adaptive Masking.')
            img_binary, (sorted_values, pixel_threshold, pixel_range, area_zplane)  = \
                adaptive_masking(img_mcherry, mm_th = mm_th, th_sel = th_sel, krn_size = krn_size, 
                exp_size = exp_size, fill_holes = fill_holes, z_threshold = z_threshold, 
                sorting = sorting, verbose = verbosity)
            status = 'Masked Image'

            # Creating overlay of binary image with GFP image
            img_overlay = img_gfp * img_binary
            
            if steps == 'Masking':
                break

            # Cropping image for further processing
            print('MAIN: Cropping data.')
            cropped_binary, cropped_image = \
                crop_image(img_binary, img_gfp, verbose = verbosity)
            status = 'Cropped Image'

            # Calculating curved worms properties 
            (mean_curved, volume_curved) = \
                calculate_worm_properties(cropped_binary, cropped_image, verbose = verbosity)
                
            if steps == 'Cropping':
                break

            # Skelotinisation and spline fitting
            print('MAIN: Skeletonisation and Spline.')
            Xinput, Yinput = create_skeleton(cropped_image, cropped_binary, verbose = verbosity)
            status = 'Skeleton created'
            X, Y, dx, dy = create_spline(Xinput, Yinput, verbose = verbosity)
            status = 'Spline created'

            # Cutting off head and tail, calculating properties
            print('MAIN: Head to tail.')
            cropped_binary_ht = head2tail_masking(X,Y,dx,dy,cropped_binary,cut_th=0.2, verbose = verbosity)
            (mean_curved_ht, volume_curved_ht) = calculate_worm_properties(cropped_binary_ht, cropped_image, verbose = verbosity)
            status = 'Cut worm'

            # Maximum intensity projection
            max_binary = np.max(cropped_binary,0)
            max_image = np.max(cropped_image,0)

            if steps == 'Skeletonisation':
                break

            # Straightening the worm
            print('MAIN: Straightening.')
            (straightened_image, straightened_binary), (xcoord, ycoord) = \
                straighten_image2D_dual((max_image, max_binary), X, Y, dx, dy, verbose = verbosity)
            status = 'Straightened worm'

            # Calculating intensity straightened worm
            (mean_straightened, volume_straightened) = calculate_worm_properties(straightened_binary,straightened_image, verbose = verbosity)

            if steps == 'Straightening':
                break

            #cutting off head and tail, calculating properties
            print('MAIN: Cutting Straigthened.')
            cutoff=int(straightened_image.shape[0]*0.2)
            straightened_binary_ht = straightened_binary[cutoff:-cutoff,:]
            straightened_image_ht = straightened_image[cutoff:-cutoff,:]
            (mean_straightened_ht, volume_traightened_ht)=calculate_worm_properties(straightened_binary_ht, straightened_image_ht, verbose = verbosity)
            status = 'Cut straightened.'
            
            if steps == 'Cut_Straightening':
                break
        
            # Storing the properties in current results
            current_res['mean_curved'] = mean_curved  #calculate volume
            current_res['volume_curved'] = volume_curved/(xdim*ydim*zdim)
            current_res['mean_curved_ht'] = mean_curved_ht
            current_res['volume_curved_ht'] = volume_curved_ht/(xdim*ydim*zdim)
            current_res['mean_straightened'] = mean_straightened
            current_res['area_straightened'] = volume_straightened/(xdim*ydim)
            current_res['mean_straightened_ht'] = mean_straightened_ht
            current_res['area_straightened_ht'] = volume_traightened_ht

            #length of worm
            print('MAIN: Calculating length.')
            length = arc_length(X,Y)/xdim
            bens_volume=np.pi*current_res['area_straightened']**2/(4*length)
            current_res['bens_volume']=bens_volume

            status = 'Completed.'

    except:
        print('MAIN: Error.')
        current_res['mean_curved'] = np.nan  #calculate volume
        current_res['volume_curved'] = np.nan 
        current_res['mean_curved_ht'] = np.nan 
        current_res['volume_curved_ht'] = np.nan 
        current_res['mean_straightened'] = np.nan 
        current_res['area_straightened'] = np.nan 
        current_res['mean_straightened_ht'] = np.nan 
        current_res['area_straightened_ht'] = np.nan 
        current_res['bens_volume'] = np.nan 

    # Saving temporal files
    if debugging:
        try:
            # Debugging and benchmarking
            start = tic()
            print('MAIN: Saving outputs. Verbose mode.', end = " ")

            # Saving status
            with open(foldername+'\A00_Status.txt', 'w') as f:
                f.write(status)

            # Presenting masking outputs
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane,
                mm_th = mm_th, scale = 'log', foldername = foldername)

            # Saving Mask
            imsave(foldername+'\A10_Mask'+'.tiff',np.float32(255*img_binary), check_contrast = False)
            # Saving Masked Data
            imsave(foldername+'\A11_Masked_data'+'.tiff',np.float32(img_overlay), check_contrast = False)
            
            # Saving Cropped Mask
            imsave(foldername+'\A20_Cropped_mask'+'.tiff',np.float32(cropped_binary), check_contrast = False)
            # Saving Cropped Data
            imsave(foldername+'\A21_Cropped_data'+'.tiff',np.float32(cropped_image), check_contrast = False)
            
            # Saving the skeleton
            plt.figure()
            plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:], cmap='gray')
            plt.plot(Yinput,Xinput,'r.')
            plt.axis('off')
            plt.title('Old Skeleton')
            plt.savefig(foldername+'\A30_Pre_Skeleton.png')

            # Saving the new skeleton
            plt.figure()
            plt.imshow(cropped_binary[int(cropped_binary.shape[0]/2),:,:], cmap='gray')
            plt.plot(Y, X, 'r.')
            plt.axis('off')
            plt.title('New Skeleton')
            plt.savefig(foldername+'\A31_Post_Skeleton.png')

            # Saving the cut mask
            plt.figure()
            plt.imshow(np.max(cropped_binary_ht,0),cmap='gray')
            plt.title('Reduced length by '+str(0.2*100)+"%")
            plt.axis('off')
            plt.savefig(foldername+'\A40_Reduced_worm.png')

            plt.figure()
            fig, axs = plt.subplots(1,2)   
            axs[0].imshow(max_binary.T)
            axs[0].plot(X,Y)
            axs[0].plot(xcoord.T,ycoord.T)
            axs[0].set_title("Original mask with sampling")
            axs[1].imshow(straightened_binary, cmap='gray')
            axs[1].set_title("Straightened mask")
            fig.savefig(foldername+'\A50_Straightened_mask.png')

            plt.figure()
            fig, axs = plt.subplots(1,2)   
            axs[0].imshow(max_image.T)
            axs[0].plot(X,Y)
            axs[0].plot(xcoord.T,ycoord.T)
            axs[0].set_title("Original image with sampling")
            axs[1].imshow(straightened_image, cmap='gray')
            axs[1].set_title("Straightened image")
            fig.savefig(foldername+'\A51_Straightened_image.png')

            stop = toc(start)

        except:
            print('MAIN: Saving Error.')
        
    if svngtm:
        print('MAIN: Saving temp file.')
        current_res_df.append(current_res)
        df = pd.DataFrame(current_res_df)
        df.to_csv(outputpath+'/Results_temp_s'+(current_res['Position'])+'_t'+(current_res['Frame'])+'.csv', index=False)

    return current_res


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
            current_res = main_function(files[0],files[1])
            results.append(current_res)
            break
        if (not image) and (k != image):
            current_res = main_function(files[0],files[1])
            results.append(current_res)


#%% Saving results
df = pd.DataFrame(results)
df.to_csv(outputpath+'/Results_parallel.csv', index = False)

#%% To run the code
# Additional notes
# import sklearn gave an error. To solve type pip install -U scikit-learn
#
#
# python tasks/Development_with_data_Lucas_parallel.py

