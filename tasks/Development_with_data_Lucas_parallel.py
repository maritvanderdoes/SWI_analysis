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

from utils import arc_length
from utils.deprecated import create_skeleton
#from utils.core_utils import create_skeleton

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing


# Saving the mask
from skimage.io import imsave

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity, debugging
from _parameters import xdim, ydim, zdim


#totest
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops

# Other coding parameters
dwnscl = False
svngfl = True
svngtm = False

# Setting the running mode
parallel = False
n_workers = 5

# Parameters in change
sorting = True
mm_th = 2.5
th_sel = 0.3
krn_size = 1
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

#%% list for all channels the stk files in folder

(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%% define the main function
def main_function(list_mcherry,list_GFP):
    # Determining the files to read
    files = (list_mcherry, list_GFP)

    # Reading the image and metadata
    print('Reading image. File selected :'+files[0])
    (img_mcherry, img_gfp), current_res = \
        read_image_and_metadata(files, data_format = data_format)

    # Downscaling (for testing purposes)
    if dwnscl:
        (img_mcherry, img_gfp) = downscaling((img_mcherry, img_gfp), verbose = verbosity)

    try:
        # Create mask
        img_binary, (sorted_values, pixel_threshold, pixel_range, area_zplane)  = \
            adaptive_masking(img_mcherry, mm_th = mm_th, th_sel = th_sel, krn_size = krn_size, 
            exp_size = exp_size, fill_holes = fill_holes, z_threshold = z_threshold, 
            sorting = sorting, verbose = verbosity)

        # Presenting outputs
        if debugging:
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane,
                mm_th = mm_th, scale = 'linear')
            plotzslides([10,11,15,18,19],img_mcherry,img_binary,img_gfp)

        # crop image for further processing
        cropped_binary, cropped_image = crop_image(img_binary, img_gfp)

        #calculating properties curved worms
        (mean_curved, volume_curved) = calculate_worm_properties(cropped_binary, cropped_image)

        print('Skeletonisation.')
        #fit spline
        Xinput,Yinput=create_skeleton(cropped_image, cropped_binary)
        print('Spline.')
        X,Y,dx,dy=create_spline(Xinput,Yinput)

        #cutting off head and tail, calculating properties
        print('Head to tail.')
        cropped_binary_ht = head2tail_masking(X,Y,dx,dy,cropped_binary,cut_th=0.2)
        (mean_curved_ht, volume_curved_ht)=calculate_worm_properties(cropped_binary_ht, cropped_image)

        #maximum intensity projection
        max_binary = np.max(cropped_binary,0)
        max_image = np.max(cropped_image,0)

        print('Straightening.')
        #straighten worm
        straightened_image=straighten_image2D(max_image,X,Y,dx,dy)
        straightened_binary=straighten_image2D(max_binary,X,Y,dx,dy)

        #calculatinte intensity straightened worm
        (mean_straightened, volume_straightened)=calculate_worm_properties(straightened_binary,straightened_image)

        #cutting off head and tail, calculating properties
        cutoff=int(straightened_image.shape[0]*0.2)
        straightened_binary_ht = straightened_binary[cutoff:-cutoff,:]
        straightened_image_ht = straightened_image[cutoff:-cutoff,:]
        (mean_straightened_ht, volume_traightened_ht)=calculate_worm_properties(straightened_binary_ht, straightened_image_ht)
   
    
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
        length=arc_length(X,Y)/xdim
        bens_volume=np.pi*current_res['area_straightened']**2/(4*length)

        current_res['bens_volume']=bens_volume

    except:
        print('Error.')
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
    if svngfl:
        # Mask
        imsave(outputpath+'\Mask_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',255*binary_image, check_contrast = False)
        # Masked_Data
        imsave(outputpath+'\Masked_data_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',np.float32(img_overlay), check_contrast = False)

    if svngtm:
        print('Saving temp file.')
        current_res_df.append(current_res)
        df = pd.DataFrame(current_res_df)
        df.to_csv(outputpath+'/Results_temp_s'+(current_res['Position'])+'_t'+(current_res['Frame'])+'.csv', index=False)

    return current_res


#%%
if parallel:
    for k,files in enumerate(zip(list_mcherry, list_GFP)):
        print('Sample selected: '+str(k))
        current_res = main_function(files[0],files[1])
        results.append(current_res)
else:
    if __name__=='__main__':
        p = multiprocessing.Pool(n_workers)
        results = p.starmap(main_function,zip(list_mcherry,list_GFP))
        p.close()

#%% Saving results
df = pd.DataFrame(results)
df.to_csv(outputpath+'/Results_parallel.csv', index=False)

# To run the code
# python tasks/Development_with_data_Lucas_parallel.py