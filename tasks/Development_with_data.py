#%% Setting up
# Loading utils package
import _setup
# Loading tools
from utils import image_lists, read_image_and_metadata
# Computing tools
from utils import adaptive_masking, calculate_worm_properties
# Plotting tools
from utils.plotting import dataset_comparison, masking_summary, plotzslides
# Benchmarking
from utils import downscaling

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt

# Saving the mask
from skimage.io import imsave
import warnings

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity

# Other coding parameters
dwnscl = False
plttng = True

# Parameters in change
mm_th = 2.5
th_sel = 0.3
krn_size = 2
exp_size = 1 # a 19 seem to be able to bridge, but slows down the 
             # code considreably. A faster and better implementation
             # is to reduce the z_threshold.
z_threshold = 0.6

# Sorting mm_th = 1.8, th_sel = 0.3 and z_threshold = 0.3
# Not Sorting mm_th = 2.5, th_sel = 0.3 and z_threshold = 0.6

#save retults
results = []
image = 5

#%% list for all channels the stk files in folder

(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%% open mcherry and segment on signal

for k,files in enumerate(zip(list_mcherry, list_GFP)):
    print('Sample selected: '+str(k))

    if True:#(k==image):
        print('File selected :'+files[0])
        # Reading the image and metadata
        images_out, current_res = \
            read_image_and_metadata(files, data_format = data_format)
        print(current_res)

        # Downscaling (for testing purposes)
        if dwnscl:
            images_out = \
                downscaling(images_out, verbose = verbosity)

        # Running the masking
        binary_mask, (sorted_values, pixel_threshold, pixel_range, area_zplane) = \
            adaptive_masking(images_out[0], mm_th = mm_th, th_sel = th_sel, krn_size = krn_size, 
            exp_size = exp_size, z_threshold = z_threshold, sorting = False,
            verbose = verbosity)

        # Presenting outputs
        if plttng:
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane,
                mm_th = mm_th, scale = 'linear')
            plotzslides([10,11,15,18,19],images_out[0],binary_mask,images_out[1])

        # Calculating properties of the segmented worm
        binary_image, metrics = calculate_worm_properties(binary_mask, images_out[1])

        # Creating overlay of binary image with GFP image
        img_overlay = images_out[0] * binary_mask

        # Storing the properties in current results
        current_res.update(dict(zip(('volume','mean_intensity','min_intensity'), metrics)))
        current_res['final_intensity'] = metrics[1] - metrics[2]  #calculate intensity

        # Save results in array
        results.append(current_res)

        # Saving the mask
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            # Mask
            imsave(outputpath+'\Mask_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',255*binary_image)
            # Masked_Data
            imsave(outputpath+'\Masked_data_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',img_overlay)

        break
# %%

