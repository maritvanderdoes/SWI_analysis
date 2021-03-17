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

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity

# Other coding parameters
dwnscl = False
plttng = False

# Parameters in change
krn_size = 5
exp_size = 5

#save retults
results = []
image = 0

#%% list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%%

#open mcherry and segment on signal
for k,files in enumerate(zip(list_mcherry, list_GFP)):
    print ('Sample selected: '+str(k))
    if (k==image):
        # Reading the image and metadata
        (img_mcherry, img_gfp), meta_out = \
            read_image_and_metadata(files, data_format = data_format)

        # Downscaling (for testing purposes)
        if dwnscl:
            (img_mcherry, img_gfp) = \
                downscaling((img_mcherry, img_gfp), verbose = verbosity)

        # Running the masking
        binary_mask, sorted_values, pixel_threshold, pixel_range, area_zplane = \
            adaptive_masking(img_mcherry, krn_size = krn_size, exp_size = exp_size, 
            verbose = verbosity)

        # Presenting outputs
        if plttng:
            masking_summary(sorted_values, pixel_threshold, pixel_range, area_zplane)
            plotzslides('Comparison',[10,11,15,18,19],img_mcherry,binary_mask,img_gfp)

        # Calculating properties of the segmented worm
        binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
        img_binary = binary_mask, img_signal = img_gfp)

        # Creating overlay of binary image with GFP image
        img_overlay = img_gfp * binary_mask

        # Storing the properties in current results
        current_res = meta_out  #get metadata
        current_res['volume'] = area  #calculate volume
        current_res['mean_intensity'] = mean_intensity
        current_res['min_intensity'] = min_intensity
        current_res['final_intensity'] = mean_intensity - min_intensity  #calculate intensity

        # Save results in array
        results.append(current_res)

        break
 # %%

