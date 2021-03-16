#%% importing
import pandas as pd
from marit_functions import image_lists_mcherry_GFP
from marit_functions import read_image
from marit_functions import img_thresholding
from marit_functions import select_zslides
from marit_functions import calculate_worm_properties
from marit_functions import get_meta_info_temp

import skimage.morphology as skimorph
from lucas_functions import generating_ball
from lucas_functions import adaptive_masking
from lucas_functions import mask_postprocessing
from lucas_functions import dataset_comparison
from lucas_functions import masking_summary
# import lucas_functions


#%% load parameters
dirpath = 'C:/Users/moraluca/Desktop/Lin28_test'
#dirpath = 'C:\Users\moraluca\Desktop\Lin28_test'
# directory = 'W:/tungsten/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210116/'

outputpath = 'C:/Users/moraluca/Desktop/Lin28_test/Output'
channel_GFP =   'w1Lucas-sim-488-561'
channel_mcherry= 'w2Lucas-sim-561-488'

#save retults
results = []
image=1

#%% to run

# list for all channels the stk files in folder
list_mcherry, list_GFP = image_lists_mcherry_GFP(dirpath, channel_mcherry, channel_GFP)
print(list_mcherry)

#%%
#open mcherry and segment on signal
for i,(file1, file2) in enumerate(zip(list_mcherry, list_GFP)):
    print (i)
    if (i==image):
        #imageprint(file1)
        img_mcherry = read_image(file1)
        img_gfp = read_image(file2)
        print(img_mcherry.shape)

        # Parameters
        # Min/Max threshold
        mm_th = 1.8
        # Threshold_selection
        th_sel = .5

        # Running the function
        THPX, SORTED, ADPT, PRETH = adaptive_masking(img_mcherry, mm_th, th_sel)

        # Mask post-processing
        # Kernel
        krn = skimorph.disk(1)
        THPX = mask_postprocessing(THPX, krn)

        # Presenting outputs
        masking_summary(PRETH, SORTED, ADPT, mm_th)

        # Plott result
        plt.figure()
        plt.imshow(THPX[10,:,:])

        plt.figure()
        plt.plot(np.sum(np.sum(THPX,axis=2),axis = 1))


        # #preprocessing and thresholding to make image binary
        # img_binary = img_thresholding(img_mcherry)

        # #select slides based on mean signal
        # img_signal, img_binary = select_zslides(img_gfp, img_binary)

        # #calculates properties of the segmented worm
        # binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
        # img_binary, img_signal)

        # #create overlay of binary image with GFP image
        # img_overlay = img_signal * img_binary

        # #add properties in current results
        # current_res = get_meta_info_temp(file2)  #get metadata
        # current_res['volume'] = area  #calculate volume
        # current_res['mean_intensity'] = mean_intensity
        # current_res['min_intensity'] = min_intensity
        # current_res[final_intensity'] = mean_intensity - min_intensity  #calculate intensity

        # #save in resultarray
        # results.append(current_res)

        break
 # %%
