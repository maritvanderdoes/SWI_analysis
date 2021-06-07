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
from _parameters import data_format, verbosity


#totest
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops

#save retults
results = []
current_res_df = []
image = 0

#list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)

xdim=0.533
ydim=0.533
zdim=2

#%% define the main function
def main_function(list_mcherry,list_GFP):
    files = (list_mcherry,list_GFP)

    # Reading the image and metadata
    print('Reading image.')
    (img_mcherry, img_gfp), current_res = read_image_and_metadata(files)

    try:
        # create mask
        img_binary, crap = adaptive_masking(img_mcherry)

        # crop image for further processing
        cropped_binary,cropped_image = crop_image(img_binary, img_gfp)

        #calculating properties curved worms
        (mean_curved, volume_curved)=calculate_worm_properties(cropped_binary, cropped_image)

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
        current_res['mean_curved'] = np.nan  #calculate volume
        current_res['volume_curved'] = np.nan 
        current_res['mean_curved_ht'] = np.nan 
        current_res['volume_curved_ht'] = np.nan 
        current_res['mean_straightened'] = np.nan 
        current_res['area_straightened'] = np.nan 
        current_res['mean_straightened_ht'] = np.nan 
        current_res['area_straightened_ht'] = np.nan 
        current_res['bens_volume'] = np.nan 

    print('Saving temp file.')
    current_res_df.append(current_res)
    df = pd.DataFrame(current_res_df)
    df.to_csv(outputpath+'/Results_temp_s'+(current_res['Position'])+'_t'+(current_res['Frame'])+'.csv', index=False)

    return current_res


#%%
if __name__=='__main__':
    p = multiprocessing.Pool(12)
    results = p.starmap(main_function,zip(list_mcherry,list_GFP))
    p.close()

#%% Saving results
df = pd.DataFrame(results)
df.to_csv(outputpath+'/Results_parallel.csv', index=False)
