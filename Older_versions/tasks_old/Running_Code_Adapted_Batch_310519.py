#%% Setting up
# Notes on updates
# Functions that have been changed:
#   adaptive_masking, new default parameters updated
#   calculate_worm_properties, from 0 to nans
# New Functions:
#   straightend_image2D_dual
# Loading utils package
import _setup
# Loading tools
from utils import image_lists, read_image_and_metadata
# Computing tools
from utils import adaptive_masking, calculate_worm_properties, crop_image
# Plotting tools
from utils.plotting import masking_summary
# Benchmarking
from utils import downscaling
from utils.benchmarking import tic, toc
# straightening
from utils import create_spline,head2tail_masking
from utils import straighten_image2D_dual

from utils import arc_length
from utils.deprecated import create_skeleton
#from utils.core_utils import create_skeleton

# Import additional libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import time
import os

# Saving the mask
from skimage.io import imsave

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import data_format, verbosity, debugging
from _parameters import xdim, ydim, zdim

# Coding parameters
steps = 'Complete'

# Other coding parameters
svngtm = False # Saving a temporal file

# Setting the parallel mode
parallel = True
n_workers = 12
batch_size = 500

#save retults
results = []
current_res_df = []
image = 0

# list for all channels the stk files in folder
(list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
#print(list_mcherry)

# Scramble
dim = np.shape(list_mcherry)
permutation_array = np.arange(0, dim[0])
permutation_array = np.random.permutation(permutation_array)
list_mcherry = [list_mcherry[k] for k in permutation_array]
list_GFP = [list_GFP[k] for k in permutation_array]

foldername = outputpath+'/Debug_Logs'
# Creating the individual debugging folder
if not os.path.exists(foldername):
    print('MAIN: Creating debugging folder.')
    os.makedirs(foldername)

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

    foldername = outputpath+'/Debug_Logs'

    status = 'File read'

    try:
        # This loop does nothing, just to allow to break the function.
        # A better practice is to return every single time.
        for it in range(0,1):
            # Create mask
            print('MAIN: Adaptive Masking')
            img_binary, additional_outputs  = \
                adaptive_masking(img_mcherry, exp_size = 5)
            
            status = 'Masked Image'

            # Calculating curved worms properties 
            (mean_curved, volume_curved) = \
                calculate_worm_properties(img_binary, img_gfp)
            
            # Storing the properties in current results
            current_res['mean_curved'] = mean_curved  #calculate volume
            current_res['volume_curved'] = volume_curved/(xdim*ydim*zdim)
            
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)

            if steps == 'Masking':
                break

            # Cropping image for further processing
            print('MAIN: Cropping data.')
            cropped_binary, cropped_image = \
                crop_image(img_binary, img_gfp)

            status = 'Cropped Image'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)

            if True:
                print('Select single plane')
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
                calculate_worm_properties(cropped_binary, cropped_image)
            
            # Storing the properties in current results
            current_res['mean_curved'] = mean_curved  #calculate volume
            current_res['volume_curved'] = volume_curved/(xdim*ydim*zdim)
            
            status = 'Cropped Image (Comp)'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
            if steps == 'Cropping':
                break

            # Skelotinisation
            print('MAIN: Skeletonisation and Spline.')
            Xinput, Yinput = create_skeleton(cropped_image, cropped_binary)
            status = 'Skeleton created'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
            
            # Spline fitting
            X, Y, dx, dy = create_spline(Xinput, Yinput)
            status = 'Spline created'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                

            # Cutting off head and tail
            print('MAIN: Head to tail.')
            cropped_binary_ht = head2tail_masking(X,Y,dx,dy,cropped_binary,cut_th=0.2)
            status = 'Cut worm'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
            
            # Calculating properties
            (mean_curved_ht, volume_curved_ht) = calculate_worm_properties(cropped_binary_ht, cropped_image)

            # Storing the properties in current results
            current_res['mean_curved_ht'] = mean_curved_ht
            current_res['volume_curved_ht'] = volume_curved_ht/(xdim*ydim*zdim)

            # Maximum intensity projection
            max_binary = np.max(cropped_binary,0)
            max_image = np.max(cropped_image,0)

            status = 'Cut worm'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
            if steps == 'Skeletonisation':
                break

            # Calculating the width of the worm as 1/10th approximately of the worm length
            length = arc_length(X,Y)

            # Straightening the worm
            print('MAIN: Straightening.')
            (straightened_image, straightened_binary), (xcoord, ycoord) = \
                straighten_image2D_dual((max_image, max_binary), X, Y, dx, dy, width_worm = int(length/10))
            status = 'Straightened worm'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                

            # Calculating intensity straightened worm
            (mean_straightened, volume_straightened) = calculate_worm_properties(straightened_binary,straightened_image)
            
            # Storing the properties in current results
            current_res['mean_straightened'] = mean_straightened
            current_res['area_straightened'] = volume_straightened/(xdim*ydim)

            status = 'Straightened worm (Comp)'
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
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
            # Saving status
            with open(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
                f.write(status)
                
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
        
        status = status+' (after error)'
        os.remove(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt')
        # Saving status
        with open(foldername+'/ERROR_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
            f.write(status)

    stop_alg = toc(start0)

    # Saving temporal files
    if debugging:
        try:
            os.remove(foldername+'/RUNNING_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt')
        except:
            print('Error in removing folder')
        # Saving status
        with open(foldername+'/COMPLETED_Sample'+'_t'+current_res['Frame']+'_s'+current_res['Position']+'_Status.txt', 'w') as f:
            f.write(status + ', with time of '+ str(np.round(stop_alg,2)) + ' (seconds)')

 
    if svngtm:
        print('MAIN: Saving temp file.')
        current_res_df.append(current_res)
        df = pd.DataFrame(current_res_df)
        df.to_csv(outputpath+'/Results_temp_s'+(current_res['Position'])+'_t'+(current_res['Frame'])+'.csv', index=False)

    return current_res


#%%
if parallel:        
    print('MAIN: Running in Parallel mode. Make sure you are running this in the terminal.')
    n_samples = len(list_mcherry)
    n_batches = int(np.ceil(n_samples/batch_size))
    start_parallel = tic()

    for batch_sel in np.arange(0,n_batches):
        print('MAIN: Running batch '+ str(batch_sel) + ' out of '+ str(n_batches))
        #batch_pos = np.arange(batch_sel*batch_size,(batch_sel+1)*batch_size)
        if batch_sel == n_batches-1:
            batch_mcherry = list_mcherry[batch_sel*batch_size:]
            batch_GFP = list_GFP[batch_sel*batch_size:]
        else:
            batch_mcherry = list_mcherry[batch_sel*batch_size:(batch_sel+1)*batch_size]
            batch_GFP = list_GFP[batch_sel*batch_size:(batch_sel+1)*batch_size]
        
        start_batch = tic()
        if __name__=='__main__':
            p = multiprocessing.Pool(n_workers)
            results = p.starmap(main_function,zip(batch_mcherry,batch_GFP))
            p.close()
            p.join()
            

            #time.sleep(10)
            #p.terminate()
            #print('TIMEOUT.', end = ' ')

        print('Batch finished. ', end = ' ')
        stop = toc(start_parallel)
        
        # Saving batch
        df = pd.DataFrame(results)
        df.to_csv(outputpath+'/Results_batch_'+str(batch_sel)+'_'+str(n_batches)+'.csv', index = False)

        print('Pause of 5 min')
        time.sleep(300)

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


# #%% Saving results
# df = pd.DataFrame(results)
# df.to_csv(outputpath+'/Results.csv', index = False)

#%% To run the code
# Additional notes
# import sklearn gave an error. To solve type pip install -U scikit-learn
#
# cd tungstenfs/nobackup/ggrossha/moraluca/SWI_analysis
#
# conda activate SWI
# pip install -e .
#
# python tasks/Running_Code_Adapted_Batch_310519.py

