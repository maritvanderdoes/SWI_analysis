

import pandas as pd
from marit_functions import *




#folder with images
dirpath= "/Users/Marit/Documents/work/HBL1gfp_worm6"

#channels
channel_GFP='1Lucas-sim-488-561'
channel_mcherry ='2Lucas-sim-561-488'

#save retults
results = []

# list for all channels the stk files in folder
list_mcherry, list_GFP=image_lists_mcherry_GFP(dirpath, channel_mcherry,channel_GFP)


#open mcherry and segment on signal
for (file1,file2) in zip(list_mcherry,list_GFP):
    print(file1)
    
    img_mcherry=read_image(file1)
    img_gfp=read_image(file2)

    #preprocessing and thresholding to make image binary
    img_binary=img_thresholding(img_mcherry)
    
    #select slides based on mean signal
    img_signal,img_binary =select_zslides(img_gfp,img_binary)
    
    
    #calculates properties of the segmented worm
    binary_image, area, mean_intensity, min_intensity=calculate_worm_properties(img_binary,img_signal)
    
    #create overlay of binary image with GFP image
    img_overlay= img_signal*img_binary
    
    #add properties in current results
    current_res = get_meta_info_temp(file2) #get metadata
    current_res['volume'] = area #calculate volume
    current_res['mean_intensity'] = mean_intensity
    current_res['min_intensity']= min_intensity
    current_res['final_intensity']= mean_intensity-min_intensity #calculate intensity
    
    #save in resultarray
    results.append(current_res)
    break
 
 
#save file as csv
df = pd.DataFrame(results)
df.to_csv(dirpath+"/results.csv", index=False)
    
    
