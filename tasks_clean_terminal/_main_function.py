#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Loading tools
from utils import read_image_and_metadata
# Computing tools
from utils import adaptive_masking, make_mask_background, calculate_worm_properties
from utils import create_skeleton, create_spline, head2tail_masking
from utils import straighten_image2D_dual_fast, arc_length

# Time out
from utils import timeout, saving_log
from utils.benchmarking import tic, toc

# Import additional libraries
from skimage.io import imsave
import numpy as np
import matplotlib.pyplot as plt
import os

# load parameters
from _parameters import data_format, debugging, sngpln
from _parameters import xdim, ydim, zdim   
from _parameters import outputpath, debugpath, data_format
from _parameters import time_out, steps

#%% define the main function
def main_function(list_mcherry,list_GFP):
    with timeout(time_out):
        # For debugging
        start0 = tic()

        # Determining the files to read
        files = (list_mcherry, list_GFP)

        # Reading the image and metadata
        print('MAIN: Reading image. File selected: '+files[0])
        (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

        # Creating debug log
        status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'File read')

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

        try:
            # Create mask
            img_binary, _  = adaptive_masking(img_mcherry)

            # Saving Mask
            if debugging:
                imsave(foldername+'/A10_Mask'+'.tiff',np.float32(255*img_binary), check_contrast = False)
                imsave(foldername+'/A11_Masked_data'+'.tiff',np.float32(img_binary*img_gfp), check_contrast = False)

            # dilate mask
            background_binary= make_mask_background(img_binary)

            # Maximum intensity projection
            img_GFP_max=np.max(img_gfp*img_binary,0)
            img_bin_max=np.max(img_binary,0)

            # Calculating curved worms properties 
            (mean_mip, volume_mip)               = calculate_worm_properties(img_bin_max, img_GFP_max) 
            (mean_curved, volume_curved)         = calculate_worm_properties(img_binary, img_gfp) 
            (mean_background, volume_background) = calculate_worm_properties(background_binary, img_gfp) 
            current_res['mean_curved']     = mean_curved  #calculate volume
            current_res['volume_curved']   = volume_curved
            current_res['mean_background'] = mean_background
            current_res['mean_mip']        = mean_mip
            current_res['volume_mip']      = volume_mip

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'], 'Masked Image', runtime = toc(start0, False))
            if steps == 'Masking': return current_res

            if sngpln:
                print('Select single plane')
                area_cropped = np.sum(np.sum(img_binary,axis=2),axis = 1)
                best_plane = np.where(area_cropped == np.max(area_cropped))[0][0]
                # Selecting the frame
                img_binary = img_binary[best_plane,:,:]
                img_gfp = img_gfp[best_plane,:,:]
                # Reshaping
                img_binary = img_binary[None,:,:]
                img_gfp = img_gfp[None,:,:]
                
            # Calculating curved worms properties 
            (mean_curved_single, area_curved_single) = calculate_worm_properties(img_binary, img_gfp)
            current_res['mean_curved_single'] = mean_curved_single  #calculate volume
            current_res['area_curved_single'] = area_curved_single/(xdim*ydim)
            
            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Cropped Image (Comp)', runtime = toc(start0, False))
            if steps == 'Cropping': return current_res

            # Skelotinisation
            Xinput, Yinput = create_skeleton(img_binary)

            if debugging:
                img_dim = img_binary.shape
                # Saving the skeleton
                fig = plt.figure()
                if img_dim[0] > img_dim[1]:
                    plt.imshow(img_binary[int(img_binary.shape[0]/2),:,:], cmap='gray')
                    plt.plot(Yinput,Xinput,'r.')
                else:
                    plt.imshow(img_binary[int(img_binary.shape[0]/2),:,:].T, cmap='gray')
                    plt.plot(Xinput,Yinput,'r.')
                plt.axis('off')
                plt.title('Old Skeleton')
                plt.savefig(foldername+'/A30_Pre_Skeleton.png')
                plt.close(fig)

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Skeleton created', runtime = toc(start0, False))
            
            # Spline fitting
            X, Y, dx, dy = create_spline(Xinput, Yinput)

            # Saving the new skeleton
            if debugging:
                fig = plt.figure()
                if img_dim[0] > img_dim[1]:
                    plt.imshow(img_binary[int(img_binary.shape[0]/2),:,:], cmap='gray')
                    plt.plot(Y, X, 'r.')
                else:
                    plt.imshow(img_binary[int(img_binary.shape[0]/2),:,:].T, cmap='gray')
                    plt.plot(X, Y, 'r.')
                plt.axis('off')
                plt.title('New Skeleton')
                plt.savefig(foldername+'/A31_Post_Skeleton.png')
                plt.close(fig)

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Spline created', runtime = toc(start0, False))
            if steps == 'Skeletonisation': return current_res

            # Cutting off head and tail
            img_binary_ht = head2tail_masking(X,Y,dx,dy,img_binary,cut_th=0.2)

            if debugging:
                # Saving the cut mask
                fig = plt.figure()
                if img_dim[0] > img_dim[1]:
                    plt.imshow(np.max(img_binary_ht,0),cmap='gray')
                else:
                    plt.imshow(np.max(img_binary_ht,0).T,cmap='gray')
                plt.title('Reduced length by '+str(0.2*100)+"%")
                plt.axis('off')
                plt.savefig(foldername+'/A40_Reduced_worm.png')
                plt.close(fig)
            
            # Calculating properties
            (mean_curved_ht, volume_curved_ht) = calculate_worm_properties(img_binary_ht, img_gfp)
            current_res['mean_curved_ht'] = mean_curved_ht  #calculate volume
            current_res['volume_curved_ht'] = volume_curved_ht/(xdim*ydim*zdim)            

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Cut worm', runtime = toc(start0, False))
            if steps == 'head2tail': return current_res

            # Calculating the width of the worm as 1/10th approximately of the worm length
            length = arc_length(X,Y)

            # Straightening the worm
            max_image=np.max(img_gfp,0)
            max_binary=np.max(img_binary,0)

            (straightened_image, straightened_binary), (xcoord, ycoord) = \
                straighten_image2D_dual_fast((max_image, max_binary), X, Y, dx, dy, width_worm = int(length/10))

            if debugging:
                vmin_val = np.min(img_gfp[img_binary])*.95
                img_dim = xcoord.shape
                img_dim2 = max_image.shape
                dwncoord = np.linspace(0,img_dim[0]-1,20).astype(int)
                
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

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Straightened worm', runtime = toc(start0, False))

            # Calculating intensity straightened worm
            (mean_straightened, area_straightened) = calculate_worm_properties(straightened_binary,straightened_image)
            current_res['mean_straightened'] = mean_straightened  #calculate volume
            current_res['area_straightened'] = area_straightened/(xdim*ydim)

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Straightened worm (Comp)', runtime = toc(start0, False))
            if steps == 'Straightening': return current_res

            #cutting off head and tail, calculating properties
            # print('MAIN: Cutting Straigthened.')
            cutoff=int(straightened_image.shape[0]*0.2)
            straightened_binary_ht = straightened_binary[cutoff:-cutoff,:]
            straightened_image_ht = straightened_image[cutoff:-cutoff,:]

            # Calculating worm properties 
            (mean_straightened_ht, area_straightened_ht)=calculate_worm_properties(straightened_binary_ht, straightened_image_ht)
            current_res['mean_straightened_ht'] = mean_straightened_ht  #calculate volume
            current_res['area_straightened_ht'] = area_straightened_ht/(xdim*ydim)

            # Saving status
            status = saving_log(debugpath, 'RUNNING', current_res['Frame'], current_res['Position'],'Cut straightened', runtime = toc(start0, False))
            if steps == 'Cutting_Straightened': return current_res
        
            #length of worm
            # print('MAIN: Calculating length.')
            length = length/xdim
            bens_volume=np.pi*current_res['area_straightened']**2/(4*length)

            current_res['bens_volume'] = bens_volume

            # Saving logs
            stop_alg = toc(start0)
            status = saving_log(debugpath, 'COMPLETED', current_res['Frame'], current_res['Position'], 'Completed', runtime = stop_alg)


        except:
            print('MAIN: Some computations have not finished.')
            status = saving_log(debugpath, 'ERROR', current_res['Frame'], current_res['Position'], status +' (after error)', runtime = toc(start0, False))


        return current_res
