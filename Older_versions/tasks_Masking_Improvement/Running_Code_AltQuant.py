#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%
# To run the code type on cmd/terminal:
# python tasks_straightening_single/Running_Code.py

# Setting the directories
from skimage.util.dtype import img_as_float
import _setup

# Loading custom libraries
from utils import image_lists, single_image_lists
from utils.benchmarking import tic, toc, saving_log

# Loading tools
from utils import image_lists, read_image_and_metadata
# Computing tools
from utils import adaptive_masking, calculate_worm_properties, crop_image
from skimage.measure import regionprops, regionprops_table
from scipy.ndimage.measurements import label 
import skimage.morphology as skimorph

# Import additional libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
import os

# load parameters
from _parameters import dirpath, outputpath, channel_GFP, channel_mcherry
from _parameters import parallel, n_workers, batch_size, data_format
from _parameters import debugpath, debugging

# Import main functions
from _main_function_alt import main_function

# %% MAIN CODE
# Initialisation
results = []
current_res_df = []
image = 0

# Selecting sample
s_sel = 4
t_sel = 40

# list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel2 = channel_GFP, s_sel = s_sel, t_sel = t_sel, channel3 = None, data_format = data_format)

#%% Running the code directly here
print('MAIN: Running in Sequential mode.')
for k,files in enumerate(zip(list_mcherry, list_GFP)):
    if k == 0:
        print('Sample selected: '+str(k))
        print('MAIN: Reading image. File selected: '+files[0])
        (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

        # img_binary, _  = adaptive_masking(img_mcherry, mm_th = 6, th_sel = 0.3, exp_size = 3)
        img_binary, _  = adaptive_masking(img_mcherry, mm_th = 4, th_sel = 0.1, exp_size = 3, sorting = True)
        # img_binary, facs  = adaptive_masking(img_mcherry, mm_th = 2.5, th_sel = .9, exp_size = 3, sorting = True)

        (mean_curved, volume_curved) = calculate_worm_properties(img_binary, img_gfp)

        cropped_binary, cropped_image = crop_image(img_binary, img_gfp)
    else:
        print('Potato')

        
#%%

ccs, num_ccs = label(img_binary)
properties=regionprops(ccs,img_gfp,['area','mean_intensity','min_intensity'])


areas = [o.area for o in properties]
selected = ccs == areas.index(max(areas))+1

# krn = skimorph.ball(3)
# surroundings = skimorph.binary_dilation(selected, krn)
# surroundings = surroundings*(np.invert(selected))

krn = skimorph.disk(3)
surroundings = np.zeros(selected.shape)
for pl in np.arange(0,selected.shape[0]):
    surroundings[pl,:,:] = skimorph.binary_dilation(selected[pl,:,:], krn)
surroundings = (surroundings*(np.invert(selected)))==1


#%%
plt.imshow(np.sum(selected,axis = 0))
plt.figure(figsize=(16,16))
plt.imshow(np.max(selected*img_gfp,axis = 0),vmax=200, vmin = 150)
plt.figure(figsize=(16,16))
plt.imshow(np.max(img_gfp,axis = 0),vmax=200, vmin = 100)

#%%
pl_sel = 14
th_sel = 170
selected_plane = selected[pl_sel,:,:]
gfp_plane = img_gfp[pl_sel,:,:]
mch_plane = img_mcherry[pl_sel,:,:]
plt.figure()
plt.plot(mch_plane,gfp_plane,'k.',markersize=3);
plt.plot(mch_plane[selected_plane],gfp_plane[selected_plane],'r.',markersize=3);
plt.plot([np.min(mch_plane),np.max(mch_plane)],[th_sel, th_sel], 'k:')
plt.xlabel('mCherry_signal')
plt.ylabel('GFP_signal')
plt.grid()

plt.figure(figsize=(16,16))
fct = .5
plt.imshow(selected_plane*gfp_plane,vmin=np.min(gfp_plane[selected_plane]),
    vmax=(np.max(gfp_plane[selected_plane]))*fct+np.min(gfp_plane[selected_plane])*(1-fct))

plt.figure(figsize=(16,16))
plt.imshow(selected_plane*gfp_plane>th_sel)

#%%
pl_sel = 13
th_sel = 170
th_sel2 = 1500

# pl_sel = 12
# th_sel = 170
# th_sel2 = 500
selected_plane = selected[pl_sel,:,:]
gfp_plane = img_gfp[pl_sel,:,:]
mch_plane = img_mcherry[pl_sel,:,:]
new_selected = (mch_plane>th_sel2)*selected_plane
plt.figure()
plt.plot(mch_plane,gfp_plane,'k.',markersize=3);
plt.plot(mch_plane[selected_plane],gfp_plane[selected_plane],'r.',markersize=3);
plt.plot([np.min(mch_plane),np.max(mch_plane)],[th_sel, th_sel], 'k:')
plt.plot([th_sel2, th_sel2],[np.min(gfp_plane),np.max(gfp_plane)], 'k:')
plt.xlabel('mCherry_signal')
plt.ylabel('GFP_signal')
plt.grid()

plt.figure(figsize=(16,16))
fct = .5
plt.imshow(selected_plane*gfp_plane,vmin=np.min(gfp_plane[selected_plane]),
    vmax=(np.max(gfp_plane[selected_plane]))*fct+np.min(gfp_plane[selected_plane])*(1-fct))

plt.figure(figsize=(16,16))
plt.imshow(new_selected*gfp_plane>th_sel)

plt.figure(figsize=(16,16))
plt.imshow((gfp_plane>th_sel)*(mch_plane>th_sel2))

#%% Quantiles
quantsel = np.linspace(0.9,1,21)
quantvals = np.quantile(img_gfp[selected],q=quantsel)

plt.figure()
plt.plot(quantvals,quantsel)

th_img = np.zeros(selected.shape)
for k in quantvals:
    # print(str(k)+' '+str(np.sum((selected*img_gfp)>=k)))
    # th_img = th_img + (img_gfp>k)
    th_img = th_img + ((selected*img_gfp)>=k)

#%%
pl_sel = 15
plt.figure(figsize=(16,16))
plt.imshow((selected*img_gfp)[pl_sel,:,:])

plt.figure(figsize=(16,16))
plt.imshow(th_img[pl_sel,:,:])

plt.figure(figsize=(16,16))
plt.imshow(th_img[pl_sel,:,:]>=quantvals[-2])

#%%

plt.figure(figsize=(16,16))
plt.imshow((selected*img_gfp)[pl_sel,:,:]>=quantvals[-1]*.8)


#%%
plt.figure(figsize=(16,16))
plt.imshow(np.max(((selected*img_gfp)>=quantvals[-2]).astype(int),axis = 0))

# %%

plt.imshow(img_binary[15,:,:])

#%%
plt.imshow(np.sum(img_binary,axis = 0))
plt.figure(figsize=(16,16))
plt.imshow(np.max(img_binary*img_gfp,axis = 0),vmax=200, vmin = 150)
plt.figure(figsize=(16,16))
plt.imshow(np.max(img_gfp,axis = 0),vmax=200, vmin = 100)

#%%
plt.plot(facs[3])

#%%
plt.imshow(img_mcherry[10,:,:]>500)

#%%
plt.imshow(img_mcherry[8,:,:])

# %%
plt.plot(np.mean(np.mean(img_mcherry,axis= 2),axis = 1))
plt.plot(img_mcherry[:,100,100])
plt.plot(img_mcherry[:,300,100])
plt.plot(img_mcherry[:,750,600])

# %%
print(np.max(np.sum(img_binary,axis= 0)))
plt.imshow(np.sum(img_binary,axis= 0))
# %%
# plt.imshow(np.sum(img_mcherry,axis= 0))
plt.imshow(np.mean(img_mcherry,axis= 0)>150)

# %%
plt.imshow(np.mean(img_mcherry,axis= 0)>1.5*np.mean(np.mean(img_mcherry,axis= 0)))

# %%
vals_proj = np.mean(img_mcherry,axis= 0)

plt.hist(vals_proj.flatten(),100)
# %%
from skimage.filters import threshold_otsu

modulator = 1
bin = threshold_otsu(vals_proj)
plt.imshow(vals_proj>modulator*bin)
# %%
modulator = 1
img_mcherry2 = img_mcherry*(vals_proj>modulator*bin)+np.min(img_mcherry)
img_binary2, facs  = adaptive_masking(img_mcherry2, mm_th = 2, th_sel = 0.1, exp_size = 3, sorting = False)

# %%
plt.imshow(img_binary2[9,:,:])

# %%
plane = 12
plt.imshow(img_mcherry2[plane,:,:])
plt.figure()
plt.imshow(img_mcherry[plane,:,:])
plt.figure()
plt.imshow(img_gfp[plane,:,:]*img_binary2[plane,:,:])


# %%
plt.semilogx(facs[2])
# %%
plt.imshow(np.sum(img_binary2,axis= 0))

# %%
plt.semilogx(np.diff(facs[2]))

# %%

#%%
plt.imshow(selected[15,:,:])

#%%
plt.imshow(np.sum(selected,axis = 0))

#%%
plt.imshow(surroundings[14,:,:])

# %%

masked = selected*img_gfp
#%%

plt.figure(figsize=(16,16))
plt.imshow(np.max(masked,axis = 0),vmax=1000, vmin = 200)
#%%
masked = selected*img_gfp

plt.figure(figsize=(16,16))
plt.imshow(masked[14,:,:],vmax=800, vmin = 200)

plt.figure(figsize=(16,16))
plt.imshow(masked[14,:,:],vmax=800)



#%%
values_hist = np.arange(80,1200,10)

hist_vals = plt.hist(img_gfp[selected],values_hist, log=True)
hist_dict = dict(zip(values_hist, hist_vals[0]))
current_res = {**current_res, **hist_dict}

plt.figure()
plt.semilogy(hist_vals[0])

#%%
plt.figure()
plt.semilogy(hist_vals[0])
plt.semilogy(hist_vals_samp)

#%%
plt.plot(hist_vals_samp/np.sum(hist_vals_samp)-hist_vals[0]/np.sum(hist_vals[0]))

#%%
plt.figure(figsize=(16,16))
plt.imshow(masked[15,:,:])
# %%
print(np.mean(img_gfp[selected]))
print(np.median(img_gfp[selected]))
print(np.min(img_gfp[selected]))
print(np.mean(img_gfp[np.invert(selected)]))
print(np.median(img_gfp[np.invert(selected)]))
print(np.min(img_gfp[np.invert(selected)]))
print(np.mean(img_gfp[surroundings]))
print(np.median(img_gfp[surroundings]))
print(np.min(img_gfp[surroundings]))

# %%
quantvals_samp = np.quantile(img_gfp[selected],q=np.linspace(0,1,201))


#%%
plt.plot(quantvals_ctrl,quantvals_samp,'x-')
plt.plot(quantvals_samp,quantvals_samp,'k:')
plt.grid()

#%%
plt.plot(quantvals_ctrl,np.linspace(0,1,201),'-')
plt.plot(quantvals_samp,np.linspace(0,1,201),'-')
plt.grid()


#%%


# %% MAIN CODE
# Initialisation
results = []
current_res_df = []
image = 0

# Selecting sample
s_sel = 8
t_sel = None

# list for all channels the stk files in folder
# (list_mcherry, list_GFP) = image_lists(dirpath, channel_mcherry, channel_GFP)
(list_mcherry, list_GFP) = single_image_lists(dirpath, channel_mcherry, channel2 = channel_GFP, s_sel = s_sel, t_sel = t_sel, channel3 = None, data_format = data_format)

#% Running the code directly here
print('MAIN: Running in Sequential mode.')
for k,files in enumerate(zip(list_mcherry, list_GFP)):
    print('MAIN: Reading image. File selected: '+files[0])
    (img_mcherry, img_gfp), current_res = read_image_and_metadata(files, data_format = data_format)

    if np.mod(int(current_res['Frame']),10) == 0:
        print('Sample selected: '+str(k))
        img_binary, _  = adaptive_masking(img_mcherry, mm_th = 2.5, th_sel = 0.3, exp_size = 3)
        # img_binary, _  = adaptive_masking(img_mcherry, mm_th = 2.5, th_sel = 0.1, exp_size = 3, sorting = True)
        # img_binary, facs  = adaptive_masking(img_mcherry, mm_th = 2.5, th_sel = .9, exp_size = 3, sorting = True)

        ccs, num_ccs = label(img_binary)
        properties=regionprops(ccs,img_gfp,['area','mean_intensity','min_intensity'])

        areas = [o.area for o in properties]
        selected = ccs == areas.index(max(areas))+1

        krn = skimorph.ball(3)
        surroundings = skimorph.binary_dilation(selected, krn)
        surroundings = surroundings*(np.invert(selected))

        current_res['mean_curved'] = np.mean(img_gfp[selected])
        current_res['median_curved'] = np.median(img_gfp[selected])
        current_res['min_curved'] = np.min(img_gfp[selected])
        current_res['max_curved'] = np.max(img_gfp[selected])
        current_res['mean_bg'] = np.mean(img_gfp[np.invert(selected)])
        current_res['median_bg'] = np.median(img_gfp[np.invert(selected)])
        current_res['min_bg'] = np.min(img_gfp[np.invert(selected)])
        current_res['max_bg'] = np.max(img_gfp[np.invert(selected)])
        current_res['mean_surr'] = np.mean(img_gfp[surroundings])
        current_res['median_surr'] = np.median(img_gfp[surroundings])
        current_res['min_surr'] = np.min(img_gfp[surroundings])
        current_res['max_surr'] = np.max(img_gfp[surroundings])
        current_res['volume'] = np.sum(selected)

        results.append(current_res)

df = pd.DataFrame(results)
df.to_csv('D:/Results_s'+str(s_sel)+'.csv', index = False)

 # %%
# Synthetic dataset for quantiles
N_cells = 5000
N_signl = int(.1*N_cells)
N_auto  = 1#int(.0*N_cells)
BG_lvl_ctrl  = 55
BG_lvl_smpl = 50
AF_lvl_ctrl = 2*BG_lvl_ctrl
AF_lvl_smpl = 2*BG_lvl_smpl
S2B_cmp = .2 #0.05
SG_lvl  = BG_lvl_smpl + S2B_cmp*BG_lvl_ctrl
nois_bg_ctrl = 4
nois_bg_smpl = 2
nois_sg = nois_bg_ctrl
nois_af_ctrl = 2
nois_af_smpl = 2

CONTROL = BG_lvl_ctrl + nois_bg_ctrl*np.random.randn(N_cells)
SAMPLE  = BG_lvl_smpl + nois_bg_smpl*np.random.randn(N_cells)
CONTROL[-N_auto:] = AF_lvl_ctrl + nois_af_ctrl*np.random.randn(N_auto)
SAMPLE[-N_auto:]  = AF_lvl_smpl + nois_af_smpl*np.random.randn(N_auto)
SAMPLE[0:N_signl]  = SG_lvl + nois_sg*np.random.randn(N_signl)

# Pseudworms
plt.figure(figsize=(8,6))
plt.plot(CONTROL)
plt.plot(SAMPLE)
plt.title('Cell intensities')
plt.xlabel('Pseudoworm cells')
plt.ylabel('Intensity levels')

# Histograms
intensity_values = np.linspace(np.min([CONTROL,SAMPLE]),np.max([CONTROL,SAMPLE]),100)
ctrl_hist = np.histogram(CONTROL,intensity_values)
smpl_hist = np.histogram(SAMPLE ,intensity_values)
plt.figure(figsize=(8,6))
plt.semilogy(intensity_values[1:],ctrl_hist[0])
plt.semilogy(intensity_values[1:],smpl_hist[0])
plt.title('Itensities histogram')
plt.xlabel('Intensity levels')
plt.ylabel('Number of counts')

# Quantiles
quantile_values = np.linspace(0,1,201)
quantvals_ctrl = np.quantile(CONTROL,q=quantile_values)
quantvals_smpl = np.quantile(SAMPLE ,q=quantile_values)
plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantile_values)
plt.plot(quantvals_smpl,quantile_values)
plt.title('Quantile plot')
plt.xlabel('Intensity levels')
plt.ylabel('Quantiles')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantvals_ctrl-80,'k:')
plt.plot(quantvals_ctrl,quantvals_smpl,'rx:')
plt.title('Quantile comparison')
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantile_values,quantvals_smpl - quantvals_ctrl,'k')
plt.title('Quantile substraction')
plt.xlabel('Quantile values')
plt.ylabel('$\Delta$')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_smpl,quantvals_smpl - quantvals_ctrl,'xk')
plt.title('Quantile substraction')
plt.xlabel('Sample')
plt.ylabel('$\Delta$')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantvals_smpl - quantvals_ctrl,'xk')
plt.title('Quantile substraction')
plt.xlabel('Control')
plt.ylabel('$\Delta$')
plt.grid()


plt.figure(figsize=(8,6))
plt.semilogy(quantvals_ctrl[1:],np.diff(quantvals_smpl)/np.diff(quantvals_ctrl),'r-')
plt.title('Quantile comparison')
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()

print((np.sum(np.diff(quantvals_smpl)/np.diff(quantvals_ctrl)>2.5)/quantvals_smpl.shape[0]).round(2))

# %%
N_tests = 50
metrics_ctrl = np.zeros([N_tests,3])
metrics_smpl = np.zeros([N_tests,3])
    
quantile_values = np.linspace(0,1,201)
quantvals_ctrl_it = np.zeros([N_tests, quantile_values.shape[0]])
quantvals_smpl_it = np.zeros([N_tests, quantile_values.shape[0]])


for k in np.arange(0,N_tests):

    CONTROL = BG_lvl_ctrl + nois_bg_ctrl*np.random.randn(N_cells)
    SAMPLE  = BG_lvl_smpl + nois_bg_smpl*np.random.randn(N_cells)
    CONTROL[-N_auto:] = AF_lvl_ctrl + nois_af_ctrl*np.random.randn(N_auto)
    SAMPLE[-N_auto:]  = AF_lvl_smpl + nois_af_smpl*np.random.randn(N_auto)
    SAMPLE[0:N_signl]  = SG_lvl + nois_sg*np.random.randn(N_signl)

    metrics_ctrl[k,0] = np.min(CONTROL)
    metrics_ctrl[k,1] = np.mean(CONTROL)
    metrics_ctrl[k,2] = np.median(CONTROL)

    metrics_smpl[k,0] = np.min(SAMPLE)
    metrics_smpl[k,1] = np.mean(SAMPLE)
    metrics_smpl[k,2] = np.median(SAMPLE)

    quantvals_ctrl_it[k,:] = np.quantile(CONTROL,q=quantile_values)
    quantvals_smpl_it[k,:] = np.quantile(SAMPLE ,q=quantile_values)

#%%
met_chk = 1
avg_metric_ctrl = np.mean(metrics_ctrl[:,met_chk])
std_metric_ctrl = np.std(metrics_ctrl[:,met_chk])
avg_metric_smpl = np.mean(metrics_smpl[:,met_chk])
std_metric_smpl = np.std(metrics_smpl[:,met_chk])
plt.figure(figsize=(8,6))
plt.plot(0,avg_metric_ctrl,'ko')
plt.plot(1,avg_metric_smpl,'ko')
plt.errorbar(0, avg_metric_ctrl,yerr = std_metric_ctrl)
plt.errorbar(1, avg_metric_smpl,yerr = std_metric_smpl)
plt.title('Direct comparison of the mean')
plt.xlabel('Sample')
plt.ylabel('Computed average')
plt.grid()

avg_metric_ctrl = np.mean(metrics_ctrl[:,1]-metrics_ctrl[:,0])
std_metric_ctrl = np.std(metrics_ctrl[:,1]-metrics_ctrl[:,0])
avg_metric_smpl = np.mean(metrics_smpl[:,1]-metrics_smpl[:,0])
std_metric_smpl = np.std(metrics_smpl[:,1]-metrics_smpl[:,0])
plt.figure(figsize=(8,6))
plt.plot(0,avg_metric_ctrl,'ko')
plt.plot(1,avg_metric_smpl,'ko')
plt.errorbar(0, avg_metric_ctrl,yerr = std_metric_ctrl)
plt.errorbar(1, avg_metric_smpl,yerr = std_metric_smpl)
plt.title('Comparison of min-substracted mean')
plt.xlabel('Sample')
plt.ylabel('Computed average')
plt.grid()

avg_metric_ctrl = np.mean(metrics_ctrl[:,1]-metrics_ctrl[:,2])
std_metric_ctrl = np.std(metrics_ctrl[:,1]-metrics_ctrl[:,2])
avg_metric_smpl = np.mean(metrics_smpl[:,1]-metrics_smpl[:,2])
std_metric_smpl = np.std(metrics_smpl[:,1]-metrics_smpl[:,2])
plt.figure(figsize=(8,6))
plt.plot(0,avg_metric_ctrl,'ko')
plt.plot(1,avg_metric_smpl,'ko')
plt.errorbar(0, avg_metric_ctrl,yerr = std_metric_ctrl)
plt.errorbar(1, avg_metric_smpl,yerr = std_metric_smpl)
plt.title('Comparison of median-substracted mean')
plt.xlabel('Sample')
plt.ylabel('Computed average')
plt.grid()

#%%
quantvals_ctrl = np.mean(quantvals_ctrl_it,axis=0)
quantvals_smpl = np.mean(quantvals_smpl_it,axis=0)
plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantile_values)
plt.plot(quantvals_smpl,quantile_values)
plt.title('Quantile plot')
plt.xlabel('Intensity levels')
plt.ylabel('Quantiles')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantvals_ctrl-20,'k:')
plt.plot(quantvals_ctrl,quantvals_smpl,'rx:')
plt.title('Quantile comparison')
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantile_values,quantvals_smpl - quantvals_ctrl,'k')
plt.title('Quantile substraction')
plt.xlabel('Quantile values')
plt.ylabel('$\Delta$')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_smpl,quantvals_smpl - quantvals_ctrl,'xk')
plt.title('Quantile substraction')
plt.xlabel('Sample')
plt.ylabel('$\Delta$')
plt.grid()

plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl,quantvals_smpl - quantvals_ctrl,'xk')
plt.title('Quantile substraction')
plt.xlabel('Control')
plt.ylabel('$\Delta$')
plt.grid()


plt.figure(figsize=(8,6))
plt.plot(quantvals_ctrl[1:],np.diff(quantvals_smpl)/np.diff(quantvals_ctrl),'r-')
plt.title('Quantile comparison')
plt.xlabel('Control')
plt.ylabel('Sample')
plt.grid()

# %%
# %%
N_tests = 50

N_cells = 1000
N_signl = int(.1*N_cells)
N_auto  = int(.1*N_cells)
BG_lvl_ctrl  = 70
BG_lvl_smpl = 50
nois_bg_ctrl = 2
nois_bg_smpl = 2
nois_sg = nois_bg_ctrl
nois_af_ctrl = 2
nois_af_smpl = 2

metrics_ctrl = np.zeros([N_tests,3])
metrics_smpl = np.zeros([N_tests,3])
    
quantile_values = np.linspace(0,1,201)
quantvals_ctrl_it = np.zeros([N_tests, quantile_values.shape[0]])
quantvals_smpl_it = np.zeros([N_tests, quantile_values.shape[0]])

S2B_vec = np.linspace(0,2,26)
A2B_vec = np.linspace(0,2,26)

avg_metric_ctrl = np.zeros([A2B_vec.shape[0],S2B_vec.shape[0],3])
avg_metric_smpl = np.zeros([A2B_vec.shape[0],S2B_vec.shape[0],3])
std_metric_ctrl = np.zeros([A2B_vec.shape[0],S2B_vec.shape[0],3])
std_metric_smpl = np.zeros([A2B_vec.shape[0],S2B_vec.shape[0],3])

for k_af, A2B_cmp in enumerate(A2B_vec):
    AF_lvl_ctrl = BG_lvl_smpl + A2B_cmp*BG_lvl_ctrl
    AF_lvl_smpl = BG_lvl_smpl + A2B_cmp*BG_lvl_smpl

    for k_sg, S2B_cmp in enumerate(S2B_vec):
        SG_lvl  = BG_lvl_smpl + S2B_cmp*BG_lvl_ctrl

        for k in np.arange(0,N_tests):

            CONTROL = BG_lvl_ctrl + nois_bg_ctrl*np.random.randn(N_cells)
            SAMPLE  = BG_lvl_smpl + nois_bg_smpl*np.random.randn(N_cells)
            CONTROL[-N_auto:] = AF_lvl_ctrl + nois_af_ctrl*np.random.randn(N_auto)
            SAMPLE[-N_auto:]  = AF_lvl_smpl + nois_af_smpl*np.random.randn(N_auto)
            SAMPLE[0:N_signl]  = SG_lvl + nois_sg*np.random.randn(N_signl)

            metrics_ctrl[k,0] = np.min(CONTROL)
            metrics_ctrl[k,1] = np.mean(CONTROL)
            metrics_ctrl[k,2] = np.median(CONTROL)

            metrics_smpl[k,0] = np.min(SAMPLE)
            metrics_smpl[k,1] = np.mean(SAMPLE)
            metrics_smpl[k,2] = np.median(SAMPLE)

            # quantvals_ctrl_it[k,:] = np.quantile(CONTROL,q=quantile_values)
            # quantvals_smpl_it[k,:] = np.quantile(SAMPLE ,q=quantile_values)
            
        avg_metric_ctrl[k_af,k_sg,0] = np.mean(metrics_ctrl[:,1])
        std_metric_ctrl[k_af,k_sg,0] = np.std(metrics_ctrl[:,1])
        avg_metric_smpl[k_af,k_sg,0] = np.mean(metrics_smpl[:,1])
        std_metric_smpl[k_af,k_sg,0] = np.std(metrics_smpl[:,1])

        avg_metric_ctrl[k_af,k_sg,1] = np.mean(metrics_ctrl[:,1]-metrics_ctrl[:,0])
        std_metric_ctrl[k_af,k_sg,1] = np.std(metrics_ctrl[:,1]-metrics_ctrl[:,0])
        avg_metric_smpl[k_af,k_sg,1] = np.mean(metrics_smpl[:,1]-metrics_smpl[:,0])
        std_metric_smpl[k_af,k_sg,1] = np.std(metrics_smpl[:,1]-metrics_smpl[:,0])

        avg_metric_ctrl[k_af,k_sg,2] = np.mean(metrics_ctrl[:,1]-metrics_ctrl[:,2])
        std_metric_ctrl[k_af,k_sg,2] = np.std(metrics_ctrl[:,1]-metrics_ctrl[:,2])
        avg_metric_smpl[k_af,k_sg,2] = np.mean(metrics_smpl[:,1]-metrics_smpl[:,2])
        std_metric_smpl[k_af,k_sg,2] = np.std(metrics_smpl[:,1]-metrics_smpl[:,2])

        
# %%

plt.figure()
plt.plot(A2B_vec,avg_metric_smpl[:,:,0]-avg_metric_ctrl[:,:,0])
plt.plot(A2B_vec,A2B_vec*0,'k:')
plt.title('Direct mean')
plt.xlabel('Autoflouroscence-to-background levels')
plt.ylabel('Relative difference of the mean')
plt.grid()

plt.figure()
plt.plot(S2B_vec,(avg_metric_smpl[:,:,0]-avg_metric_ctrl[:,:,0]).T)
plt.plot(S2B_vec,S2B_vec*0,'k:')
plt.title('Direct mean')
plt.xlabel('Signal-to-background levels')
plt.ylabel('Relative difference of the metric')
plt.grid()

plt.figure()
plt.plot(A2B_vec,avg_metric_smpl[:,:,1]-avg_metric_ctrl[:,:,1])
plt.plot(A2B_vec,A2B_vec*0,'k:')
plt.title('Min-norm')
plt.xlabel('Autoflouroscence-to-background levels')
plt.ylabel('Relative difference of the mean')
plt.grid()

plt.figure()
plt.plot(S2B_vec,(avg_metric_smpl[:,:,1]-avg_metric_ctrl[:,:,1]).T)
plt.plot(S2B_vec,S2B_vec*0,'k:')
plt.title('Min-norm')
plt.xlabel('Signal-to-background levels')
plt.ylabel('Relative difference of the metric')
plt.grid()

plt.figure()
plt.plot(A2B_vec,avg_metric_smpl[:,:,2]-avg_metric_ctrl[:,:,2])
plt.plot(A2B_vec,A2B_vec*0,'k:')
plt.title('Median-norm')
plt.xlabel('Autoflouroscence-to-background levels')
plt.ylabel('Relative difference of the mean')
plt.grid()

plt.figure()
plt.plot(S2B_vec,(avg_metric_smpl[:,:,2]-avg_metric_ctrl[:,:,2]).T)
plt.plot(S2B_vec,S2B_vec*0,'k:')
plt.title('Median-norm')
plt.xlabel('Signal-to-background levels')
plt.ylabel('Relative difference of the metric')
plt.grid()


# %%
