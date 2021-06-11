#%%
from matplotlib import image
from skimage import io
import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import skeletonize
import cv2 as cv
from fil_finder import FilFinder2D
from utils import tic, toc
import astropy.units as u


tsel = 151
ssel = 15

# Load sample
path = 'C:/Users/moraluca/Desktop/Lin28_test/Output/Single_Plane_Sample_t'+str(tsel)+'_s'+str(ssel)+'/A10_Mask.tiff'
# path = 'C:/Users/moraluca/Desktop/Lin28_test/Output/Single_Plane_Sample_t'+str(tsel)+'_s'+str(ssel)+'/A20_Cropped_mask.tiff'

img = io.imread(path)

imagesel = img[13,:,:].astype(bool)

plt.imshow(imagesel)


# %%
skeletonized_image=skeletonize(imagesel)

plt.imshow(skeletonized_image)

#%%
[X,Y]=np.where(skeletonized_image==1)
    

# %%
plt.imshow(imagesel)

plt.plot(Y,X,'rx')

# %%

dst = cv.cornerHarris(skeletonized_image.astype(np.float32),2,5,0.04)

np.min(dst)
print(np.max(dst))

[Xc,Yc]=np.where(dst*skeletonized_image>.05)

#plt.imshow(dst)


plt.imshow(imagesel)
plt.plot(Y,X,'b.')
plt.plot(Yc,Xc,'rx')


# %%


plt.imshow(imagesel)
plt.plot(Y,X,'b.')
plt.plot(Yc,Xc,'rx')
plt.xlim([1140,1170])
plt.ylim([950,980])

# %%
start = tic()
fil0 = FilFinder2D(imagesel, mask=imagesel) #creates a class, add beamwith?!
# idea, put the fill finder on one plane where the worm is in focus 
_ = toc(start)
#fil0.preprocess_image()
_ = toc(start)

#%%
fil = fil0
start = tic()
fil.create_mask(use_existing_mask=True, verbose=True)
_ = toc(start)
fil.medskel(verbose=False)
_ = toc(start)


fil.analyze_skeletons(branch_thresh=30* u.pix, skel_thresh=10 * u.pix, prune_criteria='length', max_prune_iter= 20)
_ = toc(start)

skeletonized_image=fil.skeleton_longpath
[X,Y]=np.where(skeletonized_image==1)

#
plt.imshow(imagesel)
plt.plot(Y,X,'b.')
# %%
# plt.imshow(fil.mask)
plt.imshow(fil.image)
# plt.imshow(fil.end_pts)
# %%
plt.imshow(fil.skeleton, origin='lower')
plt.contour(fil.skeleton_longpath, colors='r')
# %%
plt.imshow(fil.skeleton, origin='lower')

# %%

print(fil.filaments)

#%%
plt.imshow(fil.filaments[0])

#%%
fil1 = fil.filaments[0]
fil1.skeleton_analysis(fil.image, verbose=True, branch_thresh=5 * u.pix, prune_criteria='length')

# %%
plt.imshow(fil.skeleton, origin='lower')
plt.contour(fil.skeleton_longpath, colors='r')

# %%

# %%
fil = fil0
start = tic()
fil.create_mask(use_existing_mask=True, verbose=True)
_ = toc(start)
fil.medskel(verbose=False)
_ = toc(start)
# %%
plt.imshow(fil.skeleton, origin='lower')
print(fil.number_of_filaments)
# %%
skeletonized_image0=skeletonize(imagesel)

#imgcho = imagesel

filg = FilFinder2D(imgcho, mask=imgcho, beamwidth=0*u.pix) #creates a class, add beamwith?!
fil = filg
start = tic()
fil.create_mask(use_existing_mask=True)
_ = toc(start)
fil.medskel(verbose=False)
_ = toc(start)


fil.analyze_skeletons(branch_thresh=1* u.pix, skel_thresh=10 * u.pix, prune_criteria='length', max_prune_iter= 20)
_ = toc(start)

skeletonized_image=fil.skeleton_longpath
[X,Y]=np.where(skeletonized_image==1)

[X0,Y0]=np.where(fil.skeleton==1)


#
plt.imshow(imagesel)
plt.plot(Y0,X0,'r.')
plt.plot(Y,X,'b.')

#%%
plt.imshow(fil.skeleton, origin='lower')
plt.contour(fil.skeleton_longpath, colors='r')

#%%
long_idx = np.where(np.max(fil.lengths())==fil.lengths())[0][0]
fil1 = fil.filaments[long_idx]
fil1.skeleton_analysis(fil.image, verbose=False, branch_thresh=5 * u.pix, prune_criteria='length')

# %%

fil.number_of_filaments

#%%
fil1.longpath_pixel_coords
# %%
print(fil.output_table)
# %%
from utils import create_spline
fil1 = fil.filaments[long_idx]
[Xinput, Yinput] = fil1.longpath_pixel_coords

X, Y, dx, dy = create_spline(Xinput, Yinput)
# %%
plt.imshow(imagesel)
plt.plot(Y0,X0,'r.')
plt.plot(Yinput,Xinput,'b.')
plt.plot(Y,X,'g.')
# %%
(mean_curved, volume_curved) = calculate_worm_properties(imagesel, imagesel)
# %%
from scipy import interpolate

images_to_be_straightened = np.zeros([2,1200,1200])
images_to_be_straightened[0,:,:] = imagesel
images_to_be_straightened[1,:,:] = imagesel
verbose = True

width_worm = 50
sampling_space = 1
# Debugging and benchmarking
if verbose:
    start = tic()
    print('Straightening the worm. Verbose mode.')

straightened_images = []

# create new coordinate system
# new coord system= old coordinate system + (-dy,dx)*orthogonal coordinate
n = np.arange(-width_worm,width_worm,sampling_space)
xcoord = X[:,None] - n[None,:]*dy[:,None]
ycoord = Y[:,None] + n[None,:]*dx[:,None]

shapex=np.shape(images_to_be_straightened)[2]
shapey=np.shape(images_to_be_straightened)[1]

grid_x, grid_y = np.meshgrid(np.arange(0,shapex),np.arange(0,shapey))
coord_old = np.array([grid_x.reshape(shapex*shapey),grid_y.reshape(shapex*shapey)]).T

yreshape = ycoord.reshape(xcoord.size)
xreshape = xcoord.reshape(ycoord.size)
coord_comb = np.array([yreshape, xreshape])
coord_img = yreshape.astype(int) + 1j*xreshape.astype(int)
coord_img = np.unique(coord_img)
compress_res = np.round(np.unique(coord_img).shape[0] / coord_old.shape[0]*100,1)

# Finding the neighbours
for yk in np.array([-1,1]):
    for xk in np.array([-1,1]):
        coord_img_n = yreshape.astype(int)+yk + 1j*(xreshape.astype(int)+xk)
        coord_img = np.concatenate((coord_img, coord_img_n))
        coord_img = np.unique(coord_img)
        print(coord_img.shape)

compress_res = np.round(np.unique(coord_img).shape[0] / coord_old.shape[0]*100,1)

y_old = np.real(coord_img)
x_old = np.imag(coord_img)

criteria = (y_old<1200)*(y_old>0)*(x_old<1200)*(x_old>0)
y_old = y_old[criteria]
x_old = x_old[criteria]

#reduced_intensity_old = image_to_be_straighten[y_old.astype(int),x_old.astype(int)]
#reduced_coord_old = np.array([y_old,x_old]).T

plt.imshow(imagesel)
plt.plot(y_old,x_old,'r.')


if verbose:
    print('Creating the grid.', end = " ")
    stop = toc(start)
    start = tic()

for k, image_to_be_straighten in enumerate(images_to_be_straightened):
    # if k == 0:
    #     interp = 'nearest'
    #new image
    #intensity_old = image_to_be_straighten.reshape(shapex*shapey)
    reduced_intensity_old = image_to_be_straighten[x_old.astype(int),y_old.astype(int)]
    reduced_coord_old = np.array([y_old,x_old]).T
    _ = toc(start)
    intensity_new = interpolate.griddata(reduced_coord_old, reduced_intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)), method = interp)
    _ = toc(start)
    straightened_image= np.squeeze(intensity_new).reshape(xcoord.shape)
    straightened_image=np.nan_to_num(straightened_image)

    # Stacking
    straightened_images.append(straightened_image)

    if verbose:
        print('Straightening a worm.', end = " ")
        stop = toc(start)
        start = tic()

plt.imshow(straightened_images[1])

# %% Reduced 


# %% Comparions

f = intensity_old[None,:].astype(int)-intensity_old[:,None].astype(int)

plt.imshow(f)
# %%
#from utils.core_utils import straighten_image2D_dual_fast
straightened_images, (xcoord, ycoord) = straighten_image2D_dual_fast(images_to_be_straightened,X,Y,dx,dy)
# %%
plt.imshow(straightened_images[1])

# %%
