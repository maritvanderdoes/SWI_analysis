import timeit
import numpy as np

def tic():
    start = timeit.default_timer()
    return start

def toc(start):
    stop = timeit.default_timer()
    print('Time elapsed of ' +"{:.2f}".format((stop - start))+' seconds')
    return stop

def downscaling(images, xy_factor = 2, z0 = None, zf = None, verbose = False):
    # Checking for null parameters
    if z0 == None:
        z0 = 0

    # Single image
    if np.ndim(images) == 3:
        # Checking for null parameters
        if zf == None:
            zf = images.shape[0]
        # Downscaling
        downscaled_images = images[z0:zf,::xy_factor,::xy_factor]

    # Set of images
    elif np.ndim(images) == 4:
        downscaled_images = []
        for k in enumerate(images):
            # Checking for null parameters
            if zf == None:
                zf = images[k[0]].shape[0]
            # Downscaling
            downscaled_images.append(images[k[0]][z0:zf,::xy_factor,::xy_factor])

    # Single slide
    elif np.ndim(images) == 2:
        downscaled_images = images[::xy_factor,::xy_factor]

    # Giving output
    if verbose:
        if np.ndim(images)<4:
            print('Number of images: ' + str(1))
            print('Original size: '+str(images.shape))
            print('     New size: '+str(downscaled_images.shape))
        elif np.ndim(images)==4:
            print('Number of images: ' + str(len(images)))
            print('Original size: '+str(images[0].shape))
            print('     New size: '+str(downscaled_images[0].shape))
            print('*Assuming identical size')

    return downscaled_images