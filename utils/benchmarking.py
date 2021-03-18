import timeit
import numpy as np

def tic():
    '''
    Defines checkpoint for benchmarking and time estimation.
    '''
    start = timeit.default_timer()
    return start

def toc(start):
    '''
    Computes elapsed time from the temporal checkpoint given by start.
    Requires of the tic() function in the utils.benchmarking module.
    '''
    stop = timeit.default_timer()
    elapsed = (stop - start)
    print('Elapsed time of ' +"{:.2f}".format(elapsed)+' seconds')
    return elapsed

def downscaling(images, xy_factor = 2, z0 = None, zf = None, verbose = False):
    '''
    Downscales (3D) image(s) by a given factor for the xy plane 
    symmetrically and selects a range of z planes.

    It allows for more than an image. In this case, the output is given 
    as a tuple.

    Note that the downscaling factor is symmetrical and it is represented
    by the "side" factor. The actual number of pixels is reduced by the
    square of the factor.

    Parameters
    ----------
    images: Set of images to downscale.

    xy_factor: Downscaling factor. By defaul it downscales the image by
    a factor of two.

    [z0, zf]: range of selected z planes. By default selects the whole
    range.

    verbose: Output progression through steps and time taken. By default
        it is off.
        
    Returns
    -------
    downscaled_images: Downscale (3D) image(s).

    '''

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