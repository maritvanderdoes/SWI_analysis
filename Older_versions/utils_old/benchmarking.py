import timeit
import os
import numpy as np
import signal
from contextlib import contextmanager

def tic():
    '''
    Defines checkpoint for benchmarking and time estimation.
    '''
    start = timeit.default_timer()
    return start

def toc(start, print_elapsed = True):
    '''
    Computes elapsed time from the temporal checkpoint given by start.
    Requires of the tic() function in the utils.benchmarking module.
    '''
    stop = timeit.default_timer()
    elapsed = (stop - start)
    if print_elapsed:
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

# Debug logs
def saving_log(debugpath, state, frame, position, status, runtime = 0):

    # Debug name
    debug_name = debugpath+'/'+state+'_Sample'+'_t'+frame+'_s'+position+'_Status.txt'

    # Defining for Running
    if state == 'RUNNING':
        # Creating
        if status == 'File read':
            with open(debug_name, 'w') as f:
                f.write(status + ', with time of '+ str(np.round(runtime,2)) + ' (seconds)')

        # Appending
        else:
            with open(debug_name, 'a') as f:
                f.write('\n' + status + ', with time of '+ str(np.round(runtime,2)) + ' (seconds)')

    # Renaming if error or completed
    else:
        if os.path.exists(debug_name):
            os.remove(debug_name)
            print('File has been overwritten')

        os.rename(debugpath+'/RUNNING_Sample'+'_t'+frame+'_s'+position+'_Status.txt', debug_name)
        with open(debug_name, 'a') as f:
            f.write('\nFinal time of '+ str(np.round(runtime,2)) + ' (seconds)')

    return status

#%%
def basic_scramble(list_mcherry, list_GFP):
    # Computing the size
    dim = np.shape(list_mcherry)

    # Permuting ordering array
    permutation_array = np.arange(0, dim[0])
    permutation_array = np.random.permutation(permutation_array)

    # Shuffling
    list_mcherry = [list_mcherry[k] for k in permutation_array]
    list_GFP = [list_GFP[k] for k in permutation_array]

    return (list_mcherry, list_GFP)

#%% Generating timeout function

@contextmanager
def timeout(time):
    # Register a function to raise a TimeoutError on the signal.
    signal.signal(signal.SIGALRM, raise_timeout)
    # Schedule the signal to be sent after ``time``.
    signal.alarm(time)

    try:
        yield
    except TimeoutError:
        pass
    finally:
        # Unregister the signal so it won't be triggered
        # if the timeout is not reached.
        signal.signal(signal.SIGALRM, signal.SIG_IGN)

def raise_timeout(signum, frame):
    raise TimeoutError

