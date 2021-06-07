#%%
#import data
from utils import image_lists, read_image_and_metadata
from utils import adaptive_masking, calculate_worm_properties
from skimage.measure import regionprops_table
from scipy.ndimage.measurements import label 
import numpy as np
import pandas as pd
from skimage.morphology import skeletonize
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from scipy import interpolate
import matplotlib.pyplot as plt
import timeit




dir='/Users/Marit/Documents/work/test_images'
(list_mcherry, list_GFP) = image_lists(dir, 'w2Lucas-sim-561-488', 'w1Lucas-sim-488-561')
(img_mcherry, img_gfp), current_res = read_image_and_metadata((list_mcherry[0], list_GFP[0]))
img_binary, crap = adaptive_masking(img_mcherry)
(cropped_image, cropped_binary, volume_old, mean_intensity_old, min_intensity) = calculate_worm_properties(img_binary, img_gfp)

slidenr=4
binary_to_test=cropped_binary[slidenr,:,:]
image_to_test=cropped_image[slidenr,:,:]

#import fuctions
def create_skeleton2(mask, verbose=False):
    
    zslide_to_focus=np.int(mask.shape[0]/2)

    image=mask[zslide_to_focus,:,:]
    skeletonized_image=skeletonize(image)
    [X,Y]=np.where(skeletonized_image==1)
    
    if (verbose==True):
        plt.figure()
        plt.imshow(image, cmap='gray')
        plt.contour(skeletonized_image, colors='r')
        plt.axis('off')
        plt.title('check skeleton')
        plt.show()
    
    return X,Y

def create_spline(X,Y,smoothening, nrspline,verbose=False):
    """Creates a spline through the X and Y coordinates by a number of splinepoints.
    This gives as output the rescaled x and y, and the derivatives of the line.
    Idea explained from 
    https://stackoverflow.com/questions/37742358/sorting-points-to-form-a-continuous-line

    Args:
        X ([1D array]): xcoordinates through with spline need to be splitted
        Y ([1D array]): ycoordinates through with spline need to be splitted
        verbose (bool, optional): [description]. Defaults to False.

    Returns:
        x_new: 
        y_new:
        dx:
        dy:
    """

    #variables
    s=smoothening
    length=_arc_length(X,Y)
    nsplinepoints=int(X.shape[0]/nrspline)
    #np.int(length/(X.shape[0]*nrspline)) # optimize number of splinepoints

    inputarray = np.c_[X, Y]

    # runs a nearest neighbors algorithm on the coordinate array
    clf = NearestNeighbors(n_neighbors=2).fit(inputarray)
    G = clf.kneighbors_graph()
    T = nx.from_scipy_sparse_matrix(G)

    # sorts coordinates according to their nearest neighbors order
    order = list(nx.dfs_preorder_nodes(T, 0))
    xx = X[order]
    yy = Y[order]

    # Loops over all points in the coordinate array as origin, determining which results in the shortest path
    paths = [list(nx.dfs_preorder_nodes(T, i)) for i in range(len(inputarray))]
    mindist = np.inf
    minidx = 0

    for i in range(len(inputarray)):
        p = paths[i]           # order of nodes
        ordered = inputarray[p]    # ordered nodes
        # find cost of that order by the sum of euclidean distances between points (i) and (i+1)
        cost = (((ordered[:-1] - ordered[1:])**2).sum(1)).sum()
        #print(cost)
        if cost < mindist:
            mindist = cost
            minidx = i

    opt_order = paths[minidx]

    xxx = X[opt_order]
    yyy = Y[opt_order]

    # fits a spline to the ordered coordinates
    tckp, u = interpolate.splprep([xxx, yyy], s=s, k=3, nest=-1)
    x_new, y_new = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp)
    dx, dy = interpolate.splev(np.linspace(0,1,nsplinepoints), tckp,der=1)

    #calculating normalized derivative points
    magnitude=np.sqrt(dx**2+dy**2)
    dx=dx/magnitude
    dy=dy/magnitude

    if (verbose==True):
        plt.plot(x_new,y_new)
        plt.title('fittet spline')

    return x_new,y_new, dx,dy, nsplinepoints

def straighten_image2D(image_to_be_straighten,X,Y,dx,dy,nrspline, verbose=False):
    width_worm=100  #find optimal worm_with

    #create new coordinate system
    # new coord system= old coordinate system + (-dy,dx)*new coordinate system
    n = np.arange(-width_worm,width_worm,nrspline)
    xcoord = X[:,None] - n[None,:]*dy[:,None]
    ycoord = Y[:,None] + n[None,:]*dx[:,None]

    shapex=image_to_be_straighten.shape[1]
    shapey=image_to_be_straighten.shape[0]
    
    grid_x, grid_y = np.meshgrid(np.arange(0,shapex),np.arange(0,shapey))
    coord_old = np.array([grid_x.reshape(shapex*shapey),grid_y.reshape(shapex*shapey)]).T

    #new image
    intensity_old = image_to_be_straighten.reshape(shapex*shapey)
    intensity_new = interpolate.griddata(coord_old, intensity_old, (ycoord.reshape(xcoord.size), xcoord.reshape(ycoord.size)))
    straightened_image= np.squeeze(intensity_new).reshape(xcoord.shape)
    straightened_image=np.nan_to_num(straightened_image)

   

    if (verbose==True):
        plt.figure()
        plt.imshow(img_BIN.T)
        plt.plot(X,Y)
        plt.plot(xcoord.T,ycoord.T)
        plt.title("original image")
        plt.show()

        plt.figure()
        plt.imshow(straightened_image_BIN)
        plt.title("straightened image")
        plt.show()
    
    return straightened_image

def _arc_length(x, y):
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)

    return arc

def head2tail_masking(X,Y,dx,dy,img_binary,cut_th=0.2):
    # binary_to_test
    # Define the points
    grid_y, grid_x = np.meshgrid(np.arange(0,img_binary.shape[1]),np.arange(0,img_binary.shape[0]))

    # Find the first line
    points2mask = np.linspace(0,1,np.shape(X)[0])

    # Find the upperline
    upper = np.sum(points2mask<cut_th)
    lower = np.sum(points2mask<(1-cut_th))

    # Reference point
    ref_up = X[upper]*dx[upper]+Y[upper]*dy[upper]
    ref_low = X[lower]*dx[lower]+Y[lower]*dy[lower]

    # define the upper points to be zero
    fun_u = (grid_y*dy[upper]+grid_x*dx[upper])
    fun_l = (grid_y*dy[lower]+grid_x*dx[lower])
    binary_upper = fun_u>ref_up
    binary_lower = fun_l<ref_low

    #new image
    binary_new = binary_upper*binary_lower


    # plt.figure()
    # plt.contour((fun_u)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')
    # plt.figure()
    # plt.imshow((mat_upper)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')

    # plt.figure()
    # plt.contour((fun_l)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')
    # plt.figure()
    # plt.imshow((mat_lower)*img_binary)
    # plt.plot(Y[upper],X[upper],'ro')
    # plt.plot(Y[lower],X[lower],'rx')

    if (verbose==true)
        plt.figure()
        plt.imshow(mat_val * img_binary)
        plt.plot(Y[upper],X[upper],'ro')
        plt.plot(Y[lower],X[lower],'rx')


    return binary_new

# %%
Xinput,Yinput=create_skeleton2(cropped_binary, verbose=True)
(X,Y, dx,dy,nrsplinepoints)=create_spline(Xinput,Yinput, 100, 1,verbose=False)
binary=head2tail_masking(Xinput,Yinput,dx,dy,binary_to_test)

#%%

#calculating old properties of non straightened image
ccs, num_ccs = label(image_to_test)
properties_old=regionprops_table(ccs,image_to_test,['area','mean_intensity']) #calculates the properties of the different areas

#create spline
Xinput,Yinput=create_skeleton2(cropped_binary, verbose=True)
result = []


#loop over spline nr and smoothening numbers
for nrspline in [1]:
    print("spline_input ="+ str(nrspline))
    start = timeit.default_timer()
    for smoothening in [100]:
        print("smoothening ="+ str(smoothening))
        #create spline and straihgten image based on this spline
        (X,Y, dx,dy,nrsplinepoints)=create_spline(Xinput,Yinput,smoothening, nrspline,verbose=False)
        simage_GFP, simage_BIN=straighten_image2D(image_to_test,binary_to_test,X,Y,dx,dy, nrspline, verbose=True)

        ccs, num_ccs = label(simage_BIN)
        properties_new=regionprops_table(ccs,simage_GFP,['area','mean_intensity']) #calculates the properties of the different areas

        error_area=(properties_old['area'][0]-(properties_new['area'][0]*nrspline**2))/properties_old['area'][0]*100
        error_intensity=(properties_old['mean_intensity'][0]-properties_new['mean_intensity'][0])/properties_old['mean_intensity'][0]*100

        result.append([nrsplinepoints,nrspline,smoothening,error_area,error_intensity])

    stop = timeit.default_timer()
    elapsed = (stop - start)
    print('Elapsed time of ' +"{:.2f}".format(elapsed)+' seconds')





# %%
import pandas as pd
import seaborn as sns
df = pd.DataFrame(result, columns=['nrsplinepoints','nrspline', 'smoothing', 'error_area','error_intenisty'])
df_unmelted = df.pivot(index='nrspline', columns='smoothing')

plt.figure()
plt.title('error in area')
sns.heatmap(df_unmelted.error_area,annot=True )

plt.figure()
plt.title('error in mean_intensity')
sns.heatmap(df_unmelted.error_intenisty, annot=True)
# %%
plt.figure()
plt.title('error in mean_intensity')
sns.heatmap(df_unmelted.error_intenisty*df_unmelted.error_area, annot=True)


plt.imshow(smpl1-smpl4)



# %%

