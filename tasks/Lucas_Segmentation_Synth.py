#%% importing
import _setup

from utils import image_lists_mcherry_GFP
from utils import read_image
from utils import img_thresholding
from utils import select_zslides
from utils import calculate_worm_properties
from utils import get_meta_info_temp
from utils import adaptive_masking
from utils import mask_postprocessing

from utils.synthdata import generating_ball

from utils.plotting import dataset_comparison
from utils.plotting import masking_summary

# Import additional libraries
import pandas as pd
import skimage.morphology as skimorph
import matplotlib.pyplot as plt
import numpy as np

#%% 
# 0. Loading random data of same dimensions
# Dimensions
nz, nx, ny = (28, 120, 120)
# Offsets
sz, sx, sy = (0, .1, -.1)
# Threshold level
th_level = .5
# Background level
b_level = 10
# Noise level
n_level = 1e-1

# Parameters
# Min/Max threshold
mm_th = 1.05
# Threshold_selection
th_sel = .5
# Kernel
krn = skimorph.disk(1)

# Generating synthetic dataset
mask_data, test_data, X_ball, Y_ball, Z_ball = generating_ball(nx, ny, nz, sx, sy, sz, th_level, n_level, b_level)

# Running the function
THPX, SORTED, ADPT, PRETH = adaptive_masking(test_data, mm_th, th_sel)

# Mask post-processing
THPX = mask_postprocessing(THPX, krn)

# Comparison of synthetic datasets
dataset_comparison(mask_data, test_data, THPX)

# Presenting outputs
masking_summary(PRETH, SORTED, ADPT, mm_th)
