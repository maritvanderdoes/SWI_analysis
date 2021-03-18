#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Importing required libraries
import skimage.filters as skf
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import glob
import os
import re
from skimage import io
from scipy.ndimage.measurements import label 
from skimage.measure import regionprops

#-----------------------------------------------------------------------------
# Generating synthetic datasets

# Generating a simple ball
def generating_ball(nx, ny, nz, sx, sy, sz, th_level, n_level, b_level):
    x = np.linspace(-1, 1, nx) + sx
    y = np.linspace(-1, 1, ny) + sy
    z = np.linspace(-1, 1, nz) + sz

    Z = np.zeros([nz,nx,ny])
    X = np.zeros([nz,nx,ny])
    Y = np.zeros([nz,nx,ny])

    for z_p in range(0,nz):
        Xt, Yt = np.meshgrid(x, y)
        X[z_p,:,:] = Xt
        Y[z_p,:,:] = Yt
        Z[z_p,:,:] = z[z_p]

    mask_data = (np.sqrt(Z**2+X**2+Y**2)<th_level) 
    test_data = mask_data + n_level*np.random.randn(nz,nx,ny) + b_level

    return mask_data, test_data, X, Y, Z

