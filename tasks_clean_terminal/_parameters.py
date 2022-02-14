#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
To ensure the code is as compatible between both users, Marit and 
Lucas, this function takes care of setting that for both of them.

To ensure clarity, the code also set ups the different parameters
that do not require to be changed constantly.

Information
--------------------------------
dirpath: Directory where the files are
data_format: Order in which the position and the time are present
    in the filename.
    By default, it is 'st'.
verbosity: Do you want the functions to spit all possible information?
--------------------------------

Below are different parameters for the masking functions, but these
might no be necessary in the case of default values for functions.
'''

# Loading required packages
import os

# load parameters
dirpath         = '/tungstenfs/scratch/ggrossha/Marit/HetP_Quant/20201224_HBL1'
outputpath      = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_HetPathQuant_Results'
output_filename = 'Results_20201224_HBL1_test.csv'
debugpath       = outputpath + '/Debug'
channel_GFP     = 'w1Lucas-sim-488-561'
channel_mcherry = 'w2Lucas-sim-561-488'

# Additional parameters
n_workers = 12
data_format = 'ts'
steps = 'Complete'
debugging = True
verbosity = True
time_out = 200 # seconds

# Scaling factors
xdim = 0.533
ydim = 0.533
zdim = 2

# Single plane
sngpln = True

# # Create output folder
# outputpath = os.path.abspath(os.path.join(dirpath,'Output'))
# if os.path.exists(outputpath):
#     print('Output folder already existing')
# else:
#     print('Folder output has been created')
#     os.makedirs(os.path.abspath(os.path.join(dirpath,'Output')))

# # Create debug log
# debugpath = outputpath+'/Debug_Logs'
# if not os.path.exists(debugpath):
#     print('Debug folder has been created.')
#     os.makedirs(debugpath)

# Printing rest of summary
print('Selected data_format: ' + data_format)
print('Verbosity mode: ' + str(verbosity))
print('Debugging mode: ' + str(debugging))
print('Output directory: ' + str(outputpath))
