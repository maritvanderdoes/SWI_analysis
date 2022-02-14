#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Importing required libraries
from os import path
from os import mkdir

#print(np.size(list_mcherry))
def file_validation(dirpath, outputpath, output_filename, 
                    list_mcherry, list_GFP):

    if (path.exists(dirpath)==False):
        raise FileNotFoundError("Image path doesnt exist")
    if (path.exists(outputpath)==False):
        mkdir(outputpath)
        print("New output directory is made")
    if (path.exists(path.join(outputpath,list_mcherry[0]))==False):
        raise FileNotFoundError("Check name mCherry channel")
    if (path.exists(path.join(outputpath,list_GFP[0]))==False):
        raise FileNotFoundError("Check name GFP channel")    
    if (path.exists(path.join(outputpath,output_filename))):
        raise FileExistsError("Watch out, output file already exist. Rename outputfile")

    return