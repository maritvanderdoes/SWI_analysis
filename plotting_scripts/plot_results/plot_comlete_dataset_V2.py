#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
#Read both csv files, which contain the GFP quantification and the annotaions made previously (hatch, escape,...)
dir='/Users/Marit/Documents/programming/python/resultstoplot'

CSV_quantification = pd.read_csv(dir+"/resultsall.csv")
CSV_annotation= pd.read_csv(dir+"/goodworms_combined.csv")
min_image_taken=15