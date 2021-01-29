#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 11:57:35 2020

@author: Marit
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# dir of document
dir='/Users/Marit/Documents/results.csv'

#open resultfile and save it in different variables
result_file = pd.read_csv(dir)
Time= result_file['Frame'].values;
Pos = result_file['Position'].values
intensity= result_file['Intensity'].values
volume= result_file['Volume'].values

#create arrays with sorted data for volume and intensity
array_int = a
array_vol = np.zeros((int(max(Time)), int(max(Pos))))

for i in range(len(Time)):
    array_int[int(Time[i])-1,int(Pos[i])-1]=intensity[i]
    array_vol[int(Time[i])-1,int(Pos[i])-1]=volume[i]
    
    
#set 0 values to nan, to make plotting better
array_int[array_int == 0] = 'nan'
array_vol[array_vol == 0] = 'nan'

# gfp=range(0,1)
# negative=range(51,np.size(rearranged,1))

#wormNr=0
#start=11
#stop=len(rearranged)
#plt.plot(np.arange(0,len(rearranged[start:stop,wormNr])/2,0.5),rearranged[start:stop,wormNr])
#plt.xlabel("time (h)")
#plt.ylabel("intensity (a.u.)")

worm_array_HBL=[0,2,3,4,6,7,8,9,10]
worm_array_control=[16,19,20,21,27,28,29]


for wormNr in worm_array_HBL: #range(np.size(array_int,1)):
     plt.figure(1)
     plt.plot(array_int[:,wormNr]/array_vol[:,wormNr],'.')
     plt.xlabel("time (h)")
     plt.ylabel("intensity (a.u.)")
     plt.ylim([125,200])
     plt.title('HBL-1::GFP')
    

for wormNr in worm_array_HBL:#range(np.size(array_vol,1)):
    plt.figure(2)
    plt.plot(array_vol[:,wormNr],'.')
    plt.xlabel("time (h)")
    plt.ylabel("Volume (pixels)")
    plt.title('HBL-1::GFP')
    
    worm_array=[30]

for wormNr in worm_array_control: #range(np.size(array_int,1)):
     plt.figure(3)
     plt.plot(array_int[:,wormNr]/array_vol[:,wormNr],'.')
     plt.xlabel("time (h)")
     plt.ylabel("intensity (a.u.)")
     plt.ylim([125,200])
     plt.title('control')
    

for wormNr in worm_array_control:#range(np.size(array_vol,1)):
    plt.figure(4)
    plt.plot(array_vol[:,wormNr],'.')
    plt.xlabel("time (h)")
    plt.ylabel("Volume (pixels)")
    plt.title('control')


plt.figure(5)
pos=np.mean(array_int[:,worm_array_HBL]/array_vol[:,worm_array_HBL],1)
neg=np.mean(array_int[:,worm_array_control]/array_vol[:,worm_array_control],1)
plt.plot(pos)
plt.plot(neg)
plt.legend(['HBL-1::GFP','control'])


#np.savetxt("/Users/Marit/Downloads/output.csv", rearranged, delimiter=",")
    