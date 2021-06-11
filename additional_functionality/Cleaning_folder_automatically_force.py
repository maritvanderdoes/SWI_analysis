# Loading libraries
import os
import numpy as np

#%% Selecting the files to remove

# # Base folder
# dirpath = 'C:/Users/moraluca/Desktop/Lin28_test/To_remove'
# # Base name
# basename = 'Lin28_Project_210116'


# Base folder
# dirpath = '/tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/JB_Lin29_210604'
dirpath = '//tungsten-nas.fmi.ch/tungsten/scratch/ggrossha/Lucas/Live_Imaging/JB_Lin29_210604'
# Base name
basename = 'JB_Lin29_210604'

print(dirpath)
print(os.path.exists(dirpath))

# Channels
channel_GFP =   '_1_w1Lucas-sim-488-561'
channel_mcherry = '_1_w2Lucas-sim-561-488'
channel_BF = '_2_Lucas-Brightfield-488-561'
channels = [channel_GFP,channel_mcherry,channel_BF]

#%% Selecting file to remove
s_removal = 93
l_removal = 'pMUT'
p_removal = 33

#%% Creating the filename
k = 1
channel_sel = channels[2]

fullpath = dirpath+'/'+basename+'_t'+str(k)+'_s'+str(s_removal)+\
           '_l'+l_removal+str(p_removal)+channel_sel+'.stk'

print(fullpath)
print(os.path.exists(fullpath))
# %%
dir_log = dirpath+'/'+basename+'_s'+str(s_removal)+'.txt'

# Iterate per timepoint
for k in range(0,241):
    # Iterate per channel
    for c in range(0,3):
        # Check the filename
        filename = basename+'_t'+str(k)+'_s'+str(s_removal)+\
            '_l'+l_removal+str(p_removal)+channels[c]+'.stk'
        
        # Creating the full path
        fullpath = dirpath+'/'+filename

        # Removing files
        if os.path.exists(fullpath):
            os.remove(fullpath)
            with open(dir_log, 'a') as f:
                f.write('\n'+filename+' has been removed.')

        
# %%

# %%

for s_removal in np.arange(67,90):
    p_removal = s_removal-60

    dir_log = dirpath+'/'+basename+'_s'+str(s_removal)+'.txt'
    print('Removing '+basename+'_s'+str(s_removal))

    # Iterate per timepoint
    for k in range(0,241):
        # Iterate per channel
        for c in range(0,3):
            # Check the filename
            filename = basename+'_t'+str(k)+'_s'+str(s_removal)+\
                '_l'+l_removal+str(p_removal)+channels[c]+'.stk'
            
            # Creating the full path
            fullpath = dirpath+'/'+filename

            # Removing files
            if os.path.exists(fullpath):
                os.remove(fullpath)
                with open(dir_log, 'a') as f:
                    f.write('\n'+filename+' has been removed.')
# %%
