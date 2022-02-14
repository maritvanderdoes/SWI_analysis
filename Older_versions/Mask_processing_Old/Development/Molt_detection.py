
# %% import packages
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from sklearn.decomposition import PCA


# %%
dir='D:/GitHub/SWI_analysis/Analysis V3/'
filename = 'cleanup_Results_LJ_Lin28_Project_210305.csv'
# filename = 'Worm_segmentation_results_Michi.csv'

# Opening the file
preanalysis_csv = pd.read_csv(dir+filename)
preanalysis_csv = preanalysis_csv.sort_values(['Position','Frame'])

# Checking NAs
print(preanalysis_csv.isnull().sum())

# %%
position_sel = 76 # 50, 68-74
motion_level = 7
motion_window = 3
ec_th = 0.9
vol_change_th = 8

selection_csv = preanalysis_csv[preanalysis_csv['Position']==position_sel]
print(selection_csv.isnull().sum())

# New metrics
selection_csv['displacement'] = \
    np.sqrt((selection_csv['mask_centroid_x'].diff())**2+\
            (selection_csv['mask_centroid_y'].diff())**2)
selection_csv['volume_diff'] = selection_csv['mask_area'].diff()

#%% Processing
selection_csv['mask_area_log'] = np.log(selection_csv['mask_area'])
selection_csv['mask_area_log_average'] = selection_csv['mask_area_log'].rolling(5, center = True).median()
#%%
th = 0.1
outliers = selection_csv['mask_area_log']-selection_csv['mask_area_log_average']>th
plt.plot(selection_csv['Frame'], selection_csv['mask_area_log'])
plt.plot(selection_csv['Frame'], selection_csv['mask_area_log_average'])

plt.figure()
plt.plot(selection_csv['Frame'], outliers)

#%%
selection_csv.loc[outliers,'mask_area_log'] = selection_csv['mask_area_log_average']
volume_smooth = selection_csv['mask_area_log'].rolling(5,center = True).mean()
plt.figure()
plt.plot(selection_csv['Frame'], volume_smooth.diff())
plt.ylim([-.1, .1])

plt.figure()
plt.plot(selection_csv['Frame'], volume_smooth.diff().rolling(5,center = True).mean().diff().rolling(5,center = True).mean())
plt.ylim([-.1, .1])

#%%
window_small = 5
window_large = 15
plt.figure()
plt.plot(selection_csv['Frame'], selection_csv['displacement'].rolling(window_small, center=True).mean())

plt.figure()
plt.plot(selection_csv['Frame'], selection_csv['displacement'].rolling(window_small, center=True).mean()-selection_csv['displacement'].rolling(window_large, center=True).mean())


#%%
plt.plot(np.log(selection_csv['displacement'].rolling(5, center=True).mean()),volume_smooth.diff(), 'x' )
plt.ylim([-.1, .1])

selection = (np.log(selection_csv['displacement'].rolling(5, center=True).mean()) < 10 )&\
 (volume_smooth.diff()<.0)
plt.figure()
plt.plot(selection_csv['Frame'], selection_csv['mask_area_log'])
plt.plot(selection_csv.loc[selection,'Frame'], selection_csv.loc[selection,'mask_area_log'],'r.')

# %%

plt.plot(selection_csv['volume_diff'])
# %%
pca_csv = np.log(selection_csv[['displacement', 'mask_area']]).diff()
pca_csv['mask_eccentricty'] = selection_csv['mask_eccentricity'].diff()
pca_csv = pca_csv.dropna()
pca_csv = pca_csv - pca_csv.min()/(pca_csv.max() - pca_csv.min())
pca_trf = PCA(n_components=3).fit_transform(pca_csv)

plt.plot(pca_trf[:,0],pca_trf[:,1],'x')
plt.figure()
plt.plot(pca_trf[:,0],pca_trf[:,2],'x')
plt.figure()
plt.plot(pca_trf[:,1],pca_trf[:,2],'x')

plt.figure()
sns.histplot(pca_trf[:,0])

# %%
pca_csv = np.log(selection_csv[['displacement', 'mask_area']]).diff()

plt.plot(pca_csv['displacement'], pca_csv['mask_area'],'x ')
# %%
s_win = 1
l_win = 17
th_sel = -.01
selection_csv['mask_area_log_s'] = selection_csv['mask_area_log'].rolling(5,center = True).mean()

h_fun = volume_smooth.diff().rolling(s_win,center = True).median()
l_fun = volume_smooth.diff().rolling(l_win,center = True).median()
plt.figure()
plt.plot(selection_csv['Frame'], h_fun)
plt.plot(selection_csv['Frame'], l_fun)
plt.ylim([-.1, .1])
plt.grid()


plt.figure()
plt.plot(selection_csv['Frame'], h_fun-l_fun)
plt.axhline(th_sel,color = 'r')
plt.ylim([-.1, .1])
plt.grid()

selection = ( h_fun-l_fun)< th_sel

plt.figure()
plt.plot(selection_csv['Frame'], selection_csv['mask_area_log'])
plt.plot(selection_csv.loc[selection,'Frame'], selection_csv.loc[selection,'mask_area_log'],'r.')

# %%
