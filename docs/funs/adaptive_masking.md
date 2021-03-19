# Adaptive masking
## Function description and inputs
This function takes an 3D input image (z,x,y) and extracts the mask
```
output_mask, additional_info = adaptive_masking(
    input_image, mm_th = 3, th_sel = 0.3, krn_size = 2,
    krn_type = 'Disk', exp_size = 1, fill_holes = True, 
    z_threshold = 0.7, sorting = False, verbose = False):
```
**Inputs**
- <code>input_image</code>
- <code>mm_th</code>

**Outputs**
- <code>output_mask</code>
- <code>additional_info</code>

## Algorithm description
The way to compute the binary mask from <code>input_image</code> follows these steps:
1. If <code>sorting = True</code>, the pixels are sorted in descending order based on intensity levels for each z-plane. Otherwise, the data is just reshaped into a 2D matrix of (z, x\*y)-dimensions.
2. For each (sorted) pixel, the ratio between the maximum and minimum intensity is computed. This metric determines whether a pixel belongs to background or to actual signal (at least in a z-plane). Background pixels have a much lower intensity range, and thus by applying a direct threshold, the background pixels can be removed. Such threshold is defined by the variable <code>mm_th</code>.
3. The signal pixels are determined whether are mask or not at different z-values, based on the threshold defined by <code>th_sel</code>. If the value of intensity for the pixel at a certain z-plane is above the threshold given by <code>th_sel</code>\*min_value_pixel, then it is defined to be mask (1). Otherwise, the pixel at that z-plane is defined not to be mask (0).
4. As the algorithm might not. This is performed by the auxiliary function <code>_mask_postprocessing.py</code>.
5. 

## Data considerations
The algorithm assumes that:
- The shape of the instensity distribution is maintained.
- No major anomalies.

## Masking modes
Depending of the choice of parameters

## Supporting packages
The function requires of the <code>numpy</code>, <code>scikit-image</code> and <code>scipy</code> packages.

## Performance
