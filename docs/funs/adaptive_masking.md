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
- <code>th_sel</code>
- <code>krn_size</code>
- <code>krn_type</code>
- <code>exp_size</code>
- <code>fill_holes</code>
- <code>z_threshold</code>
- <code>sorting</code>
- <code>verbose</code>

**Outputs**
- <code>output_mask</code>
- <code>additional_info</code>
  - <code>sorted_values</code>
  - <code>pixel_threshold</code>
  - <code>pixel_range</code>
  - <code>area_zplane</code>

## Algorithm description
The way to compute the binary mask from <code>input_image</code> follows these steps:
1. If <code>sorting = True</code>, the pixels are sorted in descending order based on intensity levels for each z-plane. Otherwise, the data is just reshaped into a 2D matrix of (z, x\*y)-dimensions.
2. For each (sorted) pixel, the ratio between the maximum and minimum intensity is computed. This metric determines whether a pixel belongs to background or to actual signal (at least in a z-plane). Background pixels have a much lower intensity range, and thus by applying a direct threshold, the background pixels can be removed. Such threshold is defined by the variable <code>mm_th</code>.
3. The signal pixels are determined whether are mask or not at different z-values, based on the threshold defined by <code>th_sel</code>. If the value of intensity for the pixel at a certain z-plane is above the threshold given by <code>th_sel</code>\*min_value_pixel, then it is defined to be mask (<code>1</code>). Otherwise, the pixel at that z-plane is defined not to be mask (<code>0</code>).
4. As the algorithm might not give accurate mask or it could contain holes, a further post-processing step is performed. In this step:
   - If <code>krn_size > 1</code>, the mask is eroded and then dilated to remove spurious pixels. The kernel type is given by <code>krn_type</code>.
   - If <code>fill_holes = True</code>, the holes in the mask are filled.
   - If <code>exp_size > 1</code>, the mask is dilated and then eroded to smooth the data out. The kernel type is given by <code>krn_type</code>. <br>
This is performed by the auxiliary function <code>_mask_postprocessing.py</code>.
5. To ensure there are not spurious dispersed masking areas in the z-planes belonging to the "edges" of the worm, these z-planes are set to zero if the masking area is below the threshold defined by <code>z_threshold</code>\*max_area_in_zplanes. This is performed by the auxiliary function <code>_mask_refinement.py</code>.

## Data considerations
The algorithm assumes that:
- The shape of the instensity distribution is maintained.
- No major anomalies.

## Supporting packages
The function requires of the <code>numpy</code>, <code>scikit-image</code> and <code>scipy</code> packages.

## Performance
### Benchmarking
For an image of dimensions (28, 1200, 1200) on an alienware computer with an intel core i7 8th gen processor and 8GB ram, it takes:
- Sorting/not sorting (Step 1): ~8/1 seconds.
- Thresholding (Step 3): ~1 second.
- Pixel removal (Step 4.a):
  - <code>krn_size = 2</code> ~1 second.
- Filling holes (Step 4.b): ~1 second.
- Mask smoothing (Step 4.c):
  - <code>exp_size = 2</code> ~1 second.
- Z removal (Step 5): ~0.5 seconds.

### Masking modes
A faster performance can lead to a lower quality masking (rough edges, holes, different regions). We will assume there is no sorting, the kernels are given by a disk, <code>th_sel = 0.3</code> and <code>z_threshold = 0.7</code>. Here we can classify different modes:
- **Raw mask**: No pixel removal, no filling holes, no mask smoothing. Fast mode (~ 1.5 second)
- **Just filling holes**: No pixel removal, filling holes, no mask smoothing. Relatively fast mode (~ 3 seconds). Really jagged edges and pixels due to out-of-focus and dispersed light are accepted, giving an abnormal worm outline.
- **Decent fidelity**: Current setup. Pixel removal (size of 2), filling holes, no mask smoothing. Masking takes ~5 seconds per worm. Jagged edges, but not too critical. No spurious individual pixels.
- **Good fidelity**: Pixel removal (size of 2), filling holes, mask smoothing (size of 2). Relatively slow mode (~8 seconds).
- **Excellent fidelity**: Pixel removal (size of >3), filling holes, mask smoothing (size of >3). Slow mode (>10 seconds).
