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
This code takes an 3D input image (z,x,y) and binarises it.

## Performance
