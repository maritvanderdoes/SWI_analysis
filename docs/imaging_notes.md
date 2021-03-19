# Imaging notes

## Filename convection
Files are assumed to be of the form of <code>basename_tX_sY_additionainfo.extension</code>, where <code>X</code> is a number indicating the frame and <code>Y</code> is a number indicating the stage position. Channel and worm label are encoded in <code>additionalinfo</code>. Allowed extensions are <code>.stk</code> and <code>.tiff</code>.

Alternatively, files of the form <code>basename_sY_tX_additionainfo.extension</code> can be read, by setting the <code>data_format = 'st'</code> in <code>utils.read_image_and_metadata.py</code>.

## Recording features
The code has been developed with data acquired as follows:
- First green and red channels are recorded _simultaneously_ 🟢🔴
- Then the brightfield channel is recorded _sequentially_ ⚪

## Imaging workflow \[Grosshans' lab user\]
As a microscope and tungsten user:
1. Copy the raw data into <code>scratch</code>.
2. Quality check.
3. Compress the data using the <code>StkToTifImageCompressionWorkflow</code> workflow from <code>l1048.fmi.ch/frontend/jobs/</code>.
4. Delete the raw data once the compressed data has been checked.
5. Send the compressed data to tape.
6. If the data is now having a reduced ussage, send it to <code>nobackup</code>.

## Technical considerations
