import pandas as pd
import numpy as np
import luigi

from skimage.io import imsave

from utils import image_lists, read_image_and_metadata
from utils import adaptive_masking, calculate_worm_properties
from utils import tic, toc
  # import * is bad practice ;)


class SWIAnalysisTask(luigi.Task):

    #folder with images
    dirpath = luigi.Parameter()

    #output folder
    outputpath= luigi.Parameter()

    #channels
    channel_GFP = luigi.Parameter()  # '1Lucas-sim-488-561'
    channel_mcherry = luigi.Parameter()  # '2Lucas-sim-561-488'

    def run(self):
        #save retults
        results = []

        # list for all channels the stk files in folder
        (list_mcherry, list_GFP) = image_lists(
            self.dirpath, self.channel_mcherry, self.channel_GFP)

        #open mcherry and segment on signal
        for files in zip(list_mcherry, list_GFP):
            print(files[0])
            # Checking times
            start = tic()

            # Reading the image and metadata
            # images_out[0] = mCherry, images_out[1] = GFP
            images_out, current_res = read_image_and_metadata(files)

            # Running the masking
            binary_mask, additional_info = adaptive_masking(images_out[0])

            # Calculating properties of the segmented worm
            # Metrics are: area, mean_intensity, min_intensity and (centroid)
            binary_image, metrics = calculate_worm_properties(binary_mask, images_out[1])

            # Compute masked data
            masked_data = images_out[1] * binary_image

            # Straightening

            # Head/Tail processing

            #add properties in current results
            current_res.update(dict(zip(('volume','mean_intensity','min_intensity'), metrics[0:3])))
            current_res['final_intensity'] = metrics[1] - metrics[2]  #calculate intensity
            current_res.update(dict(zip(('centroid_z','centroid_x','centroid_y'), metrics[3])))

            # Checking time
            current_res['toc'] = toc(start)

            #save in resultarray
            results.append(current_res)

            # Save Mask
            imsave(self.outputpath+'\Mask_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',255*binary_image, check_contrast = False)
            # Save Masked data
            imsave(self.outputpath+'\Masked_data_t'+current_res['Frame']+'_s'+current_res['Position']+'.tiff',np.float32(masked_data), check_contrast = False)        

        #save file as csv
        with self.output().open('w') as out_file:
            df = pd.DataFrame(results)
            df.to_csv(out_file, index=False, line_terminator='\n')

    def output(self):
        return luigi.LocalTarget(self.outputpath + "/results.csv")


# if __name__ == '__main__':
#     luigi.build([
#         SWIAnalysisTask(
#             dirpath="/Users/Marit/Documents/work/HBL1gfp_worm6",
#             outputpath="/Users/Marit/Documents",
#             channel_GFP="w1Lucas-sim-488-561",
#             channel_mcherry="w2Lucas-sim-561-488")
#     ],
#                 local_scheduler=True)

# if __name__ == '__main__':
#     luigi.build([
#         SWIAnalysisTask(
#             dirpath="C:/Users/moraluca/Desktop/Lin28_test",
#             outputpath="C:/Users/moraluca/Desktop/Lin28_test/Output",
#             channel_GFP="w1Lucas-sim-488-561",
#             channel_mcherry="w2Lucas-sim-561-488")
#     ],
#                 local_scheduler=True)
    

# luigi --module tasks.swi_analysis_mCherry_troubleshooting SWIAnalysisTask --dirpath C:/Users/moraluca/Desktop/Lin28_test --outputpath C:/Users/moraluca/Desktop/Lin28_test/Output --channel-GFP w1Lucas-sim-488-561.stk --channel-mcherry w2Lucas-sim-561-488.stk --local-scheduler