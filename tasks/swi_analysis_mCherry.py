import pandas as pd
import luigi

from skimage.io import imsave

from utils import image_lists, read_image_and_metadata
from utils import adaptive_masking, calculate_worm_properties
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
            print('potato')

            # Reading the image and metadata
            (img_mcherry, img_gfp), meta_out = \
            read_image_and_metadata(files)

            # Running the masking
            binary_mask, sorted_values, pixel_threshold, pixel_range, area_zplane = \
            adaptive_masking(img_mcherry, krn_size = 2)

            # Calculating properties of the segmented worm
            binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
            img_binary = binary_mask, img_signal = img_gfp)

            #add properties in current results
            current_res = meta_out  #get metadata
            current_res['volume'] = area  #calculate volume
            current_res['mean_intensity'] = mean_intensity
            current_res['min_intensity'] = min_intensity
            current_res[ 'final_intensity'] = \
                mean_intensity - min_intensity  #calculate intensity

            #save in resultarray
            results.append(current_res)

            # Save
            imsave(self.outputpath+'\Mask_t'+meta_out['Frame']+'_s'+meta_out['Position']+'.tiff',binary_image)
            

        #save file as csv
        with self.output().open('w') as out_file:
            df = pd.DataFrame(results)
            df.to_csv(out_file, index=False)

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
    

# luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath C:/Users/moraluca/Desktop/Lin28_test --outputpath C:/Users/moraluca/Desktop/Lin28_test/Output --channel-GFP w1Lucas-sim-488-561.stk --channel-mcherry w2Lucas-sim-561-488.stk --local-scheduler