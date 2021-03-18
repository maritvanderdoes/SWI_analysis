import pandas as pd
import luigi

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

            # Reading the image and metadata
            # images_out[0] = mCherry, images_out[1] = GFP
            images_out, current_res = read_image_and_metadata(files)

            # Running the masking
            binary_mask, additional_info = adaptive_masking(images_out[0])

            # Calculating properties of the segmented worm
            # Metrics are: area, mean_intensity and min_intensity
            binary_image, metrics = calculate_worm_properties(binary_mask, images_out[1])

            # Straightening

            # Head/Tail processing

            
            #add properties in current results
            current_res.update(dict(zip(('volume','mean_intensity','min_intensity'), metrics[0:3])))
            current_res['final_intensity'] = metrics[1] - metrics[2]  #calculate intensity

            #save in resultarray
            results.append(current_res)


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
    

# luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath C:/Users/moraluca/Desktop/Lin28_test --outputpath C:/Users/moraluca/Desktop/Lin28_test/Output --channel-GFP w1Lucas-sim-488-561.stk --channel-mcherry w2Lucas-sim-561-488.stk --local-scheduler