import pandas as pd
import luigi

from .marit_functions import *  # import * is bad practice ;)


class SWIAnalysisTask(luigi.Task):

    #folder with images
    dirpath = luigi.Parameter()

    #channels
    channel_GFP = luigi.Parameter()  # '1Lucas-sim-488-561'
    channel_mcherry = luigi.Parameter()  # '2Lucas-sim-561-488'

    def run(self):

        #save retults
        results = []

        # list for all channels the stk files in folder
        list_mcherry, list_GFP = image_lists_mcherry_GFP(
            self.dirpath, self.channel_mcherry, self.channel_GFP)

        #open mcherry and segment on signal
        for (file1, file2) in zip(list_mcherry, list_GFP):
            print(file1)

            img_mcherry = read_image(file1)
            img_gfp = read_image(file2)

            #preprocessing and thresholding to make image binary
            img_binary = img_thresholding(img_mcherry)

            #select slides based on mean signal
            img_signal, img_binary = select_zslides(img_gfp, img_binary)

            #calculates properties of the segmented worm
            binary_image, area, mean_intensity, min_intensity = calculate_worm_properties(
                img_binary, img_signal)

            #create overlay of binary image with GFP image
            img_overlay = img_signal * img_binary

            #add properties in current results
            current_res = get_meta_info_temp(file2)  #get metadata
            current_res['volume'] = area  #calculate volume
            current_res['mean_intensity'] = mean_intensity
            current_res['min_intensity'] = min_intensity
            current_res[
                'final_intensity'] = mean_intensity - min_intensity  #calculate intensity

            #save in resultarray
            results.append(current_res)
            break

        #save file as csv
        with self.output().open('w') as out_file:
            df = pd.DataFrame(results)
            df.to_csv(out_file, index=False)

    def output(self):
        return luigi.LocalTarget(self.dirpath + "/results.csv")
