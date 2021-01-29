import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Read both csv files, which contain the GFP quantification and the annotaions made previously (hatch, escape,...)
dir='/Users/Marit/Documents/programming/python/resultstoplot'

CSV_quantification = pd.read_csv(dir+"/resultsall.csv")
CSV_annotation= pd.read_csv(dir+"/goodworms_combined.csv")
min_image_taken=15

#Indicate the different conditions, analyzed in the perforemed SWI experiment, and the colors you want to plot them
conditions=['control','hbl1','lin14','lin29']
colors = ['black','blue','green','red']
colors2= ['blue','black','green','red']
#Initialize the plot

plt.figure(1)

listtime=[]
listGFP=[]
listcondition=[]
wormtoplot=2

#Loop over the conditions indicated in the CSV_annotation file and link them to the plotting colors
for condition,colour,colour2 in zip(conditions,colors, colors2):
    
    #Create an empty array, suitable for the number of conditions and valid worms that were analyzed 
    array_int=np.zeros((np.size(np.where((CSV_annotation.condition==condition) &(CSV_annotation.quality==1))), int(max(CSV_quantification.Frame))))
    counter=0

    #Loop over worms and select the ones, having the correct condition and representing a valid chamber/worm --> Extract the time of hatching and escaping of these worms
    for i,position in enumerate(CSV_annotation.position): #loops over all the worms
        if (CSV_annotation.condition[i]==condition): #take the worm out with the specific condition
            if CSV_annotation.quality[i]==1: #checks if the worm is good
                hatch=CSV_annotation.hatch[i] #find the hatch time of this worm
                escape=CSV_annotation.escape[i] #finds the escape time of this worm
                
                #Access the GFP quantification values from the hatching time until the time the worm escapes from the valid worms, defined before
                index=np.where((CSV_quantification.Position==position) & (CSV_quantification.Frame >= hatch) & (CSV_quantification.Frame < escape)) #find the indices of this specific worm that has timepoints in between hatch and end
                time_forworm=CSV_quantification.Frame[np.squeeze(index)] #gets out the valid timepoints for this worm
                quantification_forworm=CSV_quantification.final_intensity[np.squeeze(index)] #gets out the valit quantifications for this worm
                
                
                listtime.extend(time_forworm*min_image_taken/60)
                listGFP.extend(quantification_forworm)
                listcondition.extend(np.repeat(condition,time_forworm.size,axis=0))
                
        
                #Loop over the time points and corresponding GFP quantification measured from hatching until the worm escapes and add the values to the previously initialized array
                for (timepoint,quant) in zip(time_forworm, quantification_forworm):
                     array_int[counter,int(timepoint-hatch)]=quant  #in the row a worm, and every colum corresponding timepoint the quantification
                     
                counter=counter+1 #goes to the next valid worm
        
   

    #Replace all remainig 0 in the array by nan and convert the measured time points into hours
    array_int[array_int == 0] = 'nan'
    time=np.arange(0,np.size(array_int,1)*min_image_taken/60,min_image_taken/60)
    mean=np.nanmean(array_int,0)
    std=np.nanstd(array_int,0)
    
     
    
    #Plot the individual GFP quantifications of the selected worms starting from hatch until they escape
    for wormsel in np.arange(np.size(array_int,0)):
      # if ((wormsel==wormtoplot) & (condition=='control')):
           plt.plot(time,array_int[wormsel,:],color = colour2 ,linewidth=.6, alpha=0.1)
        
    #Plot the mean GFP intensity of all selected worms from one condition
    plt.plot(time,mean,color = colour, linewidth = 2,label=condition)   
        
        

#Define title, legens and other parameters of the plot
plt.title('single worm imaging quantification for every single worm')
plt.legend()
plt.xlabel('Time of larval development (h)')
plt.ylabel('GFP intensity (a.u.)')
plt.grid()



#plot with dataframe!!!   
plt.figure(2)
sns.set_theme(style="darkgrid")
df=pd.DataFrame({'time after hatch (h)':listtime, 'GFP intensity (a.u.)':listGFP,'condition':listcondition}) 
sns.lineplot(x='time after hatch (h)',y='GFP intensity (a.u.)',hue='condition',data=df, palette=colors).set_title('GFP intensity over time')
  
 


