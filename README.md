# SWI_analysis
project for improving analysis of SWI experiments on heterochronic genes

# packages to install in conda environment:
python=3.6
numpy=1.19
matplotlib=3.3
scikit-image=0.17
scipy=1.5 
seaborn = 0.11
pandas=1.1

# main scripts present:
- swi_analysis_mcherry.py :runs the script
- marit_functions.py  : contains all the functions that are needed to run the script
- plot_complete_dataset : plots the results

# How to run things from the server:
log in to the server by:
- connecting to pulse secure if you are not at the FMI
- go to terminal, and connect with xenon7 by typing in the terminal: ssh  <code><i>username</i>@xenon7.fmi.ch</code>
- create a screen by typing: screen in the terimal
- activate/install correct conda environment
- go with ls and cd to the correct folder where the python script is located
- run the python script by typing: <code> python <i>name of the file.py</i></code>
- disattach from screen by pressing: Contr+<kbd>a</kbd>,d
- you can now close the terminal and the script will run. You can attach the screen again by first seeing the number of the screen by typing: screen -ls
- attach the screen by typing: screen -r <number>
- when the job is finished, close the screen by typing: exit
  
# help with git
https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
