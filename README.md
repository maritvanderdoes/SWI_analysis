# SWI_analysis
project for improving analysis of SWI experiments on heterochronic genes

# packages to install in conda environment:
>> python=3.6
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
- go to terminal, and connect with xenon7 by typing in the terminal: <code>ssh <i>username</i>@xenon7.fmi.ch</code>
- create a screen by typing: screen in the terimal
- activate/install correct conda environment: <code> conda activate <i>SWI</i> </code>
- go with <code>ls</code> and <code>cd</code> to the correct folder where the python script is located
- run the python script by typing: <code> python <i>name of the file.py</i></code>
- detach from screen by pressing: <kbd>ctrl</kbd>+<kbd>a</kbd>,  <kbd>ctrl</kbd>+<kbd>d</kbd>
- you can now close the terminal and the script will run. You can attach the screen again by first seeing the number of the screen by typing: <code> screen -ls </code>
- attach the screen by typing: <code>screen -r <i>number</i></code>
  - The screens will be displayed. **IMPORTANT** Typing <code>screen</code> will make a new screen.
- when the job is finished, close the screen by typing: <code>exit</code>
- If you want to stop the code, you must click <kbd>ctrl</kbd>+<kbd>c</kbd>
  
# help with git
https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
