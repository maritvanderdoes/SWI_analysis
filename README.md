# SWI_analysis
project for improving analysis of SWI experiments on heterochronic genes

# packages to install in conda environment:
> python=3.6 <br>
numpy=1.19 <br>
matplotlib=3.3 <br>
scikit-image=0.17 <br>
scipy=1.5  <br>
seaborn = 0.11 <br>
pandas=1.1 <br>

# main scripts present:
- <code>swi_analysis_mcherry.py</code> :runs the script
- <code>marit_functions.py</code>  : contains all the functions that are needed to run the script
- <code>plot_complete_dataset.py</code> : plots the results

# How to run things from the server:
log in to the server by:
- connecting to pulse secure if you are not at the FMI
- go to terminal, and connect with xenon7 by typing in the terminal: <code>ssh <i>username</i>@xenon7.fmi.ch</code>
- create a screen by typing in the terminal <code>screen</code>
- activate/install correct conda environment: <code>conda activate <i>SWI</i></code>
- go with <code>ls</code> and <code>cd</code> to the correct folder where the python script is located
- run the python script by typing: <code>python <i>name of the file.py</i></code>
- detach from screen by pressing: <kbd>ctrl</kbd>+<kbd>a</kbd>,  <kbd>ctrl</kbd>+<kbd>d</kbd>
- you can now close the terminal and the script will run. You can attach the screen again by first seeing the number of the screen by typing: <code>screen -ls</code>
The output should look like
> There are screens on:
     97684.pts-9.xenon7   (Detached)
     118940.pts-9.xenon7   (Detached)
     
- attach the screen by typing: <code>screen -r <i>number</i></code>
  - The screens will be displayed. <br>
  **IMPORTANT**: Typing <code>screen</code> will make a new screen.
- when the job is finished, close the screen by typing: <code>exit</code>
- If you want to stop the code, you must click <kbd>ctrl</kbd>+<kbd>c</kbd>
  
# help with git
https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
