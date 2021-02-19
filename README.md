# SWI_analysis
project for improving analysis of SWI experiments on heterochronic genes

# Install in conda environment:

```
conda create -n swi-analysis python=3.6
conda activate swi-analysis
pip install -e .
```

A note on this package, it seems that the version of scikit-image might be 0.15 instead of 0.17.

## Running the task

```
luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath /some/folder/with/images/ --channel-GFP pattern488 --channel-mcherry pattern566 --local-scheduler
```


# Current codes in the GitHub
- <code>swi_analysis_mcherry.py</code> :runs the script
- <code>marit_functions.py</code>  : contains all the functions that are needed to run the script
- <code>plot_complete_dataset.py</code> : plots the results

# Running the code from the (xenon7) server:
- Connect to pulse secure if you are not at the FMI
- Open the terminal and type: <code>ssh <i>username</i>@xenon7.fmi.ch</code>
- Create a new screen: <code>screen</code>
- Activate the conda environment: <code>conda activate <i>SWI</i></code>
  - To create a new conda environment: <code>conda create --name <i>SWI</i> python=3.6</code>
  - To install each package: <code>conda install <i>numpy=1.19</i></code>
- To find and set the working directory type <code>ls</code> and <code>cd /</code>
- Adjust in the python script the pathname where the images are located:
  > dirpath = "\tungstenfs\scratch\ggrossha\Lucas\Live_Imaging\LJ_Lin28_Project_210116
- To run the code of interest: <code>python <i>name of the file.py</i></code>
- To detach from screen by pressing: <kbd>ctrl</kbd>+<kbd>a</kbd>,  <kbd>ctrl</kbd>+<kbd>d</kbd>
- you can now close the terminal and the script will run. You can attach the screen again by first seeing the number of the screen by typing: <code>screen -ls</code> <br>
The output should look like:
  > <code>There are screens on:</code> <br>
  <code>  &nbsp; &nbsp;   <b>97684</b>.pts-9.xenon7   (Detached)</code><br>
  <code>  &nbsp; &nbsp;  <b>118940</b>.pts-1.xenon7   (Detached)</code>
     
- To reattach the screen: <code>screen -r <i>number</i></code> (i.e.: 97684)
  - The selected screen will now be displayed. <br>
  - **IMPORTANT**: Typing <code>screen</code> will make a new screen.
- After the job is finished, kill the screen by typing: <code>exit</code>
- If you want to stop the code, you must click <kbd>ctrl</kbd>+<kbd>c</kbd>
  
# git commands and links for help with git
First time: clone the repository from the website with git clone link
first check which branch you are, and the status of the branch: git status
If you want to download the recent github version of the code: git fetch
Check which branches are there:  
- git branch --list
create branch and go to the branch : 
- git branch <name_branch> -> git checkout <name_branch> 
or 
- git checkout -b <name_branch>
work on the code in the branch. When you are done with working, save the file, go back to the terminal and check the status of the folder.
You should see your new file in red. Add the changes by typing:
- git add <name_of_the_file>
If you want to commit it to the folder:
- git commit -m "commit message"
to upload the newest version of your branch to gitub : git push 
delete your branch git branch -d 




https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet

https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
