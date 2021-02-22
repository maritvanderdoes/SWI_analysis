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
- First time: clone the repository from the website with   <code> git clone (link of the repository) </code>
- Next times:
 check which branch you are, and the status of the branch: <code>git status </code>
- Download the recent github version of the code: <code> git fetch </code> or <code> git pull </code>
- Check which branches are there:  <code> git branch --list </code>
- create branch and go there <code>git branch (name_branch) -> git checkout (name_branch)</code>  or  <code>git checkout -b (name_branch)</code>
- When you are done with working on the branch, save the file, go back to the terminal and check the status of the folder <code>git status </code>
- Adjusted files are in red. Add the changes:  <code> git add (name_of_the_file) </code>
- When new feature is implemented, commit the changes with a message what is improved: <code> git commit -m "commit message" </code>
- Upload the newest version of your branch to gitub : <code> git push origin (namebranch) </code>
- When you are finished with the issue, merge your branch to the main branch on github by doing a pull request
- Delete your branch  <code> git branch -d (branchname) </code>



tutorial:
https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet

list with commands
https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html

https://gist.github.com/rxaviers/7360908
