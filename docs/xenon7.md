# Running the code in the xenon7 server
## Connecting to 
- Connect to pulse secure if you are not at the FMI
- Open the terminal and type: <code>ssh <i>username</i>@xenon7.fmi.ch</code>

## Working with screens
- Create a new screen: <code>screen</code>
- To find and set the working directory type <code>ls</code> and <code>cd /</code>
- To run the code of interest: <code>python <i>name of the file.py</i></code>
- To detach from screen by pressing: <kbd>ctrl</kbd>+<kbd>a</kbd>,  <kbd>ctrl</kbd>+<kbd>d</kbd>. This allows to safely close the window where the code is running.
- To see which screens are available, type <code>screen -ls</code> <br>
The output should look like:
  ```
  There are screens on: 
           97684.pts-9.xenon7   (Detached)
           118940.pts-1.xenon7   (Detached)
  ```
- To reattach the screen: <code>screen -r <i>number</i></code> (i.e.: 97684)
  - The selected screen will now be displayed. <br>
  - **IMPORTANT**: Typing <code>screen</code> will make a new screen.
- After the job is finished, kill the screen by typing: <code>exit</code>
- If you want to stop the code, you must click <kbd>ctrl</kbd>+<kbd>c</kbd>
  
## Working with Python 
- To run the code of interest: <code>python <i>name of the file.py</i></code>
- Activate the conda environment: <code>conda activate <i>SWI</i></code>
  - To create a new conda environment: <code>conda create --name <i>SWI</i> python=3.6</code>
  - To install each package: <code>conda install <i>numpy=1.19</i></code>
- Running a task with luigi:
  ```
  luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath /some/folder/with/images/ --channel-GFP pattern488 --channel-mcherry pattern566 --local-scheduler
  ```
