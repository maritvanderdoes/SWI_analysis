# First time using the SWI_analysis package
Here is a detalied tutorial con how to use the cod

First open the terminal (in Mac) or type <code>cmd</code> after clicking Windows button (in Windows).

## Running it on Xenon
If you are running the code in the xenon7, follow this section. Otherwise, jump to the next part.

0. If you are not at FMI, connect to the network by Pulse Secure.
1. Open the terminal and type:
  ```
  ssh <i>username</i>@xenon7.fmi.ch
  ```
2. The system will ask for your password. Type it and then push 
3. Setup your working directory by typing <code>cd</code> and the <code><i>address</i></code>. If you want to change to higher folders, type <code>cd /</code>. An example
  ```
  cd / # These returns to the 01 folder. From here you can go to other folders
  cd tungstenfs/nobackup/ggrossha/moraluca/
  ```
4. Once here, clone the repository. This should create a <code>SWI_analysis</code> folder.
  ```
  git clone https://github.com/maritvanderdoes/SWI_analysis.git 
  ```
5. After cloning the repository, (move to the new downloaded folder and create a virtual environment byy typing:
  ```
  cd SWI_analysis
  conda create -n SWI python=3.6
  ```
6. Create a new screen by typing
  ```
  screen
  ```
7. To activate the environment and install the package, type
  ```
  conda activate SWI
  pip install -e .
  ```
8. To execture the code type 
  ```
  luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath /tungstenfs/scratch/ggrossha/Lucas/Live_Imaging/LJ_Lin28_Project_210116 --outputpath /tungstenfs/nobackup/ggrossha/moraluca/lin28_new_210116 --channel-GFP w1Lucas-sim-488-561.stk --channel-mcherry w2Lucas-sim-561-488.stk --local-scheduler
  ```
9. If everything is working well, the filename should be displayed. You can detach the screen by pressing <kbd>ctrl</kbd>+<kbd>a</kbd>,  <kbd>ctrl</kbd>+<kbd>d</kbd>. Now you can safely close the window.

## Common issues
❓ **The code seems to run, but it does not find any image:** Either you have not put the path correctly (you might need the <code>/</code> at the beginning, or the channel extensions are not wwritten properly (due to compression, the extension can be renamed to <code>.tiff</code> from <code>.stk</code>).

❓ **I receive the error displayed below:** There is an issue with your <code>scikit-image</code> package and the image is not loaded properly (it lacks the z-dimension). You should try to install a different version of the package.
```
ERROR: [pid 49948] Worker Worker(salt=929306007, workers=1, host=f146l-f6ad55, username=moraluca, pid=49948) failed 
SWIAnalysisTask(dirpath=C:/Users/moraluca/Desktop/Lin28_test, outputpath=C:/Users/moraluca/Desktop/Lin28_test/Output, 
channel_GFP=w1Lucas-sim-488-561.stk, channel_mcherry=w2Lucas-sim-561-488.stk)
Traceback (most recent call last):
  File "c:\programdata\anaconda3\lib\site-packages\luigi\worker.py", line 191, in run
    new_deps = self._run_get_new_deps()
  File "c:\programdata\anaconda3\lib\site-packages\luigi\worker.py", line 133, in _run_get_new_deps
    task_gen = self.task.run()
  File "d:\github\swi_analysis\tasks\swi_analysis_mCherry_troubleshooting.py", line 45, in run
    binary_mask, additional_info = adaptive_masking(images_out[0])
  File "d:\github\swi_analysis\utils\core_utils.py", line 249, in adaptive_masking
    sorted_values = np.zeros([datdim[0],datdim[1]*datdim[2],3])
  ```

❓ **I do not have a results folder:** The job might:
- Not finished.
- Finish prematurely due to an error.
- The output folder migth be incorrect. Check if you can find a folder that should not be in the <code>SWI_analysis</code> folder.
- Server maintenance might have ended the job prematurely. Make sure not to queue jobs before a maintenance event. You can know these events through the fmi email. 

For further information [please follow this link](docs/xenon7.md). 
