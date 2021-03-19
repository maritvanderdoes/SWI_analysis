# 🐛 SWI_analysis package description
## General description
Project to analyse single worm imaging (SWI) experiments, focused for endogenous levels of proteins. The package allows for the quantification of a fluorescent report by means of segmenting of worms using a paralell fluorescence channel.

## Package notes
Current version of the code is 0.2. A detail description of the codes [can be found here](docs/package_notes.md)

## Imaging notes
The code has been developed with a certain properties in mind and in a specific way. Information about these [can be found here](docs/imaging_noes.md)

# 🖥️ User's guide
## 1️⃣ Running the code for the first time
**If this is the first time using the package, [please follow this link](docs/first_time.md).**

## 📆 Running the code following times
**Running the code in the xenon7**
If the code is being run locally, then skip this step. Otherwise, type in the terminal
```
ssh <i>username</i>@xenon7.fmi.ch
cd /
cd <i>working_directory_in_which_SWI_analysis_is</i>
```
To open a new screen, type
```
screen
```
For more information on the usage of the Xenon7 server and screens, [please follow this link](docs/xenon7.md).

**Running the code locally**
First activate the environment and install the package. To do so, open the terminal and type:
```
conda activate swi-analysis
pip install -e .
```
Then, to run the task, type:
```
luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath /some/folder/with/images/ --channel-GFP pattern488 --channel-mcherry pattern566 --local-scheduler
```
The output of is given by <code>results.csv</code>. These results can be analysed using the package.

Other functions beyond <code>swi_analysis_mCherry</code> can be run, but these codes are more advanced. For that, check the [package notes](docs/package_notes.md)

## 📈 Plotting the data
The file <code>results.csv</code> can be analysed by using the function <code>plot_complete_dataset.py</code>.


# 🤔 Other (useful) information
[Here is a list of useful commands for GitHub](docs/github_usage.md)
