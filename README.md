# SWI_analysis package description
## General description
Project to analyse single worm imaging (SWI) experiments, focused for endogenous levels of proteins. The package allows for the quantification of a fluorescent report by means of semgnenting of worms using a paralell fluorescence channel.

## Package notes
Current version of the code is 0.2. A detail description of the codes [can be found here](docs/package_notes.md)

# User's guide
## Running the code for the first time
**If this is the first time using the pacakge, [please follow this link](docs/first_time.md)**

## Running the code following times
If the code is being run on the xenon7 server, [please follow this link](docs/xenon7.md)

```
conda activate swi-analysis
pip install -e .
```

A note on this package, it seems that the version of scikit-image might be 0.15 instead of 0.17.

## Running the task

```
luigi --module tasks.swi_analysis_mCherry SWIAnalysisTask --dirpath /some/folder/with/images/ --channel-GFP pattern488 --channel-mcherry pattern566 --local-scheduler
```


# Git information
[Here is a list of useful commands for GitHub](docs/github_usage.md)
