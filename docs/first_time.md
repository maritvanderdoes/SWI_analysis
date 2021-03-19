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
3. Setup your working directory by typing <code>cd</code> and the <code><i>address</i></code>. If you want to change to higher folders, type <code>cd ..</code>. An example
  ```
  cd ..
  cd ..
  cd ..  # These returns to the 01 folder. From here you can go to other folders
  cd tungstenfs/nobackup/ggrossha/moraluca/SWI_analysis/
  ```
4. 

- Connect to pulse secure if you are not at the FMI
- Open the terminal and type: <code>ssh <i>username</i>@xenon7.fmi.ch</code>
- Create a new screen: <code>screen</code>

For further information [please follow this link](docs/xenon7.md). 

## Install in conda environment:

```
conda create -n swi-analysis python=3.6
```

## Activating conda environment:

```
conda activate swi-analysis
pip install -e .
```

## Cloning the repository
To clone this repository, type 
```
git clone https://github.com/maritvanderdoes/SWI_analysis.git 
```
