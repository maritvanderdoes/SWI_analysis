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
6. To activate the environment and install the package, type
  ```
  conda activate SWI
  pip install -e .
  ```

- Connect to pulse secure if you are not at the FMI
- Open the terminal and type: <code>ssh <i>username</i>@xenon7.fmi.ch</code>
- Create a new screen: <code>screen</code>

For further information [please follow this link](docs/xenon7.md). 
