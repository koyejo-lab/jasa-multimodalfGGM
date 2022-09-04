Please replace the contents of this file with relevant instructions for your repository or remove this file entirely.

This directory would generally contain source code files that contain the core code to implement the method and various utility/auxiliary functions.

Scripts/code files that execute the overall workflow to carry out an analysis and generate results for the manuscript might be placed in the main directory.

This repository implements the algorithm of multimodal functional graphical model. The program is written in both python and R langauge. 

### Setup
Build an Anaconda environment and install python and R packages 
```
#create environment from environments.yml file
conda env create --name mmfggm --file=environments.yml

#activate mmfggm
conda activate mmfggm

```
Install the R packages
```
Rscript install_packages.R
```
if there is dependency error, try to update the R version by running `conda update --all`

### Example

#### Step-by-step tutorial


#### Synthesize Toy Data
You can generate the synthetic data by either running the scripts in `./data` or download the data from

#### Run Algorithm


### Experiment
