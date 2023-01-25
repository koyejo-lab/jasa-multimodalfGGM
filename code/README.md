This repository implements the algorithm of multimodal functional graphical model. The program is written in both python and R langauge. 

# Setup
Build an Anaconda environment and install Python and R packages. The data generation processes and the initialization method are written in R. The iterative refinement procedure is implemented in Python.


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

# Simulation

## Step-by-step tutorial
Ｓtep-by-step estimation procedures of the implmentations are under the directory `./tests/notebook/`

## Synthesize Toy Data
You can generate the synthetic data by either running the scripts in `./data/`. The `datageneration_*.R` files include steps of the initialization of $A^m$.

For detailed instruction, please check the `README` under `./data/`

## Run Algorithm

### Run main experiments 
```
python ./tests/run_model3.py --init [initialization procedure] --p [dimension]  --N [sample size]  --type [graph type] --filepath [path to synthetic data]  --savepath [path to store result] --noise [noise model] --thre [threshold] --lr_initb [learning rate for initializing B] --lr_a [learning rate for A] --lr_b [learning rate for B]
```
The simulation configurations are stored in the shell scripts under the directory `./tests/scripts/` 

### Run cross-validation 
```
python ./tests/run_cv_model3.py --init [initialization procedure] --p [dimension]  --N [sample size]  --type [graph type] --filepath [path to synthetic data]  --savepath [path to store result] --noise [noise model] --thre [threshold] 
```

### Visualize the results

To visualize the results of comparison of different methods, run `./tests/notebook/plot_Comparison.ipynb`. 

To visualize the results of different sample size, run `./tests/notebook/plot_SampleComplexity.ipynb`.


# Experiment on Real Data
For the details of experiments on concurrent fMRI-EEG measurements, please check the `README` under the directory `./experiments/`

