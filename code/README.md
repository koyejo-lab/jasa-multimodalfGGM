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
Ｓtep-by-step estimation procedures of the implementations are under the directory `./tests/notebook/`

## Synthesize Toy Data
You can generate the synthetic data by either running the scripts in `./synth_data/`. The `datageneration_*.R` files include steps of the initialization of $A^m$.

For detailed instructions, please check the `README` under `./synth_data/`

## Run Algorithm

### Run main experiments 
To execute the program parallelly, use the following command. This code estimates the graph in a range of sparsity levels and is designed to run on a multi-core machine.
```
python ./tests/run_model3.py --init [initialization procedure] --p [dimension]  --N [sample size]  --type [graph type] --filepath [path to synthetic data]  --savepath [path to store result] --noise [noise model] --thre [threshold] --lr_initb [learning rate for initializing B] --lr_a [learning rate for A] --lr_b [learning rate for B]
```
The simulation configurations are stored in the shell scripts under the directory `./tests/scripts/`. The path to data and output might need to be modified. The data generation process is described in the previous section. Alternatively, the data can be downloaded [here](https://drive.google.com/drive/folders/1EbHl0Q2oE_ME3WjLWINdSlg_M9VJM0Qy?usp=share_link). 

To execute the program with a specific sparsity level, use the following command.

```
python ./tests/run_single_model3.py --init [initialization procedure] --p [dimension]  --N [sample size]  --type [graph type] --filepath [path to synthetic data]  --savepath [path to store result] --noise [noise model] --thre [threshold] --lr_initb [learning rate for initializing B] --lr_a [learning rate for A] --lr_b [learning rate for B] --s [sparsity level]
```

The estimated model will be stored in a `.npy` file and the estimation result will be stored in a `.csv` file.

### Run cross-validation 
```
python ./tests/run_cv_model3.py --init [initialization procedure] --p [dimension]  --N [sample size]  --type [graph type] --filepath [path to synthetic data]  --savepath [path to store result] --noise [noise model] --thre [threshold] 
```

### Visualize the results

Download the estimation results [here](https://drive.google.com/drive/folders/1EbHl0Q2oE_ME3WjLWINdSlg_M9VJM0Qy?usp=share_link):


To visualize the results of the comparison of different methods, run `./tests/notebook/plot_Comparison.ipynb`. 

To visualize the results of different sample sizes, run `./tests/notebook/plot_SampleComplexity.ipynb`, `./tests/notebook/plot_SampleComplexity2.ipynb`

To visualize the results of variable selection, run `./tests/notebook/plot_elbo.ipynb` and `./tests/notebook/plot_VariableSelection.ipynb`


# Experiment on Real Data
For the details of experiments on concurrent fMRI-EEG measurements, please check the `README` under the directory `./experiments/`. The experiment data is not released in public but is available upon request.

