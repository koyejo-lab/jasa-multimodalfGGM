### Experiments


#### Step 1: Preprocessing data
Please check the README in the folder `./preprocess/` to prepare the data

#### Step 2: Edge estimation
Run the following command to estimate the edge set of the latent graph
```
python run_experiments.py
```
Run `run_experiments_cv.py` to select the hyperparameters with 5-fold cross validation

#### Step 3: Graph estimation

Execute the following program to obtain the precision matrix, in-sample logliklihood, out-of-sample loglikelihood
```
python run_estimation_GGM_eeg.py     #for eeg data
python run_estimation_GGM_fmri.py    #for fmri data
python run_estimation_GGm_latent.py  #for estimating the latent graph.
```

### Step 4: Vieuslization 
run the notebook `visualize_precision.ipynb` to visualize the result

#### Neighborhood regression
To estimate the graph using neighborhood regression method, run the following command

```
python run_eeg.py     #for eeg data
python run_fmri.py    #for fmri data
```

