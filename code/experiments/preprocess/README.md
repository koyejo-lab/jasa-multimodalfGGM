### Preprocess Pipeline


#### Step 1: standard preprocessing
```
python parse_eeg.py
python parse_fmri.py
```

#### Step 2: Concatenate data across testing subjects (N)

```
python concatenate_eeg.py
python concatenate_fmri.py
```

#### Step 3:  Concatenate data from different sessions

```
python concatenate_run.py
```

#### Step 4: Compute the functional score

```
Rscript extract.score.R
```

To select the optimal basis and number of functions, please check `compute.score.R`
#### Step 5: Initiate the A operator

```
Rscript compute.correlation.R
```