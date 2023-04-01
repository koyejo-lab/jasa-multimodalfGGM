#Instruction

Data generation files are under `../code/synth_data`. Please see the instructions in `../code/synth_data/README.md` for more details.

Alternatively, pre-generated data can be downloaded [here](https://drive.google.com/drive/folders/1EbHl0Q2oE_ME3WjLWINdSlg_M9VJM0Qy?usp=share_link)



# Data format

This notebook describes the simulated data dictionary. 

## Filename
We desribe the naming of the files.
Each simulated data contains two files: one file begin with `graph*.RData` and the other begins with `data*.RData`.
The `graph*.RData` contains the ground truth and `data*.RData` contains the simulated data. A typical filename would look like
```
graph_{#graph type}_p{#dimension}_N{#sample size}{#noise model}.RData
data_{#graph type}_p{#dimension}_N{#sample size}{#noise model}.RData
```

## Data type

### The data file
The  data in `data*.RData` is a dictionary object `f`. After loading the data using function `load_Rdata()` (see the function description under `code/src/Utility/Python2R`), each data modal is acessed in the following format
```
data_m = f['data'][the modality index] #ndarray.
```
`data_m` is an ndarray with three dimension: the first dimension is the subject index $n$, the second dimension is the node index $p$, the third dimension is the basis function indez $k_m$. The initial ${\bf A}^{m (0)}$ can be acessed by `f['A'][m]`.


The  data in `graph*.RData` is a  dictionary object `g`. After loading the data using function `load_Rparams()` (see the function description under `code/src/Utility/Python2R`)
`g['Am']` is a list of ${\bf A}^{m\star}$ and `g['tildeB']` is an ndarray of dimension $(p,k,k*p)$ where `g['tildeB'][i]` is $\tilde{\bf B}_i^\star$.