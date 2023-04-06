# Workflow


##### Figure 2/Table 3
* Data preparation
  * Data generation script: `./code/synth_data/run_dgp_N2.sh`
  * Data downlaod [link](https://drive.google.com/drive/folders/1V-UalL77ybJZ8yq7mzkb1o1EeVDLR5Mw?usp=share_link). Store the files under `data_batch_N2`
* Estimation 
  * Script file: `./code/tests/script/noise_model_2/*N100.sh`.
  Modify the script file to specify the conda environment, file path, and save path.
* Visualization
  * Result download [link](https://drive.google.com/drive/folders/1praDtxiYMUVCNnnp30SG3JE5REIVCz7d?usp=sharing). The directory `/proposed` contain the results of proposed method. The directory `/comparison/` constain the results of other comparison methods. 
  * Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  * Instruction to generate table:
  to print the AUC and AUC15, set `verbose=True`

##### Figure 3
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_sample.sh`
- Estimation 
  - Script file: `./code/tests/script/sample_sample/`.
  Modify the script file to specify the conda environment, file path, and save path.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1bWxoFzd9dAZn1x7VEulELsUIbO6cAEfk?usp=sharing). Please select to download the directory `./p50`, `./p100`, `./p150` 
  - Visualization notebook:`/code/notebook/plot_SampleComplexity.ipynb`

##### Figure 5/Table 2
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_N1.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1mTZW2goCCgp7ULKmd4QIxkv95maehpSI?usp=share_link). Store the files under `data_batch_N1`
- Estimation 
  - Script file: `./code/tests/script/noise_model_1/*N100.sh`.
  Modify the script file to specify the conda environment, file path, and save path.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1jE8L4OWDd2qURhMALCo9tjrnWQbIWP8u?usp=sharing)
  The directory `/proposed` contain the results of proposed method. The directory `/comparison/` constain the results of other comparison methods. 
  - Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  - Instruction to generate table:
  to print the AUC and AUC15, set `verbose=True`

##### Figure 6/Table 4
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_N2.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1V-UalL77ybJZ8yq7mzkb1o1EeVDLR5Mw?usp=share_link). Store the files under `data_batch_N2`
- Estimation 
  - Script file: `./code/tests/script/noise_model_2/`.
  Modify the script file to specify the conda environment, file path, and save path.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1V-UalL77ybJZ8yq7mzkb1o1EeVDLR5Mw?usp=share_link)
  The directory `/proposed` contain the results of proposed method. 
  - Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  - Instruction to generate table:
  to print the AUC and AUC15, set `verbose=True`

##### Figure 7/Table 5
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_N1.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1mTZW2goCCgp7ULKmd4QIxkv95maehpSI?usp=share_link). Store the files under `data_batch_N1`
- Estimation 
  - Script file: `./code/tests/script/noise_model_1/`.
  Modify the script file to specify the conda environment, file path, and save path.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1jE8L4OWDd2qURhMALCo9tjrnWQbIWP8u?usp=sharing)
  The directory `/proposed` contain the results of proposed method. 
  - Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  - Instruction to generate table:
  to print the AUC and AUC15, set `verbose=True`

##### Figure 8
- Data preparation
  - Data generation script: `/code/synth_data/run_dgp_k.sh`,
  - Data downlaod [link](https://drive.google.com/drive/folders/151_r-ttGr2uEepZX_jkUz9F3MFl1R0Po?usp=sharing)
- Estimation 
  - Script file: `./code/tests/script/sample_k/`.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1bWxoFzd9dAZn1x7VEulELsUIbO6cAEfk?usp=sharing) 
  - Visualization notebook: `/code/notebook/plot_SampleComplexity_2.ipynb`

##### Figure 9
- Estimation 
  - Run `/code/tests/notebook/plot_elbo.ipynb` and save the result
- Visualization
  - Visualization notebook: `/code/tests/notebook/plot_elbo2.ipynb`

##### Figure 10
- Estimation 
  - Run `/code/tests/notebook/plot_elbo.ipynb` and save the result
- Visualization
  - Visualization notebook: `/code/tests/notebook/plot_elbo2.ipynb`

##### Figure 11
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_N2.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1V-UalL77ybJZ8yq7mzkb1o1EeVDLR5Mw?usp=share_link). Store the files under `data_batch_N2`
- Estimation 
  - Script file: `./code/tests/script/noise_model_2/`.
  Modify the script file to specify the conda environment, file path, and save path
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1LjQwRMvyhTJv2HB1xa3-TsCTQCPDp46H?usp=sharing)
  - Visualization notebook: `/code/tests/notebook/plot_VariableSelection.ipynb`

##### Figure 12
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_sample.sh`
- Estimation 
  - Script file: `./code/tests/script/sample_alpha/`.
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1bWxoFzd9dAZn1x7VEulELsUIbO6cAEfk?usp=sharing) 
  - Visualization notebook: `/code/notebook/plot_SampleComplexity_2.ipynb`


##### Figure 13/Table 6
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_kmk.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1TYgPNnpni8UV95MzUXu-41vs2GeW-F3i?usp=share_link)
- Estimation 
  - Script file: `./code/tests/script/noise_model1_varykmk/`
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1ItaHAojkT2SrK2Ei2bGUU6QJ_WErvej9?usp=sharing)
  - Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  - Instruction to generate table: to print the AUC and AUC15, set `verbose=True`
  
##### Figure 14/Table 7
- Data preparation
  - Data generation script: `./code/synth_data/run_dgp_kmk.sh`
  - Data downlaod [link](https://drive.google.com/drive/folders/1TYgPNnpni8UV95MzUXu-41vs2GeW-F3i?usp=share_link)
- Estimation 
  - Script file: `./code/tests/script/noise_model1_varykmk/`
- Visualization
  - Result download [link](https://drive.google.com/drive/folders/1ItaHAojkT2SrK2Ei2bGUU6QJ_WErvej9?usp=sharing)
  - Visualization notebook: `/code/notebook/plot_Comparison.ipynb`
  - Instruction to generate table: to print the AUC and AUC15, set `verbose=True`