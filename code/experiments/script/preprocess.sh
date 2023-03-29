#!/bin/bash

python ../preprocess/parse_eeg.py
python ../preprocess/parse_fmri.py

for i in 1 2 3;
do
    python ../preprocess/concatenate_eeg.py --ssid=$i
    python ../preprocess/concatenate_fmri.py --ssid=$i
done
python ../preprocess/concatenate_run.py

