from scipy.io import loadmat, savemat
from scipy.stats import zscore
from scipy.signal import detrend

import numpy as np
subjects=[
    'T02_02',
    'T03_03',
    'T04_04',
    'T05_05',
    'T06_06',
    'T07_07',
    'T08_08',
    'T09_09',
    'T10_10',
    'T11_11',
    'T12_12',
    'T13_13',
    'T14_14',
    'T15_15',
    'T16_16',
    'T17_18',
    'T19_07',
    'T20_11',
    'T21_14',
    'T22_16',
    'T23_19',
    'T24_22',
    'T25_23',
    'T27_26',
    'T28_28',
    'T29_29'
]

path="../raw_data/"
eeg_map= loadmat(path+"aparc_aseg_yeoRS7_68reg_eeg_nosubc.mat")
yeoOrder_eeg = eeg_map['yeoOrder_eeg']


with open(path+'desi_brainstorm2freesurfer.txt') as f:
    pd = f.readlines()[0].split()
desi_brainstorm2freesurfer = [int(i)-1 for i in pd]

savepath="../data/eeg/"


for s in subjects:
    
    if 'T02_02' in s:
        run_ind = [1, 3]
    elif 'T03_03' in s:
        run_ind = [1, 2]
    elif 'T28_28' in s:
        run_ind = [1, 2]
    else:
        run_ind = [1, 2, 3]

    
    for run in run_ind:
        data = loadmat("{}eeg_rest/{}/rest_{}.mat".format(path, s, run))
        ImageGridAmp = data['ImageGridAmp']
        #print(data.keys())
        ts = ImageGridAmp[:,0:150000]
        ts = ts[desi_brainstorm2freesurfer,:]
        ts = ts[(yeoOrder_eeg-1),:]

        #TR=600/150000;

        for seed in range(68):
            
            std = np.std(ts[seed,:])
            mean = np.mean(ts[seed,:])
            #print(std, mean)
            #print(max(max(abs((ts[seed,:]-mean)/std))))
            ts[seed,:]= detrend((ts[seed,:]-mean)/std)
            #print(max(max(abs(ts[seed,:]))))
        savedata = {'ts':ts}
        filename = "{}_run_{}.mat".format(s,run)
        savemat(savepath+filename, savedata)