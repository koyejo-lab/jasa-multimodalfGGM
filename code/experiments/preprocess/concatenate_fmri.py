from scipy.io import loadmat, savemat
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, default="../data/fmri/", help="data path")
parser.add_argument('--ssid', type=int, default=1, help="session id")
args = parser.parse_args()

subjects=[
    #'T02_02',
    #'T03_03',
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
    #'T28_28',
    'T29_29'
]

path = args.path
data_list = list()
for s in subjects:
    data = loadmat("{}{}_run_{}.mat".format(path,s, args.ssid))

    data_list.append(data['ts'].squeeze()[:,0:150])

data_arr = np.array(data_list)
print(data_arr.shape)
data_dict = {'ts':data_arr, 'N':data_arr.shape[0], 'p':data_arr.shape[1], 't':data_arr.shape[2]}

savemat(path+'fmri_concate_run{}_reduced_movie.mat'.format(args.ssid), data_dict)