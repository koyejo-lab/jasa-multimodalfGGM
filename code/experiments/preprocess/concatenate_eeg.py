from scipy.io import loadmat, savemat
import numpy as np
from scipy.signal import butter, lfilter, filtfilt

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y

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

path = "../resting_data/eeg/"
data_list = list()
delta_list = list()
theta_list = list()
alpha_list = list()
beta_list = list()
gamma_list = list()

total = 150000 - (3320*2)
total = 75000 - 3320
idx = np.arange(0, total, 35)
print(idx.shape)
sampling_rate = 150000 / 600

##############
# delta:1-4
# theta:4-7
# alpha:8-15
# beta:16-31
# gamma:>32
#############
#start 3320, remove head
#end 150000-3320
print(idx)
for s in subjects:
    data = loadmat("{}{}_run_3.mat".format(path,s))

    #s_data = data['ts'].squeeze()[:,3320:(150000-3320)]
    s_data = data['ts'].squeeze()[:,3320:75000]  # remove head and tail
    delta_data = np.zeros(s_data.shape)
    theta_data = np.zeros(s_data.shape)
    alpha_data = np.zeros(s_data.shape)
    beta_data = np.zeros(s_data.shape)
    gamma_data = np.zeros(s_data.shape)
    for i in range(s_data.shape[0]):    
        delta_data[i,:] = butter_bandpass_filter(s_data[i,:], 1,  4,  sampling_rate)
        theta_data[i,:] = butter_bandpass_filter(s_data[i,:], 4,  7,  sampling_rate)
        alpha_data[i,:] = butter_bandpass_filter(s_data[i,:], 8,  15, sampling_rate)
        beta_data[i,:]  = butter_bandpass_filter(s_data[i,:], 16, 31, sampling_rate)
        gamma_data[i,:] = butter_highpass_filter(s_data[i,:], 32,     sampling_rate)

    #s_data = s_data[:,idx]
    
    data_list.append(s_data[:,idx])
    delta_list.append(delta_data[:,idx])
    theta_list.append(theta_data[:,idx])
    alpha_list.append(alpha_data[:,idx])
    beta_list.append(beta_data[:,idx])
    gamma_list.append(gamma_data[:,idx])

data_arr = np.array(data_list)

data_dict = {'ts': np.array(data_list),
             'delta': np.array(delta_list),
             'theta': np.array(theta_list),
             'alpha': np.array(alpha_list),
             'beta':  np.array(beta_list),
             'gamma': np.array(gamma_list),
             'N': data_arr.shape[0], 'p':data_arr.shape[1], 't':data_arr.shape[2]}


savemat(path+'eeg_concate_run3_down35_reduced_movie.mat', data_dict)