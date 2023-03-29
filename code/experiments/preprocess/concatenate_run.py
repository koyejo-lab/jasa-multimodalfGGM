from scipy.io import loadmat, savemat
import numpy as np

path = '../data/eeg/'
mat1 = loadmat(path + 'eeg_concate_run1_down70_reduced_movie.mat')
mat2 = loadmat(path + 'eeg_concate_run2_down70_reduced_movie.mat')
mat3 = loadmat(path + 'eeg_concate_run3_down70_reduced_movie.mat')

ts    = np.concatenate((mat1['ts'], mat2['ts'], mat3['ts']), axis=0)
delta = np.concatenate((mat1['delta'], mat2['delta'], mat3['delta']), axis=0)
theta = np.concatenate((mat1['theta'], mat2['theta'], mat3['theta']), axis=0)
alpha = np.concatenate((mat1['alpha'], mat2['alpha'], mat3['alpha']), axis=0)
beta  = np.concatenate((mat1['beta'],  mat2['beta'],  mat3['beta']), axis=0)
gamma = np.concatenate((mat1['gamma'], mat2['gamma'], mat3['gamma']), axis=0)
N = ts.shape[0]
p = ts.shape[1]
t = ts.shape[2]

data_dict = {
    'ts': ts,
    'delta': delta,
    'theta': theta,
    'alpha': alpha,
    'beta':  beta,
    'gamma': gamma,
    'N': N,
    'p': p,
    't': t
}


savemat(path+'eeg_concate_all_down70_reduced_movie.mat', data_dict)



path = '../data/fmri/'
mat1 = loadmat(path + 'fmri_concate_run1_reduced_movie.mat')
mat2 = loadmat(path + 'fmri_concate_run2_reduced_movie.mat')
mat3 = loadmat(path + 'fmri_concate_run3_reduced_movie.mat')

ts    = np.concatenate((mat1['ts'], mat2['ts'], mat3['ts']), axis=0)
N = ts.shape[0]
p = ts.shape[1]
t = ts.shape[2]

data_dict = {
    'ts': ts,
    'N': N,
    'p': p,
    't': t
}


savemat(path+'fmri_concate_all_reduced_movie.mat', data_dict)