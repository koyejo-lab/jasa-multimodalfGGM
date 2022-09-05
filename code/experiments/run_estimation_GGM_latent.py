import numpy as np


import sys
import csv
import pickle

sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from evaluate import sparse_list, remove_tilde, construct_graph


from inverse_covariance import QuicGraphicalLasso
from numpy.linalg import slogdet,det
#load regression result
data = np.load('./result_experiment/experiment_result_movie_all_first23.npy', allow_pickle=True)

A1 = data.item()['Am'][0]
print(A1.shape)
A2 = data.item()['Am'][1]
print(A2.shape)
k = A1.shape[0]
tildeB = data.item()['tildeB']

graph = data.item()['graph'] #list

 


#compute AND operation
B = remove_tilde(tildeB)

print(np.max(np.abs(B)))
slist = sparse_list(B, 1e-3) 
adj = construct_graph(slist, AND=True)
for i in range(adj.shape[0]):
    adj[i,i] = 1

p = adj.shape[0]
k = A1.shape[0]
supp = np.zeros((p*k, p*k)) #pk x pk
##consturct support:
for i in range(p):
    for j in range(p):
        if adj[i,j] == 1:
            subm = np.zeros((k,k))
            sub_sup = tildeB[i,:,j*k:(j+1)*k]
            idx = np.where(sub_sup != 0)
            if idx[0].size>0:
                subm[:,:] = 1
            supp[i*k:(i+1)*k, j*k:(j+1)*k] = subm

plt.imshow(supp)
plt.savefig('supp_movie.png')
print(supp)

print(np.max(supp))
lam = (1-supp)*1e9
id = np.where(lam==0)
lam[id] = .1
#transform data to the latent space


filename= "fmri_eeg_score_movie_all.Rdata"
path = '/home/kt14/workbench/multimodal_functional_ggm/experiments/experiment_data'


name = robjects.r['load'](path+"/"+filename)
Robj = robjects.r[name[0]]
Rdict = dict(zip( Robj.names, map(list, list(Robj))))
fmri = np.array(robjects.r['data'][0])[0:23]
eeg = np.array(robjects.r['data'][1])[0:23]

tfmri = np.zeros((fmri.shape[0],fmri.shape[1]*k))
teeg = np.zeros((eeg.shape[0], eeg.shape[1]*k))

for n in range(fmri.shape[0]):
    for i in range(fmri.shape[1]):
        tfmri[n, i*k:(i+1)*k] = np.einsum('ij,j->i', A1, fmri[n,i,:])

for n in range(eeg.shape[0]):
    for i in range(eeg.shape[1]):
        teeg[n, i*k:(i+1)*k] = np.einsum('ij,j->i', A2, eeg[n,i,:])

latent_data = np.vstack((tfmri,teeg))
latent_mu = np.mean(latent_data, axis=0)

#center the data
latent_data = latent_data - latent_mu


#testdata
fmri = np.array(robjects.r['data'][0])[23:46]
eeg = np.array(robjects.r['data'][1])[23:46]

tfmri = np.zeros((fmri.shape[0],fmri.shape[1]*k))
teeg = np.zeros((eeg.shape[0], eeg.shape[1]*k))

for n in range(fmri.shape[0]):
    for i in range(fmri.shape[1]):
        tfmri[n, i*k:(i+1)*k] = np.einsum('ij,j->i', A1, fmri[n,i,:])

for n in range(eeg.shape[0]):
    for i in range(eeg.shape[1]):
        teeg[n, i*k:(i+1)*k] = np.einsum('ij,j->i', A2, eeg[n,i,:])

test_data = np.vstack((tfmri,teeg))
#center the data
test_data = test_data - latent_mu
#initiate estimator







def compute_loglikelihood(X, prec):
    N = X.shape[0]
    pk = X.shape[1]
    print(slogdet)
    _, log_det_prec = slogdet(prec)
    term2 = N*(log_det_prec - pk*np.log(2*np.pi))
    
    term1 = 0
    for n in range(N):
        temp1 = np.einsum('ij,j->i', prec, X[n])
        temp2 = np.dot(temp1, X[n])
        term1 += -0.5*temp2
    return term1 + term2

n_fold = 5
acc_lln = []
for i in range(n_fold):
    select_t = np.random.choice(46,46, replace=True)
    
    estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=lam, verbose=1)
    estimator.fit(latent_data[select_t])
    prec = estimator.precision_
    print(prec)
    print(slogdet(prec))
    np.save('result_run1_movie_first23_1e_3_3_fold{}x.npy'.format(i),prec)

    print(np.max(prec))

    insample = compute_loglikelihood(latent_data[select_t], prec)
    print('in-sample loglikelihood:', insample)
    outsample = compute_loglikelihood(test_data, prec)
    print('out-of-sample loglikelihood', outsample)
    acc_lln.append(insample)
print('average in-sample loglikelihood', np.mean(acc_lln), np.std(acc_lln))
#fit model 


#plotting.plot_matrix(adj, figure=(10, 8), labels=desikan_altas_name, grid=True,
#                     vmax=0.8, vmin=-0.8, reorder=False, colorbar=False)
#plt.tight_layout()
#plt.savefig('adj_latent_and_04.png')





