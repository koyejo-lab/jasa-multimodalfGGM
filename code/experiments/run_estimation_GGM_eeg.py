import numpy as np


import sys
import csv
import pickle

sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")
import rpy2.robjects as robjects
import matplotlib.pyplot as plt
from evaluate import sparse_list, remove_tilde, construct_graph
from functools import partial

from itertools import repeat
from inverse_covariance import QuicGraphicalLasso
from numpy.linalg import slogdet,det
from multiprocessing import Pool, cpu_count
import time

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

def compute_BIC(X, prec):
    

    term1 = -2*compute_loglikelihood(X,prec)
    #compute the sparsity pattern
    s = np.where(np.abs(prec) < 1e-5)[0].size

    N = X.shape[0]

    return term1 + s*np.log10(N)

#transform data

filename= "fmri_eeg_score_movie_all.Rdata"
path = '/home/kt14/workbench/multimodal_functional_ggm/experiments/experiment_data'


name = robjects.r['load'](path+"/"+filename)
Robj = robjects.r[name[0]]
Rdict = dict(zip( Robj.names, map(list, list(Robj))))
fmri = np.array(robjects.r['data'][0])[0:23]
eeg = np.array(robjects.r['data'][1])[0:23]
k1 = fmri.shape[2]
k2 = eeg.shape[2]

tr_eeg = np.zeros((eeg.shape[0], eeg.shape[1]*k2))



for n in range(eeg.shape[0]):
    for i in range(eeg.shape[1]):
        tr_eeg[n, i*k2:(i+1)*k2] =  eeg[n,i,:]

eeg_mean = np.mean(tr_eeg, axis = 0)
tr_eeg = tr_eeg - eeg_mean

#testdata
fmri = np.array(robjects.r['data'][0])[23:46]
eeg = np.array(robjects.r['data'][1])[23:46]

te_eeg = np.zeros((eeg.shape[0], eeg.shape[1]*k2))


for n in range(eeg.shape[0]):
    for i in range(eeg.shape[1]):
        te_eeg[n, i*k2:(i+1)*k2] =  eeg[n,i,:]


te_eeg  -= eeg_mean




##########################
#load regression result  #
#construct support set   #
##########################

data = np.load('./result_experiment/experiment_result_eeg_movie_all_first23.npy', allow_pickle=True)

tildeB = data.item()['tildeB']

graph = data.item()['graph'] #list

 
#compute AND operation
B = remove_tilde(tildeB)

print(np.max(np.abs(B)))
slist = sparse_list(B, 1e-1) 
adj = construct_graph(slist, AND=True)
for i in range(adj.shape[0]):
    adj[i,i] = 1

p = adj.shape[0]
eeg_supp = np.zeros((p*k2, p*k2)) #pk x pk
fused_supp = np.zeros((p*k2,p*k2))
##consturct support:
for i in range(p):
    for j in range(p):
        if adj[i,j] == 1:
            subm = np.zeros((k2,k2))
            sub_sup = tildeB[i,:,j*k2:(j+1)*k2]
            idx = np.where(sub_sup != 0)
            if idx[0].size>0:
                subm[:,:] = 1
            eeg_supp[i*k2:(i+1)*k2, j*k2:(j+1)*k2] = subm
            fused_supp[i*k2:(i+1)*k2, j*k2:(j+1)*k2] = subm

plt.imshow(eeg_supp)
plt.savefig('supp_movie_fmri.png')
print(np.max(eeg_supp))



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
slist = sparse_list(B, 1e-1) 
adj = construct_graph(slist, AND=True)
for i in range(adj.shape[0]):
    adj[i,i] = 1

p = adj.shape[0]
k = A1.shape[0]
supp = np.zeros((p*k2, p*k2)) #pk x pk

##consturct support:
for i in range(p):
    for j in range(p):
        if adj[i,j] == 1:
            subm = np.zeros((k2,k2))
            sub_sup = tildeB[i,:,j*k:(j+1)*k]
            idx = np.where(sub_sup != 0)
            if idx[0].size>0:
                subm[:,:] = 1
            supp[i*k2:(i+1)*k2, j*k2:(j+1)*k2] = subm
            fused_supp[i*k2:(i+1)*k2, j*k2:(j+1)*k2] = subm

plt.imshow(supp)
plt.savefig('supp_movie.png')
print(supp)

print(np.max(supp))
lam = (1-supp)*1e9
idx = np.where(lam==0)
lam[idx] = .1

fused_lam = (1-fused_supp)*1e9
idx = np.where(fused_lam==0)
fused_lam[idx] = .1


#impement cross validation to select lambda
def cv_lasso (supp, X, Kfold, lam_v):
    N = X.shape[0]

    bk_size = int(N/Kfold)
    acc_loss = 0
    lam = (1-supp)*1e9
    idx = np.where(lam==0)
    lam[idx] = lam_v 
    acc_bic = 0   
    for i in range(Kfold):
        if i != (Kfold-1):
            cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
        else:
            cv_id = np.arange(i*bk_size, N, 1)
        train_id = list(set(range(N)) - set(cv_id))
        
        Y_cv = X[cv_id,:]

        Y_train = X[train_id,:] 
        estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=lam_v, verbose=1)
        estimator.fit(Y_train)
        prec = estimator.precision_
        acc_bic += compute_BIC(Y_cv, prec)
    return acc_bic






###############################################
#                                             #
#   parameter selection for individual model  #
#                                             #
###############################################


#umcomment the following to run cv for selecting eeg_l
#lambdas = np.arange(0.1,1,0.1)
#pool = Pool(10)
#start_time1 = time.time()
#results = pool.starmap(cv_lasso, zip(repeat(eeg_supp), repeat(tr_eeg), repeat(5), lambdas))
#end_time1 = time.time()
#print("Total Execution time:", end_time1 - start_time1)
#pool.close()
#min_idx = np.argmin(results)
#eeg_l = lambdas[min_idx]

eeg_l =0.1

print('[eeg] lambda:', eeg_l)

eeg_lam = (1-eeg_supp)*1e9
idx = np.where(eeg_lam==0)
eeg_lam[idx] = eeg_l

estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=eeg_lam, verbose=1)
estimator.fit(tr_eeg)
prec = estimator.precision_
np.save('eeg_individual.npy',prec)
print(slogdet(prec))
print(np.max(prec))


insample = compute_loglikelihood(tr_eeg, prec)
print('[eeg] in-sample loglikelihood:', insample)
outsample = compute_loglikelihood(te_eeg, prec)
print('[eeg] out-of-sample loglikelihood', outsample)


n_fold = 5
acc_lln = []
occ_lln = []
for i in range(n_fold):
    select_t = np.random.choice(23,23, replace=True)
    
    estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=eeg_lam, verbose=1)
    estimator.fit(tr_eeg[select_t])
    prec = estimator.precision_
    #print(prec)
    #print(slogdet(prec))
    #np.save('result_run1_eeg_movie_first23_1e_3_3_fold{}.npy'.format(i),prec)

    #print(np.max(prec))

    insample = compute_loglikelihood(tr_eeg[select_t], prec)
    #print('in-sample loglikelihood:', insample)
    outsample = compute_loglikelihood(te_eeg, prec)
    #print('out-of-sample loglikelihood', outsample)
    acc_lln.append(insample)
    occ_lln.append(outsample)
print('[eeg] average in-sample loglikelihood', np.mean(acc_lln), np.std(acc_lln))
print('[eeg] average out-of-sample loglikelihood', np.mean(occ_lln), np.std(occ_lln))

###############################################
#                                             #
#   parameter selection for latent model      #
#                                             #
###############################################

lambdas = np.arange(0.1,1,0.05)
#umcomment the following to run cv for selecting latent_l

#pool = Pool(10)
#start_time1 = time.time()
#results = pool.starmap(cv_lasso, zip(repeat(supp), repeat(tr_eeg), repeat(5), lambdas))

#end_time1 = time.time()
#print("Total Execution time:", end_time1 - start_time1)
##pool.close()
#min_idx = np.argmin(results)
latent_l = 0.1

print('[latent] lambda:', latent_l)

lam = (1-supp)*1e9
idx = np.where(lam==0)
lam[idx] = latent_l

estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=lam, verbose=1)
estimator.fit(tr_eeg)
prec = estimator.precision_
np.save('eeg_latent.npy',prec)
print(slogdet(prec))
print(np.max(prec))

insample = compute_loglikelihood(tr_eeg, prec)
print('[latent] in-sample loglikelihood:', insample)
outsample = compute_loglikelihood(te_eeg, prec)
print('[latent] out-of-sample loglikelihood', outsample)



#fit model 

acc_lln = []
occ_lln = []
for i in range(n_fold):
    select_t = np.random.choice(23,23, replace=True)
    
    estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=lam, verbose=1)
    estimator.fit(tr_eeg[select_t])
    prec = estimator.precision_

    insample = compute_loglikelihood(tr_eeg[select_t], prec)
    
    outsample = compute_loglikelihood(te_eeg, prec)
    occ_lln.append(outsample)
    acc_lln.append(insample)
print('[latent] average in-sample loglikelihood', np.mean(acc_lln), np.std(acc_lln))
print('[latent] average out-of-sample loglikelihood', np.mean(occ_lln), np.std(occ_lln))



###############################################
#                                             #
#   parameter selection for fused model       #
#                                             #
###############################################


#lambdas = np.arange(0.1,1,0.05)
#umcomment the following to run cv for selecting fused_l
#pool = Pool(10)
#start_time1 = time.time()
#results =pool.starmap(cv_lasso, zip(repeat(fused_supp), repeat(tr_eeg), repeat(5), lambdas))
#end_time1 = time.time()
#print("Total Execution time:", end_time1 - start_time1)
#pool.close()
#min_idx = np.argmin(results)
fused_l = 0.1

print('[fused] lambda:', fused_l )

#fused_lam = (1-fused_supp)*1e9
#idx = np.where(fused_lam==0)
#print((idx[0].size))
#fused_lam[idx] = .1


estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=fused_lam, verbose=1)
estimator.fit(tr_eeg)
prec = estimator.precision_
#print(prec)
np.save('eeg_fused.npy',prec)
print(slogdet(prec))
print(np.max(prec))

insample = compute_loglikelihood(tr_eeg, prec)
print('[fused] in-sample loglikelihood:', insample)
outsample = compute_loglikelihood(te_eeg, prec)
print('[fused] out-of-sample loglikelihood', outsample)


n_fold=5

acc_lln = []
occ_lln = []
for i in range(n_fold):
    select_t = np.random.choice(23,23, replace=True)
    
    estimator = QuicGraphicalLasso(score_metric="log_likelihood", lam=fused_lam, verbose=1)
    estimator.fit(tr_eeg[select_t])
    prec = estimator.precision_
    #print(prec)
    #print(slogdet(prec))
    #np.save('result_run1_eeg_movie_first23_1e_3_3_fold{}_all.npy'.format(i),prec)

    #print(np.max(prec))

    insample = compute_loglikelihood(tr_eeg[select_t], prec)
    #print('in-sample loglikelihood:', insample)
    outsample = compute_loglikelihood(te_eeg, prec)
    #print('out-of-sample loglikelihood', outsample)
    acc_lln.append(insample)
    occ_lln.append(outsample)
print('[fused] average in-sample loglikelihood', np.mean(acc_lln), np.std(acc_lln))
print('[fused] average out-of-sample loglikelihood', np.mean(occ_lln), np.std(occ_lln))





