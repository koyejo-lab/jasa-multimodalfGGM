import sys
import csv
import pickle
sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")


from model3_experiment import flng_model
from VariableSelection import compute_BIC, compute_RSS, compute_BIC_ng, compute_RSS_ng
from Python2R import load_Rparams, load_Rdata, save_Mat
import numpy as np

import rpy2.robjects as robjects
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import time

##remember to change this
#k = 14
#km = [23,14]
#p = 68
#N = 75

##remember to change this
k = 13
km = [21, 13]
p = 68
N = 26

#implement cross validation (parse data)

eta_a = 1e-4
eta_b = 1e-3
eta_ini_b = 1e-3
tol = 1e-3
Kfold = 5
bk_size = int(N / Kfold)


#load data
filename= "fmri_eeg_score.Rdata"
path = '/home/kt14/workbench/multimodal_functional_ggm/experiments/resting_data'


name = robjects.r['load'](path+"/"+filename)
Robj = robjects.r[name[0]]
Rdict = dict(zip( Robj.names, map(list, list(Robj))))
fmri = np.array(robjects.r['data'][0])
eeg = np.array(robjects.r['data'][1])
print('shape')
print(fmri.shape)
print(eeg.shape)
X = list()
X.append(fmri.reshape((-1, p*km[0])).T)
X.append(eeg.reshape(( -1, p*km[1])).T)
#check if the reshapse is correct
print((X[0][:,10] == fmri[10,:,:].reshape(-1)).all())
print((X[1][:,10] == eeg[10,:,:].reshape(-1)).all())




#load A matrix 
filename='Amatrix.Rdata'
name = robjects.r['load'](path+"/"+filename)
A1 = np.array(robjects.r[name[0]][0])[0]
A2 = np.array(robjects.r[name[0]][1])[0]
A_list = list()
A_list.append(A1)
A_list.append(A2)



s_list = np.arange(3,34,5)
alpha_list = np.arange(3,k,4)
ub_loflist = list()
lb_loflist = list()
#s_list = [10]

id = 0

for i in range(2):
    ub_list = list()
    lb_list = list()
    for Am in A_list:
        _,s,_ = np.linalg.svd(Am)
        ubound = s[0]*(2*(i+1))
        ub_list.append(ubound)
        idx = np.where(s >= 1e-3)
        lbound = s[idx[0][-1]]/(2*(i+1))
        lb_list.append(lbound)
    ub_loflist.append(ub_list)
    lb_loflist.append(lb_list)


d1 = len(s_list)
d2 = len(alpha_list)
d3 = len(lb_loflist)
d4 = len(ub_loflist)
bic_m = np.zeros((d1,d2,d3,d4))



param_batch = list()

for s in s_list:
    for alpha in alpha_list:
        for lb_list in lb_loflist:
            for ub_list in ub_loflist:
                param = {
                    's': s,
                    'alpha': alpha,
                    'lb': lb_list,
                    'ub': ub_list
                }
                param_batch.append(param)

def run_single_loop(param):
    bic_temp = 0
    rss_temp = 0
    bic_ng = 0
    rss_ng = 0
    for i in range(Kfold):
        model = flng_model(km, k, p, eta_a,eta_b, eta_ini_b, param['s'], param['alpha'], param['lb'], param['ub'], 50)
        #split traning data
        if i != (Kfold-1):
            cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
        else:
            cv_id = np.arange(i*bk_size, N, 1)
        train_id = list(set(range(N)) - set(cv_id))
        
        Y_cv = [y[:, cv_id] for y in X]

        Y_train = [y[:, train_id] for y in X]

        
        model.fit( Y_train,tol=tol, evaluate=False, iniA= A_list, iniB=None, thre=1e-3)
        est_A = model.A
        est_tildeB = model.tildeB
        bic_temp += compute_BIC(est_A, est_tildeB, Y_cv)
        rss_temp += compute_RSS(est_A, est_tildeB, Y_cv)
        bic_ng += compute_BIC_ng(est_A, est_tildeB, Y_cv)
        rss_ng += compute_RSS_ng(est_A, est_tildeB, Y_cv)
    
    return list([bic_temp/Kfold, rss_temp/Kfold, bic_ng/Kfold, rss_ng/Kfold])


pool = Pool(cpu_count()-1)
#pool = Pool(1)
start_time1 = time.time()
results = pool.map(run_single_loop, param_batch)
end_time1 = time.time()
print("Total Execution time:", end_time1 - start_time1)
pool.close()



bic_m    = list(np.array(results)[:,0])
res_m    = list(np.array(results)[:,1])
bic_m_ng = list(np.array(results)[:,2])
rss_m_ng = list(np.array(results)[:,3])

min_bic = min(bic_m)
id = bic_m.index(min(bic_m))
print("Mininum BIC", min_bic)
print(param_batch[id])

id2 = res_m.index(min(res_m))
print("Mininum rss", min(res_m))
print(param_batch[id2])



savepath = "./"
savefile = "experiment_cv_small_332_iter50_new"
np.save(savepath+savefile+"_bic.npy", bic_m)
np.save(savepath+savefile+"_res.npy", res_m)
np.save(savepath+savefile+"_bic_ng.npy", bic_m_ng)
np.save(savepath+savefile+"_res_ng.npy", rss_m_ng)


with open(savepath+savefile+'_select_param_smaller_set_full_332_iter50_new.pkl', 'wb') as f:
    pickle.dump(param_batch[id], f)
    pickle.dump(param_batch[id2], f)    
    pickle.dump(param_batch, f)



"""
for a,s in enumerate(s_list):
    for b,alpha in enumerate(alpha_list):
        for c,lb_list in enumerate(lb_loflist):
            for d, ub_list in enumerate(ub_loflist):
                
                for i in range(Kfold):
                    model = flng_model(km, k, p, eta_a,eta_b, eta_ini_b, s, alpha, lb_list, ub_list)
                    #split traning data
                    if i != (Kfold-1):
                        cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
                    else:
                        cv_id = np.arange(i*bk_size, N, 1)
                    train_id = np.array(list(set(range(N)) - set(cv_id)))
                    Y_cv = [y[:, cv_id] for y in X]
                    Y_train = [y[:, train_id] for y in X]
                    
                    model.fit( Y_train,tol=tol, evaluate=False, iniA= A_list, iniB=None, thre=1e-3)
                    est_A = model.A
                    est_tildeB = model.tildeB
                    bic_m[a,b,c,d] += compute_BIC(est_A, est_tildeB, Y_cv)
"""
