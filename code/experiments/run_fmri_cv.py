import sys
import csv
import pickle
sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")


from neighborhood import ng_model
from VariableSelection import compute_BIC_single, compute_RSS_single, compute_BIC_ng_single, compute_RSS_ng_single
from Python2R import load_Rparams, load_Rdata, save_Mat
import numpy as np
from evaluate import sparse_list, remove_tilde
import rpy2.robjects as robjects
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import time


#implement cross validation (parse data)
k = 21
km = [21, 13]
p = 68
N = 26
eta_b = 1e-3

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

print(fmri.shape)
print(eeg.shape)
X = list()
X.append(fmri.reshape((-1, p*km[0])).T)
#X.append(eeg.reshape(( -1, p*km[1])).T)
#check if the reshapse is correct
print((X[0][:,10] == fmri[10,:,:].reshape(-1)).all())
#print((X[1][:,10] == eeg[10,:,:].reshape(-1)).all())





#s_list = [10]

id = 0



s_list = np.arange(3,34,4)
alpha_list = np.arange(3,21,4)


param_batch = list()

for s in s_list:
    for alpha in alpha_list:
        param = {
            's': s,
            'alpha': alpha,
        }
        param_batch.append(param)


def run_single_loop(param):
    bic_temp = 0
    rss_temp = 0
    bic_ng = 0
    rss_ng = 0
    for i in range(Kfold):
        model = ng_model(km, k, p, eta_b, param['s'], param['alpha'], 200)

        if i != (Kfold-1):
            cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
        else:
            cv_id = np.arange(i*bk_size, N, 1)
        train_id = list(set(range(N)) - set(cv_id))
        
        Y_cv = [y[:, cv_id] for y in X]

        Y_train = [y[:, train_id] for y in X]
        
        model.fit( X, tol=tol, evaluate=False, iniB=None, thre=1e-3)
        est_tildeB = model.tildeB
        bic_temp += compute_BIC_single(est_tildeB, Y_cv[0])
        rss_temp += compute_RSS_single(est_tildeB, Y_cv[0])
        bic_ng += compute_BIC_ng_single(est_tildeB, Y_cv[0])
        rss_ng += compute_RSS_ng_single(est_tildeB, Y_cv[0])

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



savepath = "./variable_selection/"
savefile = "experiment_cv_fmri"
np.save(savepath+savefile+"_bic.npy", bic_m)
np.save(savepath+savefile+"_res.npy", res_m)
np.save(savepath+savefile+"_bic_ng.npy", bic_m_ng)
np.save(savepath+savefile+"_res_ng.npy", rss_m_ng)


with open(savepath+savefile+'_select_param.pkl', 'wb') as f:
    pickle.dump(param_batch[id], f)
    pickle.dump(param_batch[id2], f)    
    pickle.dump(param_batch, f)