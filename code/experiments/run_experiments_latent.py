import sys
import csv
import pickle
sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")


from model3_experiment import flng_model
from VariableSelection import compute_BIC, compute_RSS, compute_BIC_ng, compute_RSS_ng
from Python2R import load_Rparams, load_Rdata, save_Mat
import numpy as np
from evaluate import sparse_list, remove_tilde
import rpy2.robjects as robjects
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import time


#implement cross validation (parse data)
k = 17
km = [17,17]
p = 68
N = 23
eta_a = 1e-5
eta_b = 1e-4
eta_ini_b = 1e-5
tol = 1e-3
Kfold = 5
bk_size = int(N / Kfold)



#load data
filename= "fmri_eeg_score_rest_all.Rdata"
path = './resting_data'


name = robjects.r['load'](path+"/"+filename)
Robj = robjects.r[name[0]]
Rdict = dict(zip( Robj.names, map(list, list(Robj))))
fmri = np.array(robjects.r['data'][0])[0:23] #use the first 23 subjects for training
eeg = np.array(robjects.r['data'][1])[0:23]

print(fmri.shape)
print(eeg.shape)
X = list()
X.append(fmri.reshape((-1, p*km[0])).T)
X.append(eeg.reshape(( -1, p*km[1])).T)
#check if the reshapse is correct
print((X[0][:,10] == fmri[10,:,:].reshape(-1)).all())
print((X[1][:,10] == eeg[10,:,:].reshape(-1)).all())


#load A matrix 
filename='Amatrix_rest_all_first23.Rdata'
name = robjects.r['load'](path+"/"+filename)
A1 = np.array(robjects.r[name[0]][0])[0]
A2 = np.array(robjects.r[name[0]][1])[0]
A_list = list()
A_list.append(A1)
A_list.append(A2)



s_list = np.arange(3,int(p/2),4)
alpha_list = np.arange(3,k,3)
ub_loflist = list()
lb_loflist = list()
#s_list = [10]

id = 0

for i in range(4):
    ub_list = list()
    lb_list = list()
    for Am in A_list:
        _,s,_ = np.linalg.svd(Am)
        ubound = s[0]*(i+1)
        ub_list.append(ubound)
        idx = np.where(s >= 1e-3)
        lbound = s[idx[0][-1]]/(i+1)
        lb_list.append(lbound)
    ub_loflist.append(ub_list)
    lb_loflist.append(lb_list)

s_list = [18]
alpha = 7
ub = ub_loflist[-1]
lb = lb_loflist[-1]


def run_single_loop(s):

    model = flng_model(km, k, p, eta_a,eta_b, eta_ini_b, s, alpha, lb, ub, 1000)


    
    model.fit( X, tol=tol, evaluate=False, iniA= A_list, iniB=None, thre=1e-3)
    B = remove_tilde(model.tildeB)
    select=sparse_list(B,tol=1e-3)

    result = {
        'Am': model.A,
        'tildeB': model.tildeB,
        'graph': select
    }
    return result


pool = Pool(1)
#pool = Pool(1)
start_time1 = time.time()
results = pool.map(run_single_loop, s_list)
end_time1 = time.time()
print("Total Execution time:", end_time1 - start_time1)
pool.close()

np.save('experiment_result_rest_all_first23.npy', results[0])
keys = results[0].keys()
with open("./experiment_result_rest_all_first23.csv", "w", newline='') as output_file:
#with open('people.csv', 'w', newline='') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(results)