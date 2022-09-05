import sys
import csv
import pickle
sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")


from model3 import flng_model
from VariableSelection import compute_BIC, compute_RSS, compute_BIC_ng, compute_RSS_ng
from Python2R import load_Rparams, load_Rdata, save_Mat
import numpy as np

import rpy2.robjects as robjects
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import time

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--s', type=int, default=10,
                    help='specigy sparsity, must be smaller than p')
parser.add_argument('--alpha', type=float, default=.1,
                    help='specify alpha')  
parser.add_argument('--init', type=str, default=.5,
                    help='est or noinit')                    
parser.add_argument('--type', type=str, default='tridiag1')
parser.add_argument('--N', type=int, default=100)
parser.add_argument('--p', type=int, default=50)
parser.add_argument('--savepath', type=str, default='../results')
parser.add_argument('--noise', type=str, default='')
parser.add_argument('--filepath', type=str, default='../data3')
parser.add_argument('--thre', type=float, default=1e-2)
args = parser.parse_args()


path=args.filepath
filename= "graph_{}_p{}_N{}{}.RData".format(args.type, args.p, args.N, args.noise)
true_param = load_Rparams(path, filename)
filename= "data_{}_p{}_N{}{}.RData".format(args.type, args.p, args.N, args.noise)
data = load_Rdata(path,filename)

p = args.p
k = 9
km = [9,9] #M=2
eta_a = 1e-4
eta_b = 1e-2
eta_ini_b = 1e-3
s = 30
alpha = min(int(true_param['alpha'][0])+2,k)
tau = 5
tol = 1e-3

X = list()
#folding data to the right shape
print("N=", data['data'][0].shape[0])
k1 = data['data'][0].shape[2]
X.append(data['data'][0].reshape((-1,p*k1)).T)
k2 = data['data'][1].shape[2]
X.append(data['data'][1].reshape((-1,p*k2)).T)

s_list = np.arange(3, int(args.p/2),3)
alpha_list = np.arange(1,k,2)
ub_loflist = list()
lb_loflist = list()
#s_list = [10]

id = 0

for i in range(2):
    ub_list = list()
    lb_list = list()
    for Am in data['A']:
        _,s,_ = np.linalg.svd(Am)
        ubound = s[0]*(2*(i+1))
        ub_list.append(ubound)
        idx = np.where(s >= 1e-3)
        lbound = s[idx[0][-1]]/(2*(i+1))
        lb_list.append(lbound)
    ub_loflist.append(ub_list)
    lb_loflist.append(lb_list)

N = args.N
Kfold = 5

d1 = len(s_list)
d2 = len(alpha_list)
d3 = len(lb_loflist)
d4 = len(ub_loflist)
bk_size = int(N / Kfold)
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
        model = flng_model(km, k, p, eta_a,eta_b, eta_ini_b, param['s'], param['alpha'], param['lb'], param['ub'])
        #split traning data
        if i != (Kfold-1):
            cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
        else:
            cv_id = np.arange(i*bk_size, N, 1)
        train_id = list(set(range(N)) - set(cv_id))
        
        Y_cv = [y[:, cv_id] for y in X]

        Y_train = [y[:, train_id] for y in X]

        
        model.fit( Y_train,tol=tol, true_A=true_param['Am'], true_tildeB = true_param['tildeB'], evaluate=True, iniA= data['A'], iniB=None, thre=args.thre)
        est_A = model.A
        est_tildeB = model.tildeB
        bic_temp += compute_BIC(est_A, est_tildeB, Y_cv)
        rss_temp += compute_RSS(est_A, est_tildeB, Y_cv)
        bic_ng += compute_BIC_ng(est_A, est_tildeB, Y_cv)
        rss_ng += compute_RSS_ng(est_A, est_tildeB, Y_cv)
    
    return list([bic_temp/Kfold, rss_temp/Kfold, bic_ng/Kfold, rss_ng/Kfold])

#specify the number of cpu
pool = Pool(cpu_count()-2)
#pool = Pool(1)
start_time1 = time.time()
results = pool.map(run_single_loop, param_batch)
end_time1 = time.time()
print("Total Execution time:", end_time1 - start_time1)
pool.close()
"""
for a,s in enumerate(s_list):
    for b,alpha in enumerate(alpha_list):
        for c,lb_list in enumerate(lb_loflist):
            for d, ub_list in enumerate(ub_loflist):
                
                for i in range(Kfold):
                    model = flng_model(km, k, p, eta_a,eta_b, eta_ini_a, eta_ini_b, s, alpha, lb_list, ub_list)
                    #split traning data
                    if i != (Kfold-1):
                        cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
                    else:
                        cv_id = np.arange(i*bk_size, N, 1)
                    train_id = np.array(list(set(range(N)) - set(cv_id)))
                    Y_cv = [y[:, cv_id] for y in Y]
                    Y_train = [y[:, train_id] for y in Y]
                    
                    model.fit( Y_train, evaluate=False, initial = False, iniA= initialA, iniB=None, thre=1e-3)
                    est_A = model.A
                    est_tildeB = model.tildeB
                    bic_m[a,b,c,d] += compute_BIC(est_A, est_tildeB, Y_cv)
"""
bic_m    = list(np.array(results)[:,0])
res_m    = list(np.array(results)[:,1])
bic_m_ng = list(np.array(results)[:,2])
rss_m_ng = list(np.array(results)[:,3])

print(bic_m)
print(res_m)
min_bic = min(bic_m)
id = bic_m.index(min(bic_m))
print("Mininum BIC", min_bic)
print(param_batch[id])

id2 = res_m.index(min(res_m))
print("Mininum BIC", min(res_m))
print(param_batch[id2])





savepath = "../results_variable/"
savefile = "model3_{}_p{}_N{}_I{}_cv_run2".format(args.type, args.p, args.N, args.init)
np.save(savepath+savefile+"_bic.npy", bic_m)
np.save(savepath+savefile+"_res.npy", res_m)
np.save(savepath+savefile+"_bic_ng.npy", bic_m_ng)
np.save(savepath+savefile+"_res_ng.npy", rss_m_ng)


with open(savepath+savefile+'_select_param.pkl', 'wb') as f:
    pickle.dump(param_batch[id], f)