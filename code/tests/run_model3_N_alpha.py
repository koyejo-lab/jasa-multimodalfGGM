import sys
import csv

sys.path.append("../src/Estimation")
sys.path.append("../src/Utility")

from model3 import flng_model
from Python2R import load_Rparams, load_Rdata, save_Mat
import numpy as np

import rpy2.robjects as robjects
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
import time

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--s', type=int, default=0,
                    help='specigy sparsity, must be smaller than p')
parser.add_argument('--alpha', type=int, default=0,
                    help='specify alpha')  
                   
parser.add_argument('--type', type=str, default='tridiag1')
parser.add_argument('--p', type=int, default=50)
parser.add_argument('--filepath', type=str, default='../data4')
parser.add_argument('--savepath', type=str, default='../results')
parser.add_argument('--noise', type=str, default='')
parser.add_argument('--thre', type=float, default=1e-3)
parser.add_argument('--k', type=int, default=9)
parser.add_argument('--lr_b', type=float, default=1e-2)
parser.add_argument('--lr_a', type=float, default=1e-3)
parser.add_argument('--lr_initb', type=float, default=1e-4)
parser.add_argument('--id', type=int, default=1)

args = parser.parse_args()




p = args.p
k = args.k
km = [9,9] #M=2
eta_a = args.lr_a
eta_b = args.lr_b
eta_ini_b = args.lr_initb
if args.alpha == 0.:
    alpha =  min(int(true_param['alpha'][0])+2,k)
else:
    alpha = args.alpha




if args.s == 0:
    s_list = np.arange(1,args.p,1)
else:
    s_list = [args.s]

#s_list = [10]
tau = 5
tol = 1e-3


id = 0
N_list = [275, 305, 344, 393, 458, 550, 688, 917, 1376, 2752, 13761]
#N_list = [352, 391, 440, 503, 587, 705, 881, 1175, 1762, 3525, 17626]
#N_list = [324, 360, 405, 462, 540, 648, 810, 1080, 1620, 3240,16200]
#N_list = [401, 445, 501, 573, 668, 802, 1003, 1337, 2006, 4012]
tpr_fpr = np.zeros((len(N_list), 2))

def run_single_model_N(N):

    path=args.filepath
    filename= "graph_{}_p{}_N{}{}_k9.RData".format(args.type, args.p, N, args.noise)
    true_param = load_Rparams(path, filename)
    filename= "data_{}_p{}_N{}{}_k9.RData".format(args.type, args.p, N, args.noise)
    data = load_Rdata(path,filename)

    X = list()
    #folding data to the right shape
    print("N=", data['data'][0].shape[0])
    k1 = data['data'][0].shape[2]
    X.append(data['data'][0].reshape((-1,p*k1)).T)
    k2 = data['data'][1].shape[2]
    X.append(data['data'][1].reshape((-1,p*k2)).T)

    ub_list = list()
    lb_list = list()
    for Am in data['A']:
        _,s,_ = np.linalg.svd(Am)
        ubound = s[0]*(2*(2+1))
        ub_list.append(ubound)
        idx = np.where(s >= 1e-3)
        lbound = s[idx[0][-1]]/(2*(2+1))
        lb_list.append(lbound)


    model = flng_model(km, k, p, eta_a, eta_b, eta_ini_b, s_list[0], alpha, lb_list, ub_list)
    tpr, fpr = model.fit( X, tol, true_A=true_param['Am'], true_tildeB = true_param['tildeB'], evaluate=True, iniA=data['A'], iniB=None, thre=args.thre)

    result = {
              'N': N,
              'tpr': tpr,
              'fpr': fpr,
              'maxAnorm': model.maxAnorm,
              'maxBnorm': model.maxBnorm,
              'maxAnorm_nor': model.maxAnorm_nor,
              'maxBnorm_nor': model.maxBnorm_nor,
              'avgAnorm_nor': model.avgAnorm_nor,
              'avgBnorm_nor': model.avgBnorm_nor,
              'avgAnorm': model.avgAnorm,
              'avgBnorm': model.avgBnorm,
              'distance': model.distance,
              'fn_res': model.fn_res
            }
    return result
#s_list = [20]
#s_list = [i+1 for i in range(4)]
#Parallelization
#pool = Pool(min(cpu_count()-10, len(s_list)))
#test single element

pool = Pool(min(15, len(N_list)))

#pool = Pool(1)
start_time1 = time.time()
results = pool.map(run_single_model_N, N_list)
end_time1 = time.time()
print("Total Execution time:", end_time1 - start_time1)
pool.close()


for id, result in enumerate(results):
    tpr_fpr[id, 0] = result['tpr']
    tpr_fpr[id, 1] = result['fpr']


savepath = args.savepath
savefile = "model3_{}_p{}_{}_alpha{}_rate{}{}{}_ss{}_ss_run2".format(args.type, args.p, args.noise, alpha,
                                                        np.log10(1/args.lr_initb).astype(int)
                                                        , np.log10(1/args.lr_a).astype(int)
                                                        , np.log10(1/args.lr_b).astype(int)
                                                        , args.id)
                                                                            
np.save(savepath+"/"+savefile+"_tprfpr.npy", tpr_fpr )

keys = results[0].keys()
with open(savepath+"/"+savefile+"_res.csv", "w", newline='') as output_file:
#with open('people.csv', 'w', newline='') as output_file:
    dict_writer = csv.DictWriter(output_file, keys)
    dict_writer.writeheader()
    dict_writer.writerows(results)