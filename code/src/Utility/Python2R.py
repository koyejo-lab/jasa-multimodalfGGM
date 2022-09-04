import rpy2.robjects as robjects
import numpy as np
import scipy.io

def load_Rparams(path,filename):
    name = robjects.r['load'](path+"/"+filename)
    Robj = robjects.r[name[0]]
    Rdict = dict(zip( Robj.names, map(list, list(Robj))))

    Am_list = list()
    for Am in Rdict['Am']:
        Am_list.append(np.array(Am))
    Rdict['Am'] = Am_list
    tildeB_list = list()
    for i, Bi in enumerate(Rdict['B']):
        Bi_array = np.array(Bi)
        k = Bi_array.shape[0]
        p = int(Bi_array.shape[1] / k) + 1
        temp_Bi = np.zeros((k, p*k))
        temp_Bi[:, 0:(i*k)] = -Bi_array[:,0:(i*k)]
        temp_Bi[:, (i*k):((i+1)*k)] = np.eye(k)
        temp_Bi[:,(i+1)*k::] = -Bi_array[:,(i*k)::]
        tildeB_list.append(temp_Bi)
    Rdict['tildeB'] = np.array(tildeB_list)
    del Rdict['B']

    return Rdict
    

def load_Rdata(path,filename):
    name = robjects.r['load'](path+"/"+filename)
    Robj = robjects.r[name[0]]
    Rdict = dict(zip( Robj.names, map(list, list(Robj))))

    Am_list = list()
    for Am in Rdict['A']:
        Am_list.append(np.array(Am))
    Rdict['A'] = Am_list

    score_list = list()

    for Xm in Rdict['data']:
        score_list.append(np.array(Xm)) 
    Rdict['data'] = score_list

    return Rdict   

def save_Mat(dict_data, path, filename):
    # data must be in dictionary format
    scipy.io.savemat(path+"/"+filename+".mat", dict_data)
    print("Finish saving mat data")
    
