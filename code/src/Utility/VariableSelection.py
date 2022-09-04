#select s, alpha, tau (selection of model parameter)
#select k,km (selection of basis numer)
import numpy as np
from evaluate import sparse_list, remove_tilde, sparse_idx
from multiprocessing import Pool
import time
import sys
sys.path.append("../Estimation")


from model3 import flng_model

def logdet(estA, tilde_estB, X):
    """
    compute the function

    Parameters
    ----------
    X: a list of data, length M, each element is a ndarray (km*p, N); 
    estA: a list of A_m, length M, each element is a ndarray (k,km);
    tilde_estB: a  ndarray (p, k,k*p):
    """
    p = tilde_estB.shape[0]

    def mini_operation(Am,Xm):
        big_Am = np.kron(np.eye(p), Am)
        AmX = np.einsum('ij,jk->ik', big_Am, Xm)
        Nm = Xm.shape[1]
        BiAmX = np.einsum('ijk,kl->ijl', tilde_estB, AmX)
        empCov = np.einsum('ijk,iqk->ijq', BiAmX, BiAmX) / Nm
        
        sum_logdet = 0
        for cov in empCov:
            w, _ = np.linalg.eigh(cov)
            idx = np.where(w > 1e-5)
            
            sum_logdet += np.sum(np.log10(w[idx]))
        
        #return sum([np.log10(np.linalg.det(cov) + 1e-5) for cov in empCov])*Nm
        
        return sum_logdet *Nm
    
    return sum(list(map(mini_operation, estA, X)))

def logdet_single(tilde_estB, Xm):
    """
    compute the function

    Parameters
    ----------
    X: a list of data, length M, each element is a ndarray (km*p, N); 
    tilde_estB: a  ndarray (p, k,k*p):
    """

    Nm = Xm.shape[1]
    BiX = np.einsum('ijk,kl->ijl', tilde_estB, Xm)
    empCov = np.einsum('ijk,iqk->ijq', BiX, BiX) / Nm
    
 
    w, _ = np.linalg.eigh(empCov)
    idx = np.where(w > 1e-5)
    
    sum_logdet = np.sum(np.log10(w[idx]))
        
        #return sum([np.log10(np.linalg.det(cov) + 1e-5) for cov in empCov])*Nm
        
    return sum_logdet *Nm
    
    

def res(estA, tilde_estB, X):
    """
    compute the function

    Parameters
    ----------
    X: a list of data, length M, each element is a ndarray (km*p, N); 
    estA: a list of A_m, length M, each element is a ndarray (k,km);
    tilde_estB: a  ndarray (p, k,k*p):
    """
    p = tilde_estB.shape[0]

    def mini_operation(Am,Xm):
        big_Am = np.kron(np.eye(p), Am)
        AmX = np.einsum('ij,jk->ik', big_Am, Xm)
        Nm = Xm.shape[1]
        BiAmX = np.einsum('ijk,kl->ijl', tilde_estB, AmX)
        
        
        return np.sum(BiAmX**2)/(2)
    
    return sum(list(map(mini_operation, estA, X)))

def res_single(tilde_estB, Xm): 
    
    BiAmX = np.einsum('ijk,kl->ijl', tilde_estB, Xm)
    return np.linalg.norm(BiAmX)**2 / 2


def compute_BIC(estA, tilde_estB, X):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = 0
    for Xm in X:
        N += Xm.shape[1]
    #compute the residual

    term1 = logdet(estA, tilde_estB, X)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    #s_list = sparse_list(est_B, tol=1e-3)
    #for arr in s_list:
    #    s += arr.size
    s = sparse_idx(est_B, tol=1e-3)[0].size
    a = 0
    for Am in estA:
        a += sparse_idx(Am, tol=1e-3)[0].size
    print('BIC scale',term1, (s+a)*np.log10(N))
    return term1 + (s+a)*np.log10(N)

def compute_BIC_single(tilde_estB, Xm):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = Xm.shape[1]
    #compute the residual

    term1 = logdet_single(tilde_estB, Xm)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)

    s = sparse_idx(est_B, tol=1e-3)[0].size

    print('BIC scale',term1, s*np.log10(N))
    return term1 + s*np.log10(N)

def compute_BIC_ng(estA, tilde_estB, X):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = 0
    for Xm in X:
        N += Xm.shape[1]
    #compute the residual

    term1 = logdet(estA, tilde_estB, X)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    #s_list = sparse_list(est_B, tol=1e-3)
    #for arr in s_list:
    #    s += arr.size
    s_list = sparse_list(est_B, tol=1e-3)[0]
    s = sum([x.size for x in s_list])
    print('BIC ng', term1, (s)*np.log10(N))   
    return term1 + (s)*np.log10(N)

def compute_BIC_ng_single(tilde_estB, Xm):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = Xm.shape[1]
    #compute the residual

    term1 = logdet_single(tilde_estB, Xm)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    #s_list = sparse_list(est_B, tol=1e-3)
    #for arr in s_list:
    #    s += arr.size
    s_list = sparse_list(est_B, tol=1e-3)[0]
    s = sum([x.size for x in s_list])
    print('BIC ng', term1, (s)*np.log10(N))   
    return term1 + (s)*np.log10(N)


def compute_RSS(estA, tilde_estB, X):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = 0
    for Xm in X:
        N += Xm.shape[1]
    #compute the residual

    term1 = res(estA, tilde_estB, X)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    #s_list = sparse_list(est_B, tol=1e-3)
    #for arr in s_list:
    #    s += arr.size
    s = sparse_idx(est_B, tol=1e-3)[0].size
    a = 0
    for Am in estA:
        a += sparse_idx(Am, tol=1e-3)[0].size
    print('rss', term1, (s+a)*np.log10(N))
    return term1 + (s+a)*np.log10(N)


def compute_RSS_single(tilde_estB, Xm):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = Xm.shape[1]
    #compute the residual

    term1 = res_single(tilde_estB, Xm)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    #s_list = sparse_list(est_B, tol=1e-3)
    #for arr in s_list:
    #    s += arr.size
    s = sparse_idx(est_B, tol=1e-3)[0].size
    print('rss', term1, s*np.log10(N))
    return term1 + s*np.log10(N)


def compute_RSS_ng(estA, tilde_estB, X):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = 0
    for Xm in X:
        N += Xm.shape[1]
    #compute the residual

    term1 = res(estA, tilde_estB, X)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    s_list = sparse_list(est_B, tol=1e-3)[0]
    s = sum([x.size for x in s_list])
    print('rss ng', term1, (s)*np.log10(N))
    return term1 + (s)*np.log10(N)

def compute_RSS_ng_single(tilde_estB, Xm):
    """
    estA: a list of ndarray, length M
    tilde_estB: a ndarray (p,k,p*k)
    X: a list of data, length M, (km*p, N)
    """
    N = Xm.shape[1]
    #compute the residual

    term1 = res_single(tilde_estB, Xm)
    

    #compute the sparsity pattern
    s = 0
    est_B = remove_tilde(tilde_estB)
    s_list = sparse_list(est_B, tol=1e-3)[0]
    s = sum([x.size for x in s_list])
    print('rss ng', term1, (s)*np.log10(N))
    return term1 + (s)*np.log10(N)

def CVBIC(Y, initialA, Kfold, N, km, k, p, eta_a, eta_b, eta_ini_a, eta_ini_b, s_list, alpha_list, lb_loflist, ub_loflist):
    ####

    #K: K-fold SCV
    ####
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
        for i in range(Kfold):
            model = flng_model(km, k, p, eta_a,eta_b, eta_ini_a, eta_ini_b, param['s'], param['alpha'], param['lb'], param['ub'])
            #split traning data
            if i != (Kfold-1):
                cv_id = np.arange(i*bk_size, (i+1)*bk_size, 1)
            else:
                cv_id = np.arange(i*bk_size, N, 1)
            train_id = np.array(list(set(range(N)) - set(cv_id)))
            Y_cv = [y[:, cv_id] for y in Y]
            Y_train = [y[:, train_id] for y in Y]
            
            model.fit( Y_train, tol=1e-3, evaluate=False, initial = False, iniA= initialA, iniB=None, thre=1e-3)
            est_A = model.A
            est_tildeB = model.tildeB
            bic_temp += compute_BIC(est_A, est_tildeB, Y_cv)

        return bic_temp

    pool = Pool(8)
    start_time1 = time.time()
    bic_m = pool.map(run_single_loop, param_batch)
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
    bic_m = list(bic_m)
    min_bic = min(bic_m)
    id = bic_m.index(min(bic_m))
    print("Mininum BIC", bic_m)
    print(param_batch[id])
    return min_bic, param_batch[id]


