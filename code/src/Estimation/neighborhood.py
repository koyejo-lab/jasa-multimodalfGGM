#Description:
#This file implements the neighborhood regression method for estimating functional graphical models
#The model is similar to model3.py except that there is no latent unobserved processes.
#Hence only the estimation of B is implemented

import numpy as np
from projection import proj_A_lu, proj_B
from cov_operator import extract_blkdiag, extract_offdiag, extract_sub
from scipy import linalg
from itertools import groupby

import sys
sys.path.append("../src/Utility")
from evaluate import tpr_fpr

import time
from pathos.multiprocessing import ProcessingPool as Pool

def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def big_pseudo_invt(B,p):
    u,s,vh = linalg.svd(B)
    idx = np.where(np.abs(s) > 1e-5)

    newu = u[:,idx]
    newvh = vh[idx,:]
    news = s[idx]
    return np.kron(np.eye(p-1), np.squeeze(newu @ np.diag(1./news) @ newvh).T)



class ng_model:
    def __init__(self, km, k, p, eta_b, s, alpha, maxiter):
        
        self.tildeB = np.zeros((p, k, k*p))
        self.maskB  = np.ones((p, k, k*p))
        self.addB   = np.zeros((p, k, k*p))
        for i in range(p):
            self.tildeB[i, :, i*k:(i+1)*k] = np.eye(k)
            self.addB[i, :, i*k:(i+1)*k] = np.eye(k)
            self.maskB[i, :, i*k:(i+1)*k] = 0.
        

        self.eta_b = eta_b
        self.s = s
        self.alpha = alpha
   
      
        self.p = p
        self.k = k
        self.maxiter = maxiter
        
        #for error evaluation

        self.res = list()
        

        
        
        

    def update_B(self, scov, eta=None):
        if eta is None:
            eta = self.eta_b
        #s_1 = time.time()
       
        #s_2 = time.time()
        new_B = self.tildeB -  eta * (self.tildeB @ scov)
        #s_5 = time.time()
        new_B = new_B * self.maskB
        #s_6 = time.time()
        
        #block sparse projection, disperse projection
        new_B = proj_B(new_B, self.s, self.alpha)
        return new_B + self.addB

         


    def set_B(self, B):
        self.tildeB = B

    def fn(self,X): 

        
        BiAmX = np.einsum('ijk,kl->ijl', self.tildeB, X[0])
        return linalg.norm(BiAmX)**2 /(2*X[0].shape[1]*self.p)
        
        #return sum(list(map(mini_operation, X)))/(self.p)



    def fit(self, X, tol,  true_tildeB = None, evaluate=False, iniB=None, thre=1e-2):
        """
        fit the model from X

        Parameters
        ----------
        X: a list of data, length M, each element is a ndarray (km*p, N); 
        true_A: a list of ground truths A, length M, each element is a ndarray (km,k);
        true_tildeB:  ndarray (p, k,k*p);   
        """

        ###check dimension compatibility
        #assert(pool is not None)
        
        print("construct sample covariance")
        #construct sample covariance 
        scov = [np.einsum('ij,kj->ik',x, x) / x.shape[1] for x in X]


        res = self.fn(X)

        
        start_time = time.time()
        itr = 0
        self.res.append(res)
        start_res = res
        pre_res = 1e-5
        
        while(np.abs(pre_res-res)/(pre_res) > tol and start_res*1.5 >=res):
            #the second condition is to ensure that the algorithm terminate before going unbounded
            itr += 1
 
            self.tildeB = self.update_B(scov)
            #t3 = time.time()
            pre_res = res

            res =  self.fn(X)


            self.res.append(res)
            #print(itr, res, maxA_norm)
            if itr%1 == 0:
                print(itr, res)
            if res >= pre_res:
                break
            if itr >=self.maxiter:
                break
        
        end_time = time.time()

        print("*** Finish Optimization ***")
        print("Execution time: {}".format(end_time - start_time))

        if true_tildeB is not None:
            tpr, fpr=tpr_fpr(true_tildeB,self.tildeB, thre)
            print("TPR:", tpr,"FPR:", fpr)
            return tpr, fpr

        



##test the function class

if __name__ == '__main__':
    #build model
    #setting parameter
    p = 50
    k = 9
    km = [9, 9, 9] #M=3
    
    eta_a = 1e-5
    eta_b = 1e-5
    eta_ini_a = 1e-2
    eta_ini_b = 1e-2
    s = 20
    alpha = 0.5
    tau =1

    model = flng_model(km, k, p, eta_a,eta_b, eta_ini_a, eta_ini_b, s, alpha, tau)
    
    #generate X

    N = 1000
    X = [np.random.randn(p*k_m,N) for k_m in km]
    scov = [np.einsum('ij,kj->ik',x,x) / N for x in X]
    tol = 1e-3
    #test initialization
    #model.initialization(X, scov, tol)
    
    A_list = list()
    for k_m in km:
        A_list.append(np.random.randn(k,k_m))
    #test fitting
    model.fit(X, tol,  true_tildeB = None, evaluate=True, initial = False, iniA=A_list)


