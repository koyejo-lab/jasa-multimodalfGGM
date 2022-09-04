
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



class flng_model:
    def __init__(self, km, k, p, eta_a, eta_b, eta_ini_b, s, alpha, lbound, ubound, maxiter):
        
        self.A = [np.zeros((k, km[i])) for i in range(len(km))]
        self.tildeB = np.zeros((p, k, k*p))
        self.maskB  = np.ones((p, k, k*p))
        self.addB   = np.zeros((p, k, k*p))
        for i in range(p):
            self.tildeB[i, :, i*k:(i+1)*k] = np.eye(k)
            self.addB[i, :, i*k:(i+1)*k] = np.eye(k)
            self.maskB[i, :, i*k:(i+1)*k] = 0.
        
        self.eta_a = eta_a
        self.eta_b = eta_b
        self.eta_ini_b = eta_ini_b
        self.s = s
        self.alpha = alpha
   
        self.M = len(km)       
        self.p = p
        self.k = k
        self.ub = ubound
        self.lb = lbound
        self.maxiter = maxiter
        
        #for error evaluation
        self.ini_res = list()
        self.res = list()
        self.fn_res = list()
        self.distance = list()
        self.Adist = list()
        self.maxAnorm = list()
        

        
        
        

    def update_B(self, scov, eta=None):
        if eta is None:
            eta = self.eta_b
        #s_1 = time.time()
       
        #s_2 = time.time()
        def compute_Am_scov(b_Am,sc):
            temp1 = np.einsum('ij,jk->ik', b_Am, sc)
            temp2 = np.einsum('ij,kj->ik', temp1, b_Am)
            return temp2
        Am_scov = list(map(compute_Am_scov, self.big_Am, scov))
        #s_3 = time.time()
        A_scov = sum(Am_scov)
        #s_4 = time.time()
        new_B = self.tildeB -  eta * (self.tildeB @ A_scov)
        #s_5 = time.time()
        new_B = new_B * self.maskB
        #s_6 = time.time()
        
        #block sparse projection, disperse projection
        new_B = proj_B(new_B, self.s, self.alpha)
        return new_B + self.addB

    def update_B_fixA(self, A_scov, eta=None):
        if eta is None:
            eta = self.eta_ini_b
        #s_4 = time.time()
        new_B = self.tildeB -  eta * (self.tildeB @ A_scov)
        #s_5 = time.time()
        new_B = new_B * self.maskB
        #s_6 = time.time()
        
        #block sparse projection, disperse projection
        new_B = proj_B(new_B, self.s, self.alpha)
        #s_7 = time.time()

        #print("time update", s_5-s_4, "time mask",s_6-s_5, "time project",s_7-s_6)
        return new_B + self.addB

    def update_A(self, scov):
   
        #reshape from (p,k,k*p) to (p,k,k,p)
        tildeB_res = np.reshape(self.tildeB, (self.p, self.k, self.k, self.p), 'F')
        def update_per_Am(big_Am, Am, sc, lb, ub):
            #s_1 = time.time()
            #big_Am = np.kron(np.eye(self.p), Am)
            #s_2 = time.time()
            Bibig_Am = np.einsum('ijk,kl->ijl', self.tildeB, big_Am)

            Bibig_Amsc = np.einsum('ijl,lq->ijq', Bibig_Am, sc)
            #reshpae the matrix from (p,k,km*p) to (p,k,km,p)
            #s_3 = time.time()
            Bibig_Amsc_res = np.reshape(Bibig_Amsc, (self.p, self.k, -1, self.p),'F')
            #s_4 = time.time()
            new_Am = Am - self.eta_a * np.einsum('ikjl,ikql->jq', tildeB_res, Bibig_Amsc_res)
            #s_5 = time.time()

            #print(s_2-s_1, s_3-s_2, s_4-s_3, s_5-s_4)

            return proj_A_lu(new_Am,  lb, ub)
        return list(map(update_per_Am, self.big_Am, self.A, scov, self.lb, self.ub))
            

        
    
    def initialize_B(self, scov, X, tol, thre):
        #self.tildeB = np.random.normal(0., .001, size=(self.p, self.k, self.k*self.p))
        self.tildeB = np.zeros((self.p, self.k, self.k*self.p))
        # set diagonal entries
        k = self.tildeB.shape[1]
        for i in range(self.p):
            if (self.tildeB[i, :, k*i:k*(i+1)] == np.eye(k)).all() != True:
                self.tildeB[i, :, k*i:k*(i+1)] = np.eye(k)
        #update B

        big_Am = list(map(lambda x: np.kron(np.eye(self.p), x), self.A))
        #s_2 = time.time()
        def compute_Am_scov(b_Am,sc):
            temp1 = np.einsum('ij,jk->ik', b_Am, sc)
            temp2 = np.einsum('ij,kj->ik', temp1, b_Am)
            return temp2
        Am_scov = list(map(compute_Am_scov, big_Am, scov))
        #s_3 = time.time()
        A_scov = sum(Am_scov)

        res = self.fn(X)
        pre_res = 1e-5
        iter = 0
        self.ini_res.append(res)
        start_res = res
        
        while((np.abs(pre_res-res)/(pre_res)) > tol and start_res >=res):
            iter += 1
            self.tildeB = self.update_B_fixA( A_scov)
            pre_res = res
            res =  self.fn(X)
            self.ini_res.append(res)
            if iter % 1 == 0:
                print(iter, res)

            if (iter > self.maxiter):
                break         





    def set_A(self, A):
        self.A = A
    def set_B(self, B):
        self.tildeB = B

    def fn(self,X):
        """
        compute the function

        Parameters
        ----------
        X: a list of data, length M, each element is a ndarray (km*p, N); 
        A: a list of A_m, length M, each element is a ndarray (k,km);
        tildeB: a  ndarray (p, k,k*p):
        """
 

        def mini_operation(Am,Xm):
            big_Am = np.kron(np.eye(self.p), Am)
            AmX = np.einsum('ij,jk->ik', big_Am, Xm)
            Nm = Xm.shape[1]
            BiAmX = np.einsum('ijk,kl->ijl', self.tildeB, AmX)
            return linalg.norm(BiAmX)**2 /(2*Nm)
        
        return sum(list(map(mini_operation, self.A, X)))/(self.p*self.M)

    
    def residual(self, true_A, true_tildeB):

        resB = linalg.norm(self.tildeB-true_tildeB)**2
        #print("resB", resB)
        resA = sum([linalg.norm(self.A[m]-true_A[m])**2 for m in range(self.M)])
        #print("resA", resA)
        self.Adist.append(resA)
        return resA + resB


    def fit(self, X, tol, true_A=None, true_tildeB = None, evaluate=False, iniA=None, iniB=None, thre=1e-2):
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
        
        if true_A is not None and true_tildeB is not None:
            assert(all_equal([Am.shape[0] for Am in true_A]))
        print("construct sample covariance")
        #construct sample covariance 
        scov = [np.einsum('ij,kj->ik',X[i], X[i]) / X[i].shape[1] for i in range(self.M)]
        #scov = list(map(lambda x: np.einsum('ij,kj->ik',x,x)/x.shape[1], X))
        print("start initilization")

        if iniA is not None:
            print('initial B')
            self.set_A(iniA)
            
            self.initialize_B(scov,X, tol, thre)
            #self.tildeB = np.random.normal(0., 1., size=(self.p, self.k, self.k*self.p))
        else: 
            #random initialization
            tempA = [np.random.normal(0., 1, size=A.shape) for A in self.A]
            self.tildeB = np.zeros((self.p, self.k, self.k*self.p))
            
            self.A = tempA
        
        ##sanity check B_ii is an identity matrix 
        k = self.tildeB.shape[1]
        for i in range(self.p):
            if (self.tildeB[i, :, k*i:k*(i+1)] == np.eye(k)).all() != True:
                self.tildeB[i, :, k*i:k*(i+1)] = np.eye(k)

        print("Finish Initialization")
        

        if evaluate == False:
            res = self.fn(X)
        else:
            if true_A is not None and true_tildeB is not None:
                res = self.residual(true_A, true_tildeB)
                self.distance.append(res)
            else:
                raise ValueError("true A, B not provided")
        
        start_time = time.time()
        itr = 0
        self.res.append(res)
        start_res = res
        pre_res = 1e-5
        
        while(np.abs(pre_res-res)/(pre_res) > tol and start_res*1.5 >=res):
            #the second condition is to ensure that the algorithm terminate before going unbounded
            self.big_Am = list(map(lambda x: np.kron(np.eye(self.p), x), self.A))
            itr += 1
            t1 = time.time()
            self.A = self.update_A(scov)
            t2 = time.time()
            self.tildeB = self.update_B(scov)
            t3 = time.time()
            pre_res = res
            if evaluate:
                res = self.residual(true_A, true_tildeB)
                self.distance.append(res)
            else:
                res =  self.fn(X)
            maxA_norm = max([np.linalg.norm(A) for A in self.A])
            self.maxAnorm.append(maxA_norm) 
            t4 = time.time()
            #print('time emission', t2-t1, t3-t2, t4-t3)

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
    model.fit(X, tol, true_A=None, true_tildeB = None, evaluate=True, initial = False, iniA=A_list)


