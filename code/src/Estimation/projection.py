#Description:
##This file implements the projection functions

import numpy as np
from scipy import linalg



def proj_A(A, bound):
    """
    perform projection to satisfy the incoherence property
    _____Input________
    A: ndarray, size=(k,km)
    bound: scaler, upper bound
    _____Output_______
    projA: ndarray, size=(k,km)

    """
    row_norm = linalg.norm(A, axis=1)+1e-5 #k 
    idx = np.where(row_norm > bound)
    scale_m = np.ones(A.shape[0])
    scale_m[idx] = bound/row_norm[idx]

    return np.einsum('ij,jk->ik', np.diag(scale_m), A)


def proj_A_lu(A, lbound, ubound):
    """
    perform projection to satisfy the incoherence property
    _____Input________
    A: ndarray, size=(k,km)
    ubound: scaler, upper bound
    lbound: scaler, lower bound
    _____Output_______
    projA: ndarray, size=(k,km)

    """

    row_norm = linalg.norm(A, axis=1)+1e-5 #k 
    idx = np.where(row_norm > ubound)
    scale_m = np.ones(A.shape[0])
    scale_m[idx] = ubound/row_norm[idx]
    idx = np.where(row_norm < lbound)
    scale_m[idx] = lbound/row_norm[idx]
    return np.einsum('ij,jk->ik', np.diag(scale_m), A)

def proj_B(B, s, alpha):
    """
    perform block sparse thresholding, and disperse operator
    _____Input_____
    B: ndarray, size=(p,k,k*p)
    s: integer, assert(s <= p)
    alpha: double, (0,1]
    _____Output_____
    projB: ndarray, size=(p,k,k*p)
    """



    p = B.shape[0]
    k = B.shape[1]
    #alphak = int(alpha*k)
    if alpha < 1.:
        alphak = int(alpha*k)
    else:
        alphak = alpha
    #print(alphak)
    #perform block sparse thresholding
    B_rs= np.copy(np.reshape(B, (p,k,k,p), 'F'))
    norm_B_rs = linalg.norm(B_rs, axis=(1,2))
    sort_B_rs = np.argsort(-1*norm_B_rs, axis=1)[:,s::] #(p,p-s)
    #big_idx = np.empty((p,sort_X_rs.shape[1]*k), dtype=np.int8)
    re_idx =  np.repeat(sort_B_rs[...,None], k*k, axis=2).transpose((0, 2, 1))
    re_idx = np.expand_dims(re_idx, axis=2).reshape((p, k, k, -1))
    np.put_along_axis(B_rs, re_idx, 0., axis=3)
    

    #perform disperse operator

    col_sort = np.argsort(-1*np.abs(B_rs), axis = 1)[:, alphak::, :, :] #col sort
    row_sort = np.argsort(-1*np.abs(B_rs), axis = 2)[:, :, alphak::, :] #row_sort
    np.put_along_axis(B_rs, col_sort, 0., axis = 1)
    np.put_along_axis(B_rs, row_sort, 0., axis = 2)

    
    #reshape the matrix back
    return np.reshape(B_rs, (p, k, k*p),'F')



if __name__ == '__main__':
    ##test projection to C
    p = 100
    k = 20
    k_m = 30



    ##test projection to B
    alpha = 0.3
    alphak = int(alpha*k)

    test_B = np.concatenate([np.arange((p-1)*k*k), np.zeros(k*k)])
    test_B = np.reshape(test_B, (k,k*p), 'F')
    test_B = np.repeat(test_B[None,...], p, axis=0)
    valid_B = np.copy(test_B)
    
    valid_B[:,:,0:(p-s-1)*k] = 0.

    #truncate row

    valid_B[:,0:(k-alphak),:] = 0.
    valid_B_re = np.reshape(valid_B, (p,k,k,p), 'F')
    valid_B_re[:,:,0:(k-alphak),:] = 0.

    valid_B = np.reshape(valid_B_re,(p,k,k*p),'F')

    return_B = proj_B(test_B, s, alpha)

    
    assert (return_B == valid_B).all() == True
    print("pass test for func projB")





