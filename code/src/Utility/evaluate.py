#tpr, fpr
import numpy as np

def remove_tilde(B):
    p = B.shape[0]
    k = B.shape[1]

    m = np.zeros((p, k, k*(p-1)))
    for i in range(p):
        if i==0:
            yidx = np.arange(k,k*p,1)
        elif i==p-1:
            yidx = np.arange(0,k*(p-1),1)
        else:
            idx1 = np.arange(0,k*i,1)
            idx2 = np.arange(k*(i+1),k*p,1)
            yidx = np.hstack((idx1,idx2))

        m[i,:,:] = B[i,:,yidx].T

    return(m)
    
def sparse_list(B, tol=1e-1):
    p = B.shape[0]
    k = B.shape[1]
    B2 = np.reshape(B, (p,k,k,p-1), order='F')
    Bnorm = np.linalg.norm(B2, axis=(1,2))
    #s_list = [np.where(Bnorm[i,:]>tol)[0] for i in range(p)]
    s_list = []
    for i in range(p):
        idx = np.where(Bnorm[i,:]>tol)[0]
        #correct the index
        idx2 = np.where(idx >= i)[0]
        idx[idx2] += 1
        s_list.append(idx)
    return s_list

def sparse_idx(B, tol=1e-1):
    return np.where(B>tol)

def construct_graph(s_list, AND=True):
    p = len(s_list)
    adj_m = np.zeros((p,p), dtype=np.int8)
    if AND:
        for id, s_id in enumerate (s_list):
            for x in s_id:
                if id in s_list[x]:
                    adj_m[x, id] = 1
                    adj_m[id, x] = 1

    else:
        for id, s_id in enumerate (s_list):
            for x in s_id:
                adj_m[id, x] = 1
                adj_m[x, id] = 1
    return adj_m

def tpr_fpr(tilde_trueB, tilde_estB, tol=1e-2):
    trueB = remove_tilde(tilde_trueB)
    estB = remove_tilde(tilde_estB)
    
    trueB_list = sparse_list(trueB, tol)
    estB_list = sparse_list(estB, tol)

    est_graph = construct_graph(estB_list, AND =True)
    true_graph = construct_graph(trueB_list, AND=True)
    
    TP = 0
    FP = 0
    TN = 0
    FN = 0
    p = len(trueB_list)
    for i in range(p):
        yidx = np.delete(np.arange(p), i-1)

        est_edge = est_graph[i, yidx]
      
        true_edge = true_graph[i, yidx]

        TP += np.where((est_edge & true_edge) == 1)[0].size
        FP += np.where(est_edge == 1)[0].size - np.where((est_edge & true_edge) == 1)[0].size
        TN += yidx.size - np.where((est_edge|true_edge) ==1)[0].size
        FN += np.where(true_edge  == 1)[0].size - np.where((est_edge & true_edge) == 1)[0].size
    
    TPR = TP / (TP+FN)
    FPR = FP / (FP+TN)

    return TPR, FPR



