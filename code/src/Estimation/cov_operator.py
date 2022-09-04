
import numpy as np


## comput the sum of block diagonal \sum_{i=1}^pX_i^mX_i^{m T} shape (km, km)
def extract_blkdiag(Y, p):
    """
    sum of block diagonal matrix
    --------
    Y: ndarray (km*p,km*p)
    """
    km = int(Y.shape[0] / p)
    temp = np.zeros((km, km)) 
    for i in range (p):
        ext_m = Y[i*km : (i+1)*km, i*km : (i+1)*km]
        #ensure that the matrix is symmetric with small tolerance
        assert np.allclose(ext_m, ext_m.T)
        temp += ext_m 
    return temp

# a list of length M of X_{\{i}}^mX_i^{m T} shape(p,km(p-1),km)
def extract_offdiag(Y, p):

    assert np.allclose(Y, Y.T)
    km = int(Y.shape[0] / p)
    mask = np.kron(np.eye(p), np.ones((km, km)))
    idx = np.where(mask == 0.)
    ext = np.reshape(Y[idx], (p, km, km*(p-1))).transpose(0,2,1)

    return ext

# a list of length M of X_{\{i}}^mX_{\{i}}^{m T} shape(p,km(p-1),km(p-1))
def extract_sub(Y, p):
    km = int(Y.shape[0] / p)
    temp = np.zeros((p, km*(p-1), km*(p-1)))
    idx = np.arange(km * p)
    for i in range(p):
        s_idx = np.delete(idx, np.arange(i*km, (i+1)*km, 1))
        temp[i,:,:] = Y[np.ix_(s_idx, s_idx)]
    return temp


if __name__ == '__main__':
    p = 8
    km = 5
    X = np.random.randn(p*km, 10)
    scov = X @ X.T
    
    test_scov = np.copy(scov)
    test_scov *= 1- np.kron(np.eye(p), np.ones((km, km)))

    valid_scov = np.zeros((km,km))
    for i in range(p):
        valid_scov += scov[i*km:(i+1)*km, i*km:(i+1)*km]
    assert (extract_blkdiag(test_scov, p) == 0.).all()
    assert (extract_blkdiag(scov, p) == valid_scov).all()
    print("*** Pass the test for func extrac_blkdiag ***")

    ext_scov = extract_offdiag(scov, p)
    assert (ext_scov[0] == scov[km::, 0:km]).all()
    for i in range(1,p,1):
        assert (ext_scov[i, 0:km*i, :] == scov[0:km*i, km*i : km*(i+1)]).all()
        assert (ext_scov[i, km*i::, :] == scov[km*(i+1)::, km*i : km*(i+1)]).all()

    print("*** Pass the test for extract_offdiag ***")
    
    scov_2 = np.array([
        [1, 2, 3, 4, 5, 6],
        [2, 3, 4, 5, 6, 7], 
        [3, 4, 5, 6, 7, 8],
        [4, 5, 6, 7, 8, 9],
        [5, 6, 7, 8, 9, 10],
        [6, 7, 8, 9, 10,11]])
    
    test_1 = np.array([
        [5, 6, 7, 8],
        [6, 7, 8, 9],
        [7, 8, 9, 10],
        [8, 9, 10,11]])
    
    test_2 = np.array([
        [1, 2, 5, 6],
        [2, 3, 6, 7], 
        [5, 6, 9, 10],
        [6, 7, 10, 11]])
    
    test_3 = np.array([
        [1, 2, 3, 4],
        [2, 3, 4, 5],
        [3, 4, 5, 6],
        [4, 5, 6, 7]])
    sub_scov = extract_sub(scov_2, 3)
    assert (sub_scov[0] == test_1).all()
    assert (sub_scov[1] == test_2).all()
    assert (sub_scov[2] == test_3).all()

    print("*** Pass the test for extract_sub ***")   