import numpy as np
import matplotlib.pyplot as plt 
from scipy.interpolate import make_interp_spline

def remove_duplicate(x,y):
    unique_v, idx, count = np.unique(x, return_inverse=True, return_counts=True)
    m = np.zeros((unique_v.size, 2))
    m[:,0] = unique_v
    duplicate_v = np.where(count>1)[0]
    for i in duplicate_v:
        j = np.where(idx==i)[0]
        max_val = np.max(y[j])
        m[ i, 1] = max_val
    single_v = np.where(count==1)[0]
    for i in single_v:
        j = np.where(idx==i)[0]
        m[i,1] = y[j]
    return m

def plot_multiple_roc(file_list, name_list, title='ROC curve', nknot=10, scatter=False):
    for file, name in zip(file_list, name_list):
        arr = np.load(file)
        sort_fpr = np.argsort(arr[:,1])

        arr2 = np.zeros(arr.shape)
        arr2[:,0] = arr[sort_fpr,0]
        arr2[:,1] = arr[sort_fpr,1]

        m = remove_duplicate(arr2[:,1], arr2[:,0])
        #m = np.vstack((np.array([0,0]), m))
        X_Y_Spline = make_interp_spline(m[:,0], m[:,1])
        X_ = np.linspace(0, np.max(arr2[:, 1]), nknot)
        Y_ = X_Y_Spline(X_)
        plt.plot(np.hstack((0.,X_)), np.hstack((0.,Y_)), '-o', label=name)
        plt.legend()
        if scatter:
            plt.scatter(arr2[:,1], arr2[:,0])
        #plt.scatter(m[:,0], m[:,1], )

    
    x = np.arange(0,1.1,.1)
    plt.plot(x,x,color='k')
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.title(title)
    plt.show()