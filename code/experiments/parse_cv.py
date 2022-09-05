import numpy as np
import matplotlib.pyplot as plt

x = np.arange(3,34,4)

obj = np.load('./variable_selection/experiment_cv_fmri_res_ng.npy')
print(obj.shape)
for i in range(8):
    plt.plot(obj[5*i:(i+1)*5],label= "s={}".format(x[i]))
plt.legend()
plt.title('rss')
plt.savefig('fmri_rss_ng.png')
plt.close()

print(np.argmin(obj))

x = np.arange(3,34,4)
obj = np.load('./variable_selection/experiment_cv_fmri_bic_ng.npy')
print(np.min(obj))
for i in range(8):
    plt.plot(obj[5*i:(i+1)*5], label= "s={}".format(x[i]))
plt.title('BIC')

plt.legend()
plt.savefig('fmri_bic_ng.png')
print(np.argmin(obj))