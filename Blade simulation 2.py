import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import pandas as pd
import scipy.stats as stats

import time

l2error= []
meancrps= []
comptime= []

dim=2
active_dimensions = np.arange(0,dim)
 
Xtest = np.array(pd.read_table("Rmatlab_files/generate_text/Xtest.txt", sep=","))
Exact = np.array(pd.read_table("Rmatlab_files/generate_text/ytest.txt", sep=","))

#########################################

# X1 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
# Y1 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
# Y1


X2 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
Y2 = np.array(pd.read_table("Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
Y2


 ''' Train level 1 '''
start = time.time()
k1 = GPy.kern.RBF(dim, ARD = True)
m1 = GPy.models.GPRegression(X=X1, Y=Y1, kernel=k1)

m1[".*Gaussian_noise"] = m1.Y.var()*0.01
m1[".*Gaussian_noise"].fix()

m1.optimize(max_iters = 500)

m1[".*Gaussian_noise"].unfix()
m1[".*Gaussian_noise"].constrain_positive()

m1.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)

mu1, v1 = m1.predict(X2)


 ''' Train level 2 '''
XX = np.hstack((X2, mu1))

k2 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
+ GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

m2 = GPy.models.GPRegression(X=XX, Y=Y2, kernel=k2)

m2[".*Gaussian_noise"] = m2.Y.var()*0.01
m2[".*Gaussian_noise"].fix()

m2.optimize(max_iters = 500)

m2[".*Gaussian_noise"].unfix()
m2[".*Gaussian_noise"].constrain_positive()

m2.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)


''' Predict at test points '''
# sample f_1 at xtest
nsamples = 100
ntest = Xtest.shape[0]
mu1, C1 = m1.predict(Xtest, full_cov=True)
Z = np.random.multivariate_normal(mu1.flatten(),C1,nsamples)

# push samples through f_2
tmp_m = np.zeros((nsamples,ntest))
tmp_v = np.zeros((nsamples,ntest))
for i in range(0,nsamples):
 mu, v = m2.predict(np.hstack((Xtest, Z[i,:][:,None])))
 tmp_m[i,:] = mu.flatten()
 tmp_v[i,:] = v.flatten()

# get posterior mean and variance
mean = np.mean(tmp_m, axis = 0)[:,None]
var = np.mean(tmp_v, axis = 0)[:,None]+ np.var(tmp_m, axis = 0)[:,None]
var = np.abs(var)
end = time.time()

# Exact = Exact[:,None]
error = np.sqrt(np.mean((mean-Exact)**2)) 
crps = np.mean(-np.sqrt(var)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mean)/np.sqrt(var))-(Exact-mean)/np.sqrt(var)*(2*stats.norm.cdf((Exact-mean)/np.sqrt(var))-1)))
ctime = (end - start)

l2error.append(error)
meancrps.append(crps)
comptime.append(ctime)


