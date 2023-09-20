import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import scipy.stats as stats
import rpy2.robjects as robjects

import time

''' function definitions '''

def high(x):
 x1 = x[:,0]
 x2 = x[:,1]
 return (1 - np.exp(-1/(2*x2))) * (2300*x1**3 + 1900*x1**2 + 2092*x1 + 60) / (100*x1**3 + 500*x1**2 + 4*x1 +20)

def low(x):
 x1 = x[:,0]
 x2 = x[:,1]
 xx11 = x1+0.05
 xx12 = x1-0.05
 xx21 = x2+0.05
 xx22 = np.maximum([0], [x2-0.05])
 return 1/4*(high(np.vstack((xx11,xx21)).T)+high(np.vstack((xx11,xx22)).T)+high(np.vstack((xx12,xx21)).T)+high(np.vstack((xx12,xx22)).T))

''' load data '''
r = robjects.r
loaded_data = r.readRDS("/Users/junoh/Downloads/tmp_data_currin.rds")
X1=np.array(loaded_data[0])
X2=np.array(loaded_data[1])
Y1=np.array(loaded_data[2])[:,None]
Y2=np.array(loaded_data[3])[:,None]
Xtest=np.array(loaded_data[4])
Exact=np.array(loaded_data[5])[:,None]

''' Define training and test points '''
dim = 2
plot = 1
N1 = 20
N2 = 10
ensemble = 1

Nts = 1000
active_dimensions = np.arange(0,dim)

''' Train level 1 '''
start = time.time()
k1 = GPy.kern.RBF(dim, ARD = True)
m1 = GPy.models.GPRegression(X=X1, Y=Y1, kernel=k1)

m1[".*Gaussian_noise"] = m1.Y.var()*0.01
m1[".*Gaussian_noise"].fix()

m1.optimize(max_iters = 500)

m1[".*Gaussian_noise"].unfix()
m1[".*Gaussian_noise"].constrain_positive()

m1.optimize_restarts(20, optimizer = "bfgs",  max_iters = 1000)

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

m2.optimize_restarts(20, optimizer = "bfgs",  max_iters = 1000)


''' Predict at test points '''
#sample f_1 at xtest
nsamples = 1000
mu1, C1 = m1.predict(Xtest, full_cov=True)
Z = np.random.multivariate_normal(mu1.flatten(),C1,nsamples)

# push samples through f_2
tmp_m = np.zeros((nsamples,Nts))
tmp_v = np.zeros((nsamples,Nts))
for i in range(0,nsamples):
 mu, v = m2.predict(np.hstack((Xtest, Z[i,:][:,None])))
 tmp_m[i,:] = mu.flatten()
 tmp_v[i,:] = v.flatten()

# get posterior mean and variance
mean = np.mean(tmp_m, axis = 0)[:,None]
var = np.mean(tmp_v, axis = 0)[:,None]+ np.var(tmp_m, axis = 0)[:,None]
var = np.abs(var)
end = time.time()

error = np.sqrt(np.mean((mean-Exact)**2))
score = np.mean(-(Exact-mean)**2/var-np.log(var))
crps = np.mean(-np.sqrt(var)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mean)/np.sqrt(var))-(Exact-mean)/np.sqrt(var)*(2*stats.norm.cdf((Exact-mean)/np.sqrt(var))-1)))
ctime = (end - start)
