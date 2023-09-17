import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import scipy.stats as stats

import time

l2error= []
meanscore= []
meancrps= []
comptime= []

for kk in range(1, 101): 
 np.random.seed(kk)
 ''' function definitions '''
 def high(x):
  return (x-np.sqrt(2))*low(x)**2
 # return low(x)**2

 def low(x):
  return np.sin(8.0*np.pi*x)

 ''' Define training and test points '''
 dim = 1
 s = 2
 plot = 1
 N1 = 13
 N2 = np.array([8])
 ensemble = 1

 Nts = 1000
 Xtest = np.linspace(0,1, Nts)[:,None]
 Exact= high(Xtest)
 Low = low(Xtest)

 for ii in range(0,N2.shape[0]):
  for jj in range(0,ensemble):

   X1 = np.linspace(0,1, N1)[:,None]
   perm = np.random.permutation(N1)
   X2 = X1[perm[0:N2[ii]]]

   Y1 = low(X1)
   Y2 = high(X2)

   ''' Train level 1 '''
   start = time.time()
   k1 = GPy.kern.RBF(1)
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

   k2 = GPy.kern.RBF(1, active_dims = [1])*GPy.kern.RBF(1, active_dims = [0]) \
   + GPy.kern.RBF(1, active_dims = [0])

   m2 = GPy.models.GPRegression(X=XX, Y=Y2, kernel=k2)

   m2[".*Gaussian_noise"] = m2.Y.var()*0.01
   m2[".*Gaussian_noise"].fix()

   m2.optimize(max_iters = 500)

   m2[".*Gaussian_noise"].unfix()
   m2[".*Gaussian_noise"].constrain_positive()

   m2.optimize_restarts(20, optimizer = "bfgs",  max_iters = 1000)


   ''' Predict at test points '''
   # sample f_1 at xtest
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
   # print( "N1 = %d, N2 = %d, sample = %d, error = %e" % (N1, N2[ii], jj+1, error))

   l2error.append(error)
   meanscore.append(score)
   meancrps.append(crps)
   comptime.append(ctime)

l2error
np.mean(l2error) 
np.sort(l2error) 

meanscore 
np.mean(meanscore) 
np.sort(meanscore) 

meancrps
np.mean(meancrps) 
np.sort(meancrps) 

comptime

