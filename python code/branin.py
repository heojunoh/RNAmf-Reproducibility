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
 
 def high(x):
    x1 = x[:,0]
    x2 = x[:,1]
    return (-1.275*x1**2 / np.pi**2 + 5.0*x1/np.pi + x2 - 6.0)**2 + (10.0 - 5.0/(4.0*np.pi))*np.cos(x1) + 10.0

 def medium(x):
    x1 = x[:,0]
    x2 = x[:,1]
    return 10.0*np.sqrt(high(x-2.0)) + 2.0*(x1-0.5)-3.0*(3.0*x2-1.0) - 1.0

 def low(x):
    x1 = x[:,0]
    x2 = x[:,1]
    return medium(1.2*(x+2.0)) - 3.0*x2 + 1.0

 def scale_range(x,ub,lb):
    Np = x.shape[0]
    dim = x.shape[1]
    for i in range(0,Np):
        for j in range(0,dim):
            tmp = ub[j] -lb[j]
            x[i][j] = tmp*x[i][j] + lb[j]
    return x

 def rmse(pred, truth):
    pred = pred.flatten()
    truth = truth.flatten()
    return np.sqrt(np.mean((pred-truth)**2))


 ''' Create training set '''
 N1 = 20
 N2 = 15
 N3 = 10

 plot = 1
 save = 0

 dim = 2
 lb = np.array([-5.0, 0.0])
 ub = np.array([10.0, 15.0])

 tmp = np.random.rand(1000,dim)
 Xtrain = scale_range(tmp,ub,lb)
 idx = np.random.permutation(1000)
 X1 = Xtrain[idx[0:N1], :]
 X2 = Xtrain[idx[0:N2], :]
 X3 = Xtrain[idx[0:N3], :]

 Y1 = low(X1)[:,None]
 Y2 = medium(X2)[:,None]
 Y3 = high(X3)[:,None]

 nn = 40
 lb = np.array([-5.0, 0.0])
 ub = np.array([10.0, 15.0])
 x1 = np.linspace(lb[0], ub[0], 50)
 x2 = np.linspace(lb[1], ub[1], 50)
 X, Y = np.meshgrid(x1, x2)

 tmp = np.random.rand(1000,2)
 Xtest = scale_range(tmp,ub,lb)

 Exact = high(Xtest)
 Medium = medium(Xtest)
 Low = low(Xtest)

 active_dimensions = np.arange(0,dim)

 start = time.time()

 ''' Train level 1 '''
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


 # Prepare for level 3: sample f_1 at X3
 nsamples = 100
 ntest = X3.shape[0]
 mu0, C0 = m1.predict(X3, full_cov=True)
 Z = np.random.multivariate_normal(mu0.flatten(),C0,nsamples)
 tmp_m = np.zeros((nsamples,ntest))
 tmp_v = np.zeros((nsamples,ntest))

 # push samples through f_2
 for i in range(0,nsamples):
    mu, v = m2.predict(np.hstack((X3, Z[i,:][:,None])))
    tmp_m[i,:] = mu.flatten()
    tmp_v[i,:] = v.flatten()

 # get mean and variance at X3
 mu2 = np.mean(tmp_m, axis = 0)
 v2 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
 mu2 = mu2[:,None]
 v3 = np.abs(v2[:,None])


 ''' Train level 3 '''
 XX = np.hstack((X3, mu2))

 k3 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
     + GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

 m3 = GPy.models.GPRegression(X=XX, Y=Y3, kernel=k3)

 m3[".*Gaussian_noise"] = m3.Y.var()*0.01
 m3[".*Gaussian_noise"].fix()

 m3.optimize(max_iters = 500)

 m3[".*Gaussian_noise"].unfix()
 m3[".*Gaussian_noise"].constrain_positive()

 m3.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)


 # Compute posterior mean and variance for level 3 evaluated at the test points

 # sample f_1 at Xtest
 nsamples = 100
 ntest = Xtest.shape[0]
 mu0, C0 = m1.predict(Xtest, full_cov=True)
 Z = np.random.multivariate_normal(mu0.flatten(),C0,nsamples)

 # push samples through f_2 and f_3
 tmp_m = np.zeros((nsamples**2,ntest))
 tmp_v = np.zeros((nsamples**2,ntest))
 cnt = 0
 for i in range(0,nsamples):
    mu, C = m2.predict(np.hstack((Xtest, Z[i,:][:,None])), full_cov=True)
    Q = np.random.multivariate_normal(mu.flatten(),C,nsamples)
    for j in range(0,nsamples):
        mu, v = m3.predict(np.hstack((Xtest, Q[j,:][:,None])))
        tmp_m[cnt,:] = mu.flatten()
        tmp_v[cnt,:] = v.flatten()
        cnt = cnt + 1


 # get f_2 posterior mean and variance at Xtest
 mu3 = np.mean(tmp_m, axis = 0)
 v3 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
 mu3 = mu3[:,None]
 v3 = np.abs(v3[:,None])

 end = time.time()


 Exact = Exact[:,None]

 error = np.sqrt(np.mean((mu3-Exact)**2))
 score = np.mean(-(Exact-mu3)**2/v3-np.log(v3))
 crps = np.mean(-np.sqrt(v3)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mu3)/np.sqrt(v3))-(Exact-mu3)/np.sqrt(v3)*(2*stats.norm.cdf((Exact-mu3)/np.sqrt(v3))-1)))
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



# ### RMSE ### 20, 15, 10
# c(16.51563706,  17.61260983,  19.27095676,  20.03083803,
#         21.17312365,  21.49043772,  21.97356509,  22.07240125,
#         22.54659889,  23.11322589,  24.05164811,  24.24957211,
#         24.38184806,  25.37281651,  26.35344171,  26.50723802,
#         28.36562405,  28.68204586,  28.82430544,  29.0372357 ,
#         29.16104224,  29.48510606,  29.68764215,  29.71943289,
#         29.98071657,  30.18067373,  30.3015576 ,  30.38552358,
#         30.5151775 ,  30.98234873,  31.19926802,  31.28714229,
#         31.59409527,  31.59560226,  31.91475028,  32.21009853,
#         32.64194056,  32.99039166,  34.47904725,  34.65733423,
#         34.77926843,  35.10888961,  36.07956805,  36.19937603,
#         36.45751832,  36.58075663,  36.74966396,  37.3700588 ,
#         38.05483154,  38.08946518,  38.11235063,  38.45652428,
#         39.93109198,  41.07009993,  41.24307847,  42.5570682 ,
#         42.91638998,  43.08970545,  43.33807062,  43.39204447,
#         43.67477756,  44.16319876,  44.42081343,  44.65708265,
#         45.32749562,  46.12288939,  47.11973168,  47.19152415,
#         47.49545135,  47.60377949,  47.87656701,  48.00777565,
#         48.08254957,  48.15646375,  48.53462565,  49.30408324,
#         49.85945489,  50.32123495,  51.69290225,  51.71658022,
#         51.80019439,  52.09135228,  52.57640215,  52.83387251,
#         53.65493215,  53.81735274,  53.95296813,  54.1873426 ,
#         54.53113248,  57.39591913,  57.51211135,  57.82645775,
#         59.24087794,  59.33988191,  59.65029112,  59.76397285,
#         61.93991579,  62.17425702,  62.94331531, 100.68267997)
# 
# ### mean CRPS result.branin.meancrps ### The smaller, the better
# c(7.93069595,  9.43343044,  9.54151377,  9.72126493, 10.39037927,
#        10.49107708, 10.51889256, 10.78811531, 10.83320515, 11.01893776,
#        11.04393738, 11.09730405, 11.67525392, 11.96963657, 11.97959774,
#        12.07384909, 12.69938011, 13.23650216, 13.36906042, 13.81575332,
#        14.46785737, 14.56691336, 14.62767005, 14.69378076, 14.77545607,
#        14.79852786, 14.86748539, 15.07782292, 15.20813245, 15.39326357,
#        15.51186074, 15.57391663, 15.68483006, 15.84166467, 15.93522859,
#        16.41943072, 16.47444469, 16.65862628, 16.8445593 , 17.10562516,
#        17.16838018, 17.67307722, 18.03696358, 18.22482711, 18.27345742,
#        18.62117588, 18.67404549, 18.98408725, 19.26929307, 19.52907175,
#        19.7236741 , 19.89496869, 20.74817045, 20.90615823, 20.98521786,
#        21.18327277, 21.24209368, 21.27275237, 21.36112172, 21.36218409,
#        21.92558969, 22.52076681, 22.99033822, 23.01103413, 23.06385702,
#        23.17814015, 23.90019583, 23.9649098 , 24.08005458, 24.35688269,
#        25.2657552 , 25.70082389, 25.8761768 , 25.93351156, 26.4556661 ,
#        26.58978344, 26.74993292, 28.075578  , 28.21124157, 28.30089178,
#        28.52870706, 28.63183036, 28.80295901, 28.98688465, 29.23877078,
#        30.51799923, 30.54782418, 30.99171635, 31.10885636, 31.25312795,
#        31.70788774, 32.2403172 , 32.70030603, 32.95489386, 33.11272083,
#        33.14405864, 36.87733569, 38.90765285, 39.42471151, 70.28433906)
# 
# ### computation time result.branin.comptime ### The smaller, the better
# c(139.29871487617493, 162.51807618141174, 146.56250405311584, 140.66339707374573, 158.13290071487427, 
# 147.68331122398376, 155.9636242389679, 145.16405200958252, 150.1841320991516, 147.00824093818665, 
# 160.47584414482117, 150.96788907051086, 140.0498218536377, 147.64123916625977, 170.92029285430908, 
# 147.69771003723145, 142.63090920448303, 147.36768412590027, 158.81940603256226, 155.16283583641052, 
# 146.09358882904053, 156.11737513542175, 141.87225008010864, 147.26390290260315, 149.22877407073975, 
# 143.82359290122986, 139.06521797180176, 151.26320695877075, 159.24053812026978, 159.7500557899475, 
# 146.2730691432953, 137.30180287361145, 146.83052515983582, 138.1788101196289, 149.07607412338257, 
# 152.51751899719238, 158.44029903411865, 149.99088191986084, 153.45190000534058, 146.9187879562378, 
# 135.79412007331848, 137.6956911087036, 131.6539227962494, 121.37452006340027, 130.94265818595886, 
# 133.0022599697113, 131.50036096572876, 130.9718210697174, 133.38457798957825, 112.36569213867188, 
# 132.10127425193787, 115.22399830818176, 139.21120619773865, 130.9202868938446, 131.9804151058197, 
# 141.5452208518982, 156.39409589767456, 160.71244382858276, 158.8587532043457, 157.62313866615295, 
# 151.53280806541443, 158.65682983398438, 145.93891191482544, 152.2540988922119, 128.97311997413635, 
# 140.41571307182312, 124.02498006820679, 122.42486214637756, 130.61541986465454, 124.27168083190918, 
# 129.37533283233643, 134.17281985282898, 130.25748801231384, 131.57602286338806, 129.56485891342163, 
# 113.14144897460938, 109.5228500366211, 100.79314517974854, 100.6189169883728, 97.40959692001343, 
# 101.4079077243805, 96.97514390945435, 103.08406090736389, 102.71305871009827, 106.89587378501892, 
# 106.34330081939697, 99.79533576965332, 97.57815408706665, 106.1473662853241, 94.20387697219849, 
# 102.51707100868225, 101.41523885726929, 93.11203002929688, 93.22160816192627, 97.22728395462036, 
# 89.63261795043945, 95.05528521537781, 90.62356877326965, 93.9550347328186, 95.87637901306152)
# 
# 
