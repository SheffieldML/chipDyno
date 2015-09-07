import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close()

pb.ion()
#%pylab inline

#OK Done
ker_Lin = GPy.kern.Linear(input_dim=1)
ker_Br = GPy.kern.Brownian(input_dim=1)
ker_RBF = GPy.kern.RBF(input_dim=1, variance = 2.2, lengthscale=2.5)
ker_Cos = GPy.kern.Cosine(input_dim=1, variance = 2.2, lengthscale=2.5)
ker_Exp = GPy.kern.Exponential(input_dim=1, variance = 2.2, lengthscale=2.5)
ker_PerExp = GPy.kern.PeriodicExponential(input_dim=1, variance = 2.2, lengthscale=2.5)

#ker_Lin = GPy.kern.Bias(input_dim=1)
#ker_Br = GPy.kern.White(input_dim=1)
#ker_RBF = GPy.kern.RBF(input_dim=1, variance = 2.2, lengthscale=2.5)
#ker_Cos = GPy.kern.Cosine(input_dim=1, variance = 2.2, lengthscale=2.5)
#ker_Mat32 = GPy.kern.PeriodicExponential(input_dim=1, variance = 4.2, lengthscale=4.5)
#ker_Mat52 = GPy.kern.Exponential(input_dim=1, variance = 4.2, lengthscale=4.5)


X = np.linspace(-10,10,501)[:,None]
#Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
#Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
#Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)


fig = pb.figure(figsize=(14.,6.5))

plt.subplot(2, 3, 1)
K = ker_Lin.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(a). Linear Kernel')

plt.subplot(2, 3, 2)
K = ker_Br.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(b). Brownian Kernel')

plt.subplot(2, 3, 3)
K = ker_RBF.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(c). RBF Kernel')

plt.subplot(2, 3, 4)
K = ker_Cos.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(d). Cosine Kernel')

plt.subplot(2, 3, 5)
K = ker_Exp.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(e). Ornstein-Uhlenbeck Kernel')

plt.subplot(2, 3, 6)
#kern = ker1 + ker2
K = ker_PerExp.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.xlabel('(f). Periodic Exponential Kernel')
pb.draw()

