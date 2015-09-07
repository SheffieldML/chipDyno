import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt

plt.close('all')

ax=plt.subplot2grid((3,3), (1,0), rowspan=2)

ker1 = GPy.kern.Cosine(input_dim=1, variance = .7, lengthscale=0.5)
ker2 = GPy.kern.Cosine(input_dim=1, variance = 0.5, lengthscale=1.)
ker3 = GPy.kern.Cosine(input_dim=1, variance = 0.2, lengthscale=2.)
ker1.plot(ax=ax)
ker2.plot(ax=ax)
ker3.plot(ax=ax)
plt.title('(a) Kernels')


plt.subplot2grid((3,3), (0,1), colspan=2, rowspan=3)
#k1 = GPy.kern.RBF(1,.1,1)
#k1 = GPy.kern.RBF(1,variance=.1, lengthscale=1.)
# Simulate sample paths
X = np.linspace(-10,10,501)[:,None]

#ax = plt.subplot(111)
Y = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
plt.plot(X,Y[0,:], label= '$l=0.5; a=0.7$')
Y = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
plt.plot(X,Y[0,:], label= '$l=1.0; a=0.5$')
Y = np.random.multivariate_normal(np.zeros(501),ker3.K(X),1)
plt.plot(X,Y[0,:], label= '$l=2.0; a=0.2$')
#plt.ylim([-1.2,1.2])
plt.title('(b) Sample functions')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=1, fancybox=True, shadow=True)

plt.legend(loc='upper center', bbox_to_anchor=(-0.35, 1.05), ncol=1, fancybox=True, shadow=True)
plt.show()