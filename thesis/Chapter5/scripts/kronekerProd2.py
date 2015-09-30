#import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close('all')

pb.ion()
#%pylab inline

#a = np.arange(4).reshape((2,2))
a = np.arange(4).reshape((2,2))
b = np.arange(4).reshape((2,2))

#a = np.linspace(1.0, 4.0, num=9).reshape((3,3))
b = np.linspace(1.0, 4.0, num=9).reshape((3,3))

a = np.linspace(1.0, 4.0, num=4).reshape((2,2))
#b = np.linspace(1.0, 4.0, num=4).reshape((2,2))

a = np.matrix('1 0.15; .5 .7')
b = np.matrix('1 0.7 0.35; 0.6 0.9 0.6; 0.2 0.5 0.8')
#b = np.matrix('1 0.7 0.4; 0.6 0.9 0.5; 0.2 0.5 0.7')

X = np.kron(a,b)
#Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
#Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
#Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)

col = 'BuGn'
col = 'summer'
col = 'winter'
col = 'copper'
col = 'coolwarm'

plt.matshow(a, cmap= col)
plt.colorbar()
plt.xlabel('(a) ')

plt.matshow(b, cmap= col)
plt.colorbar()
plt.xlabel('(b) ')

plt.matshow(X, cmap= col)
plt.colorbar()
plt.xlabel('(c) ')


#fig = pb.figure(figsize=(14.,6.5))
#fig = pb.figure()

#plt.subplot(1, 3, 1)
#K = ker_Lin.K(X,X)
#plt.matshow(a)


pb.draw()

