#import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close()

pb.ion()
#%pylab inline

#a = np.arange(4).reshape((2,2))
a = np.arange(2)
b = np.arange(4).reshape((2,2))


X = np.kron(a,b)
#Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
#Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
#Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)


fig = pb.figure(figsize=(14.,6.5))
#fig = pb.figure()

plt.subplot(1, 3, 1)
#K = ker_Lin.K(X,X)
plt.matshow(a)
plt.colorbar()
plt.xlabel('(a). Linear Kernel')


plt.subplot(1, 3, 2)
#K = ker_Lin.K(X,X)
plt.matshow(b)
plt.colorbar()
plt.xlabel('(b). Linear Kernel')

plt.subplot(1, 3, 3)
#K = ker_Lin.K(X,X)
plt.matshow(X)
plt.colorbar()
plt.xlabel('(c). Linear Kernel')

pb.draw()

