import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close()

pb.ion()
#%pylab inline

X = np.linspace(10,20,100)[:,None]
#print X
#X_star = np.linspace(10,20,51)[:,None]
#X_pred = np.linspace(8,20,82)[:,None]
#x_pred = linspace(0, 4, num_pred_data)[:, None]
#print X_pred

ker = GPy.kern.RBF(input_dim=1, variance = 0.2, lengthscale=1.)
K = ker.K(X,X)
#K_star = ker.K(X,X_pred)
#K_starstar = ker.K(X_pred,X_pred)
#full_K = np.vstack([np.hstack([K, K_star]), np.hstack([K_star.T, K_starstar])])
#plt.imshow(full_K)
#plt.colorbar()
#plt.show()

#pb.savefig("/home/muhammad/Dropbox/transferReport/tex/diagrams/Cov_Structure.eps")
#pb.savefig("Cov_Structure.eps")

#pb.figure()
#plt.imshow(K)
#plt.colorbar()
#plt.show()



ker1 = GPy.kern.Matern32(input_dim=1, variance = 0.1, lengthscale=0.5)
ker2 = GPy.kern.Cosine(input_dim=1, variance = 0.2, lengthscale=2.)

kern = ker1 + ker2

X = np.linspace(-10,10,501)[:,None]
Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)


fig = pb.figure()
plt.subplot(2, 3, 1)
plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')
#plt.colorbar()
#plt.title('(a). Exponentiated Quadratic')


plt.subplot(2, 3, 2)
plt.plot(X,Y2[0,:], label= '$l=0.5; a=0.7$')
#plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')

plt.subplot(2, 3, 3)
plt.plot(X,Yn[0,:], label= '$l=0.5; a=0.7$')
#plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')

plt.subplot(2, 3, 4)
K = ker1.K(X,X)
plt.imshow(K)
#plt.colorbar()

plt.subplot(2, 3, 5)
K = ker2.K(X,X)
plt.imshow(K)
#plt.colorbar()

plt.subplot(2, 3, 6)
#kern = ker1 + ker2
K = kern.K(X,X)
plt.imshow(K)
#plt.colorbar()


pb.draw()



kern = ker1 * ker2

#X = np.linspace(-10,10,501)[:,None]
#Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
#Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)


fig = pb.figure()
plt.subplot(2, 3, 1)
plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')
#plt.colorbar()
#plt.title('(a). Exponentiated Quadratic')


plt.subplot(2, 3, 2)
plt.plot(X,Y2[0,:], label= '$l=0.5; a=0.7$')
#plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')

plt.subplot(2, 3, 3)
plt.plot(X,Yn[0,:], label= '$l=0.5; a=0.7$')
#plt.plot(X,Y1[0,:], label= '$l=0.5; a=0.7$')

plt.subplot(2, 3, 4)
K = ker1.K(X,X)
plt.imshow(K)
#plt.colorbar()

plt.subplot(2, 3, 5)
K = ker2.K(X,X)
plt.imshow(K)
#plt.colorbar()

plt.subplot(2, 3, 6)
#kern = ker1 + ker2
K = kern.K(X,X)
plt.imshow(K)
#plt.colorbar()


pb.draw()