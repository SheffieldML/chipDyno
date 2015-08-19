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

pb.figure()
plt.imshow(K)
plt.colorbar()
plt.show()


fig = pb.figure()
plt.subplot(2, 2, 1)
ker = GPy.kern.RBF(input_dim=1, variance = 0.2, lengthscale=1.)
K = ker.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.title('(a). Exponentiated Quadratic')

plt.subplot(2, 2, 2)
ker = GPy.kern.Matern52(input_dim=1, variance = 0.2, lengthscale=1.)
K = ker.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.title('(b). Matern52 kernel')

plt.subplot(2, 2, 3)
ker = GPy.kern.OU(input_dim=1, variance = 0.2, lengthscale=1.)
K = ker.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.title('(c). Ornstein-Uhlenbeck')

plt.subplot(2, 2, 4)
ker = GPy.kern.Cosine(input_dim=1, variance = 0.2, lengthscale=1.)
K = ker.K(X,X)
plt.imshow(K)
plt.colorbar()
plt.title('(d). Cosine kernel')


pb.draw()