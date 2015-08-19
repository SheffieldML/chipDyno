import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close()

pb.ion()
#%pylab inline

X = np.linspace(10,15,18)[:,None]
#print X
#X_star = np.linspace(10,20,51)[:,None]
X_pred = np.linspace(8,20,82)[:,None]
#x_pred = linspace(0, 4, num_pred_data)[:, None]
#print X_pred

ker1 = GPy.kern.OU(input_dim=1, variance = .7, lengthscale=0.5)
ker2 = GPy.kern.OU(input_dim=1, variance = 0.5, lengthscale=1.8)
ker3 = GPy.kern.RBF(input_dim=1, variance = 0.2, lengthscale=1.)

ker = ker2

K = ker.K(X,X)
K_star = ker.K(X,X_pred)
K_starstar = ker.K(X_pred,X_pred)
full_K = np.vstack([np.hstack([K, K_star]), np.hstack([K_star.T, K_starstar])])
plt.imshow(full_K)
plt.colorbar()
plt.show()
#pb.savefig("/home/muhammad/Dropbox/transferReport/tex/diagrams/Cov_Structure.eps")
#pb.savefig("Cov_Structure.eps")
