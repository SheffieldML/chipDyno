import numpy as np
import pylab as pb
import GPy
import scipy.io
import matplotlib.pyplot as plt

#%matplotlib inline

#pb.close('all')

data = GPy.util.datasets.spellman_yeast_cdc15()
Y = data['Y'].fillna(0) # Replace missing values with zero following Sanguinetti et al.
t = data['t']
print data['info'], data['details']
#print Y

#def plot_gene(gene_name='YBR265W'): # Test gene
def plot_gene(gene_name='YER124C'): # Gene used in the paper
    plt.plot(data['t'], data['Y'][gene_name], 'rx')
    plt.title('Gene: ' + gene_name)
    plt.xlabel('time/minutes')
plot_gene('YER124C')
#plot_gene('YBR265W')

data = GPy.util.datasets.lee_yeast_ChIP()
# set S to find relationships where p-value is less than 1e-3
S = data['Y'].T<1e-3
#S.shape

intersect = list(set(S.columns) & set(Y.columns))
# this is the new value for n
print 'There are n=', len(intersect), 'overlapping genes.'
# Now reset S and Y to have the overlap of genes
Y = Y[intersect].T # now an n by T matrix
S = S[intersect].T # now an n by q matrix
S.shape
#S.shape

# step 1, find the SVD of S.
n, q = S.shape
T = Y.shape[1]
R, Lambda, V = np.linalg.svd(S)
# Extract first q columns for Q
Q = R[:, :q]
# remaining columns for U
U = R[:, q:]


# Find sigma2 by looking at variance of y_u
Y_u = np.dot(U.T, Y)
sigma2 = 1./(T*(n-q))*(Y_u*Y_u).sum()
print "sigma2 found as", sigma2


# Prepare the data for processing in GPy
Y_q = np.dot(Q.T, Y) # project data onto the principal subspace of X 

# Generate the input associated with each Y, the TF and the time point.
x0, x1 = np.asarray(np.meshgrid(t.flatten(),np.arange(q)))
X = np.hstack([x0.flatten()[:, None], x1.flatten()[:, None]])
y = Y_q.flatten()[:, None]


kern = GPy.kern.RBF(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=4, active_dims=[1])
#kern = GPy.kern.OU(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=5, active_dims=[1])
m = GPy.models.GPRegression(X, y, kern)
m.mul.rbf.variance = 0.1
m.mul.rbf.lengthscale = 50
m.Gaussian_noise.variance = sigma2
#m.Gaussian_noise.variance.constrain_fixed()
#()#kern, t, y_q)
#m.sigma2 = constrained_fixed(sigma2)
#m.optimize(messages=True)


print m


m.optimize(messages=True)


m.plot(fixed_inputs=[(1, 1)]) # this would plot ACE2.
#pb.savefig('/home/muhammad/Desktop/YBR265W.png')
#pb.savefig('YER124C_rbf_r7_01_50.eps')
print m

