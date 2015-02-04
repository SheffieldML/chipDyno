import numpy as np
import pylab as pb
import GPy
import scipy.io
import matplotlib.pyplot as plt

#%matplotlib inline

# download the software
import urllib

urllib.urlretrieve('https://github.com/sods/ods/archive/master.zip', 'master.zip')

# unzip the software
import zipfile
zip = zipfile.ZipFile('./master.zip', 'r')
for name in zip.namelist():
    zip.extract(name, '.')

# add the module location to the python path.    
import sys
sys.path.append("./ods-master/")


import pods
data = pods.datasets.spellman_yeast_cdc15()
Y = data['Y'].fillna(0) # Replace missing values with zero following Sanguinetti et al.
t = data['t']
print data['info'], data['details']
#print Y

print Y.shape


#def plot_gene(gene_name='YBR265W'): # Test gene
def plot_gene(gene_name='YER124C'): # Gene used in the paper
    plt.plot(data['t'], data['Y'][gene_name], 'rx')
    plt.title('Gene: ' + gene_name)
    plt.xlabel('time/minutes')
plot_gene('YER124C')
#plot_gene('YBR265W')


data = pods.datasets.lee_yeast_ChIP()
# set S to find relationships where p-value is less than 1e-3
S = data['Y'].T<1e-3

#print S.columns
#S=S[:,0:5]
#TF=transcription_Factors[0:7]


intersect = list(set(S.columns) & set(Y.columns))
# this is the new value for n
print 'There are n=', len(intersect), 'overlapping genes.'
# Now reset S and Y to have the overlap of genes
Y = Y[intersect].T # now an n by T matrix
S = S[intersect].T # now an n by q matrix
#S=S[:,0:5]
print S.shape



S = S[['ACE2', 'FKH2', 'NDD1', 'SWI4', 'MCM1', 'MBP1', 'SKN7', 'YAP1', 'MSN4']]
#S = S[['ACE2', 'FKH2', 'NDD1', 'SWI4','SKN7']]
#S = S[['ACE2', 'FKH2', 'NDD1', 'SWI4']]
#S = S[['ACE2']]
#print S


#print S.as_matrix(columns=None).shape
Snew= S.as_matrix(columns=None)
S = Snew.astype(int)
#plt.imshow(S[0:113,:])
#S=S[:,0:20] # Select the number of transcription factor here
#print S
#plt.imshow(S[0:7,:])
Snew.shape



# step 1, find the SVD of S.
n, q = S.shape
T = Y.shape[1]
R, Lambda1, V = np.linalg.svd(S)
# Extract first q columns for Q
Q = R[:, :q]
# remaining columns for U
U = R[:, q:]
Lambda=np.zeros(q*q).reshape(q,q)
for i in range(q):
    Lambda[i,i]=Lambda1[i]

Q.shape
U.shape
Y.shape

# Find sigma2 by looking at variance of y_u
Y_u = np.dot(U.T, Y)
sigma2 = 1./(T*(n-q))*(Y_u*Y_u).sum()
print "sigma2 found as", sigma2

# Prepare the data for processing in GPy
Y_q = np.dot(Q.T, Y) # project data onto the principal subspace of X 

# Generate the input associated with each Y, the TF and the time point.
x0, x1 = np.asarray(np.meshgrid(t.flatten(),np.arange(q)))
X = np.hstack([x0.flatten()[:, None], x1.flatten()[:, None]])
np.save("Xt.npy",X)
y = Y_q.flatten()[:, None]

print y.shape, T,q
#print m.Y.shape

kern1 = GPy.kern.RBF(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=q, active_dims=[1])
kern1.rbf.lengthscale = 70
kern1.rbf.variance.constrain_fixed()
kern2 = GPy.kern.White(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=q, active_dims=[1])
kern1 = GPy.kern.OU(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=5, active_dims=[1])
#kern1 = GPy.kern.Matern32(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=5, active_dims=[1])


#kern2 = GPy.kern.RBF(1, active_dims=[0])*GPy.kern.Coregionalize(1,q,rank=q, active_dims=[1])
#kern2.rbf.lengthscale = 20
#kern2.rbf.variance.constrain_fixed()

m = GPy.models.GPRegression(X, y, kern1+kern2)

m.Gaussian_noise.variance = sigma2
m.Gaussian_noise.variance.constrain_fixed()

#()#kern, t, y_q)
#m.sigma2 = constrained_fixed(sigma2)
#m.optimize(messages=True)

#m.constrain_positive('')
m.optimize(messages=True)
#m.optimize(messages=True, optimizer='scg')

from IPython.display import display
display(m)
#print m


#x = np.ones((20,2))
#x[:,0] = np.linspace(0,100,20)
pb.figure()
k = m.kern.K(m.X)
import pylab as pb
pb.imshow(k)
pb.colorbar()
#pb.savefig("kern_4TF")


x = np.ones((20,2))
x[:,0] = np.linspace(0,100,20)
ys = m._raw_predict(m.X,full_cov=True)

fig = pb.figure()
import pylab as pb
pb.imshow(ys[1])
pb.colorbar()


m.plot(fixed_inputs=[(1, 1)],fillcol='g', linecol='g') # this would plot ACE2.
#print m
#pb.savefig("ACE2_OU_Wh_9TF2.png")

m.plot_f(fixed_inputs=[(1, 1)],fillcol='g', linecol='g') # this would plot ACE2.
#print m
#pb.savefig("ACE2_OU_Wh_9TF2.png")

fig = pb.figure()
ax = fig.add_subplot(111)
m.plot(fixed_inputs=[(1, 0)],ax=ax) # this would plot ACE2.
m.plot(fixed_inputs=[(1, 1)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
pb.draw()
#print m


fig = pb.figure()
ax = fig.add_subplot(111)
m.plot_f(fixed_inputs=[(1, 0)],ax=ax) # this would plot ACE2.
m.plot_f(fixed_inputs=[(1, 1)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
pb.draw()
#print m

fig = pb.figure()
ax = fig.add_subplot(111)
m.plot_f(fixed_inputs=[(1, 0)],ax=ax) # this would plot ACE2.
m.plot_f(fixed_inputs=[(1, 1)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
m.plot_f(fixed_inputs=[(1, 2)],fillcol='g', linecol='g',ax=ax) # this would plot ACE2.
m.plot_f(fixed_inputs=[(1, 3)],fillcol='y', linecol='y',ax=ax) # this would plot ACE2.
m.plot_f(fixed_inputs=[(1, 4)],fillcol='c', linecol='c',ax=ax) # this would plot ACE2.
pb.draw()
#print m

#pb.close('all')

fig = pb.figure()
ax = fig.add_subplot(331)
m.plot_f(fixed_inputs=[(1, 0)],ax=ax) # this would plot ACE2.
plt.title('FKH22')

ax = fig.add_subplot(332)
m.plot_f(fixed_inputs=[(1, 1)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
plt.title('ACE2')

ax = fig.add_subplot(333)
m.plot_f(fixed_inputs=[(1, 2)],fillcol='y', linecol='y',ax=ax) # this would plot ACE2.
plt.title('NDD1')

ax = fig.add_subplot(334)
m.plot_f(fixed_inputs=[(1, 3)],fillcol='g', linecol='g',ax=ax) # this would plot ACE2.
plt.title('SWI4')

ax = fig.add_subplot(335)
m.plot_f(fixed_inputs=[(1, 4)],fillcol='c', linecol='c',ax=ax) # this would plot ACE2.
plt.title('MCM1')

ax = fig.add_subplot(336)
m.plot_f(fixed_inputs=[(1, 5)],fillcol='b', linecol='b',ax=ax) # this would plot ACE2.
plt.title('MBP1')

ax = fig.add_subplot(337)
m.plot_f(fixed_inputs=[(1, 6)],fillcol='#FF8800', linecol='#FF8800',ax=ax) # this would plot ACE2.
plt.title('SKN7')

ax = fig.add_subplot(338)
m.plot_f(fixed_inputs=[(1, 7)],fillcol='#FF0088', linecol='#FF0088',ax=ax) # this would plot ACE2.
plt.title('YAP1')

ax = fig.add_subplot(339)
m.plot_f(fixed_inputs=[(1, 8)],fillcol='#8800FF', linecol='#8800FF',ax=ax) # this would plot ACE2.
plt.title('MSN4')

pb.draw()


fig = pb.figure()
ax = fig.add_subplot(541)
m.plot_f(fixed_inputs=[(1, 0)],ax=ax) # this would plot ACE2.
ax = fig.add_subplot(542)
m.plot_f(fixed_inputs=[(1, 1)],fillcol='y', linecol='y',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(543)
m.plot_f(fixed_inputs=[(1, 2)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(544)
m.plot_f(fixed_inputs=[(1, 3)],fillcol='g', linecol='g',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(545)
m.plot_f(fixed_inputs=[(1, 4)],fillcol='c', linecol='c',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(546)
m.plot_f(fixed_inputs=[(1, 5)],fillcol='b', linecol='b',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(547)
m.plot_f(fixed_inputs=[(1, 6)],fillcol='#FF8800', linecol='#FF8800',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(548)
m.plot_f(fixed_inputs=[(1, 7)],fillcol='#FF0088', linecol='#FF0088',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(549)
m.plot_f(fixed_inputs=[(1, 8)],fillcol='#8800FF', linecol='#8800FF',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,10)
m.plot_f(fixed_inputs=[(1, 9)],ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,11)
m.plot_f(fixed_inputs=[(1, 10)],fillcol='y', linecol='y',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,12)
m.plot_f(fixed_inputs=[(1, 11)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,13)
m.plot_f(fixed_inputs=[(1, 12)],fillcol='g', linecol='g',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,14)
m.plot_f(fixed_inputs=[(1, 13)],fillcol='c', linecol='c',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,15)
m.plot_f(fixed_inputs=[(1, 14)],fillcol='b', linecol='b',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,16)
m.plot_f(fixed_inputs=[(1, 15)],fillcol='#FF8800', linecol='#FF8800',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,17)
m.plot_f(fixed_inputs=[(1, 16)],fillcol='#FF0088', linecol='#FF0088',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,18)
m.plot_f(fixed_inputs=[(1, 17)],fillcol='g', linecol='g',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,19)
m.plot_f(fixed_inputs=[(1, 18)],fillcol='#8800FF', linecol='#8800FF',ax=ax) # this would plot ACE2.
ax = fig.add_subplot(5,4,20)
m.plot_f(fixed_inputs=[(1, 19)],fillcol='r', linecol='r',ax=ax) # this would plot ACE2.
pb.draw()