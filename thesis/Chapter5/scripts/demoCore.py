import GPy
import numpy as np
import pylab as pb
import matplotlib.pyplot as plt
plt.close()


#This functions generate data corresponding to two outputs
f_output1 = lambda x: 4. * np.cos(x/5.) - .4*x - 70. + np.random.rand(x.size)[:,None] * 2.
f_output2 = lambda x: 6. * np.cos(x/5.) + .2*x + 70. + np.random.rand(x.size)[:,None] * 8.


#{X,Y} training set for each output
X1 = np.random.rand(100)[:,None]; X1=X1*75
X2 = np.random.rand(100)[:,None]; X2=X2*70 + 30
Y1 = f_output1(X1)
Y2 = f_output2(X2)
#{X,Y} test set for each output
Xt1 = np.random.rand(100)[:,None]*100
Xt2 = np.random.rand(100)[:,None]*100
Yt1 = f_output1(Xt1)
Yt2 = f_output2(Xt2)

################
xlim = (0,100); ylim = (0,50)


K = GPy.kern.Matern32(1)

m1 = GPy.models.GPRegression(X1,Y1,kernel=K.copy())
m1.optimize()
m2 = GPy.models.GPRegression(X2,Y2,kernel=K.copy())
m2.optimize()
fig = pb.figure(figsize=(12,8))
#Output 1
ax1 = fig.add_subplot(121)
m1.plot(plot_limits=xlim,ax=ax1)
ax1.plot(Xt1[:,:1],Yt1,'r+',mew=1.5)
ax1.set_title('Example 1: Without coregionalization')
#Output 2
ax2 = fig.add_subplot(122)
m2.plot(plot_limits=xlim,ax=ax2)
ax2.plot(Xt2[:,:1],Yt2,'r+',mew=1.5)
ax2.set_title('Example 2: Without coregionalization')

#################
def plot_2outputs(m,xlim,ylim):
    fig = pb.figure(figsize=(12,8))
    #Output 1
    ax1 = fig.add_subplot(121)
    ax1.set_xlim(xlim)
    ax1.set_title('Example 1: with coregionalization')
    m.plot(plot_limits=xlim,fixed_inputs=[(1,0)],which_data_rows=slice(0,100),ax=ax1,linecol='green', fillcol='green')
    ax1.plot(Xt1[:,:1],Yt1,'r+',mew=1.5)
    #Output 2
    ax2 = fig.add_subplot(122)
    ax2.set_xlim(xlim)
    ax2.set_title('Example 2: with coregionalization')
    m.plot(plot_limits=xlim,fixed_inputs=[(1,1)],which_data_rows=slice(100,200),ax=ax2,linecol='green', fillcol='green')
    ax2.plot(Xt2[:,:1],Yt2,'r+',mew=1.5)



K1 = GPy.kern.Bias(1)
K2 = GPy.kern.Linear(1)
K3 = GPy.kern.Matern32(1)
lcm = GPy.util.multioutput.LCM(input_dim=1,num_outputs=2,kernels_list=[K1,K2,K3])

m = GPy.models.GPCoregionalizedRegression([X1,X2],[Y1,Y2],kernel=lcm)
m['.*ICM.*var'].unconstrain()
m['.*ICM0.*var'].constrain_fixed(1.)
m['.*ICM0.*W'].constrain_fixed(0)
m['.*ICM1.*var'].constrain_fixed(1.)
m['.*ICM1.*W'].constrain_fixed(0)
m.optimize()
plot_2outputs(m,xlim=(0,100),ylim=(-20,60))

plt.show()

"""
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


X = np.linspace(-10,10,11)[:,None]
#Y1 = np.random.multivariate_normal(np.zeros(501),ker1.K(X),1)
#Y2 = np.random.multivariate_normal(np.zeros(501),ker2.K(X),1)
#Yn = np.random.multivariate_normal(np.zeros(501),kern.K(X),1)


fig = pb.figure(figsize=(14.,6.5))

plt.subplot(2, 3, 1)
K = ker_Lin.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(a). Linear Kernel')

plt.subplot(2, 3, 2)
K = ker_Br.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(b). Brownian Kernel')

plt.subplot(2, 3, 3)
K = ker_RBF.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(c). RBF Kernel')

plt.subplot(2, 3, 4)
K = ker_Cos.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(d). Cosine Kernel')

plt.subplot(2, 3, 5)
K = ker_Exp.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(e). OU Kernel')

plt.subplot(2, 3, 6)
#kern = ker1 + ker2
K = ker_PerExp.K(X,X)
plt.matshow(K)
plt.colorbar()
plt.xlabel('(f). Periodic Exponential Kernel')
pb.draw()

"""