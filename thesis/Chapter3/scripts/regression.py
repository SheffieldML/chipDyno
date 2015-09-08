from numpy.random import seed
from infpy.gp import GaussianProcess, gp_1D_X_range, gp_plot_samples_from
from pylab import plot, savefig, title, close, figure, xlabel, ylabel

from matplotlib import pyplot as plt
import pylab as pb



# seed RNG to make reproducible and close all existing plot windows
seed(12)
close('all')

#
# Kernel
#
from infpy.gp import SquaredExponentialKernel as SE
from infpy.gp import Matern32Kernel as M32

kernel = SE([1])
#kernel = M32([1])

#
# Part of X-space we will plot samples from
#
support = gp_1D_X_range(-8., 8.1, .1)

#
# Plot samples from prior
#
#figure()
#plt.figure(figsize=(18,18))
fig = pb.figure(figsize=(9,17))
ax = fig.add_subplot(511)
gp = GaussianProcess([], [], kernel)
gp_plot_samples_from(gp, support, num_samples=30)
xlabel('a). t')
ylabel('f')
#title('Samples from the prior')
#savefig('samples_from_prior.png')
#savefig('samples_from_prior.eps')


#
# Data
#
X = [[-2.5]]
Y = [3.]
ax = fig.add_subplot(512)
#
# Plot samples from posterior
#
#figure()
plot([x[0] for x in X], Y, '+',mew=3)
gp = GaussianProcess(X, Y, kernel)
gp_plot_samples_from(gp, support, num_samples=15)
xlabel('b). t = -2.5')
ylabel('f')
#title('Samples from the posterior')
#savefig('samples_from_posterior.png')
#savefig('samples_from_posterior.eps')

#
# Data
#
X = [[-2.5], [4.5]]
Y = [3., -2.9]
ax = fig.add_subplot(513)
#
# Plot samples from posterior
#
#figure()
plot([x[0] for x in X], Y, '+',mew=3)
gp = GaussianProcess(X, Y, kernel)
gp_plot_samples_from(gp, support, num_samples=10)
xlabel('c). t = -2.5, 4.5')
ylabel('f')
#title('Samples from the posterior')
#savefig('samples_from_posterior.png')
#savefig('samples_from_posterior.eps')


#
# Data
#
X = [[-4.5], [-2.5], [1.3], [4.5]]
Y = [-1.9,  3., .5, -2.9]
ax = fig.add_subplot(514)
#
# Plot samples from posterior
#
#figure()
plot([x[0] for x in X], Y, '+',mew=3)
gp = GaussianProcess(X, Y, kernel)
gp_plot_samples_from(gp, support, num_samples=8)
xlabel('d). t = -4.5, -2.5, 1.3, 4.5 ')
ylabel('f')
#title('Samples from the posterior')
#savefig('samples_from_posterior.png')
#savefig('samples_from_posterior.eps')



#
# Data
#
X = [[-4.5], [-2.5], [-.5], [1.3], [2.8],[4.5]]
Y = [-1.9,  3., -3., .5, .7, -2.9]
ax = fig.add_subplot(515)
#
# Plot samples from posterior
#
#figure()
plot([x[0] for x in X], Y, '+',mew=3)
gp = GaussianProcess(X, Y, kernel)
gp_plot_samples_from(gp, support, num_samples=4)
xlabel('e). t = -4.5, -2.5, -.5, 1.3, 2.8, 4.5')
ylabel('f')
#title('Samples from the posterior')
#savefig('samples_from_posterior.png')
#savefig('samples_from_posterior.eps')


pb.show()