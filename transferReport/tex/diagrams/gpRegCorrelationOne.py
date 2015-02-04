import numpy as np
import pickle
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import decimal


filename = "/home/muhammad/Dropbox/CElegans/ranking/pyGpReg/OutRanScores_Sample1"
file = open(filename,'rb')
rs1 = pickle.load(file)

filename = "/home/muhammad/Dropbox/CElegans/ranking/pyGpReg/OutRanScores_Sample2"
file = open(filename,'rb')
rs2 = pickle.load(file)

filename = "/home/muhammad/Dropbox/CElegans/ranking/pyGpReg/OutRanScores_Sample3"
file = open(filename,'rb')
rs3 = pickle.load(file)

#plt.plot(rs1)
resCorr12 = pearsonr(rs1, rs2)
print 'Pearson correlation between sample 1 and 2:', resCorr12[0]

resCorr13 = pearsonr(rs1, rs3)
print 'Pearson correlation between sample 1 and 3:', resCorr13[0]

resCorr23 = pearsonr(rs2, rs3)
print 'Pearson correlation between sample 2 and 3:', resCorr23[0]

plt.close('all')


N = rs2.shape[0]
colors = np.random.rand(N)


#plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=1)
#plt.subplot(311)
ax = plt.gca()
ax.scatter(rs2,rs1, c=colors)
#plt.title('RS1 Vs RS2')
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlim((-10,1000))
plt.ylim((-10,1000))


resCorr12 =round(resCorr12[0],3)

plt.xlabel('Ranking Score 2')
plt.ylabel('Ranking Score 1')
plt.title('Scatter plot of ranking score of different genes \n from replication 1 and replication2')
plt.text(150,400,'$\\rho_{12}$ =%1.30s'%resCorr12,fontsize=18, color='red')

#plt.text(1200,400,'$\\rho$ =%1.30s'%resCorr12,fontsize=18)
plt.show()

