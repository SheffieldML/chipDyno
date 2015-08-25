#%config
#%matplotlib inline
#%config InlineBackend.figure_format = 'png'  # 'svg'

#import time
#start_time = time.time()

import sys
sys.path.append('/home/muhammad/mlprojects/GPclust')
import GPclust

from matplotlib import pyplot as plt
import numpy as np
#import colvb
import GPy
import pylab as pb
import pickle 
import os

#expression=np.load("expression.npy")
#file = '/home/muhammad/Dropbox/PGWork/SuraData/GSE46298_RAW/Resnorm_01_64.csv'
#data = np.genfromtxt(file, delimiter=',')
#expressionAll = data[1:,1:]
#expression = expressionAll[:45037,]
#expression = np.loadtxt(file, delimiter=',', usecols=range(1, 57))
#gene_names=np.load("gene_name_probeID.npy")

filename = "/home/muhammad/Dropbox/clusterALS/noCoregion/mModelnoCore_10000genes_M.pickle"
#filename = "model/m1Model_10000genes.pickle"
#filename = "model/m1Model_500genes.pickle"
#filename = "model/m1Model_15000genes.txt"
file = open(filename,'rb')
object_file = pickle.load(file)

m1 = object_file

#%matplotlib inline
#j=[30,30,30,30,30,30,30,30,60,60,60,60,60,60,60,60,90,90,90,90,90,90,90,90,120,120,120,120,120,120,120,120]
#j=[30,30,30,30,60,60,60,60,90,90,90,120,120,120,120]
#j=[30,60,90,120]

plt.figure(figsize=(12,18));
m1.plot(on_subplots=True, colour=True, in_a_row=False, newfig=False, min_in_cluster=1, joined=False, ylim=(-4.,4.))
pb.savefig('noCore10000genesNew.pdf')
#pb.savefig('model/10000genesPDF.pdf')

#print "Total Execution time in Second :", time.time() - start_time
#print "Total Execution time in Hour :",  (time.time() - start_time)/3600.
