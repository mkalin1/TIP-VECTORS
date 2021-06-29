from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import radius_neighbors_graph
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csgraph
from sklearn.cluster import SpectralClustering
import glob
from scipy.stats import gaussian_kde
import pandas as pd
from sklearn.neighbors import KernelDensity
files=glob.glob('LEU-VAL*.txt')
for filename in files:
    #newarr=[]
    newarr2=[]
    binned = np.loadtxt(filename,dtype=int)[:, 20:180]
    x=[]
    y=[]
    z=[]
    for i in range(0,360):
        for j in range(0,160):
            angle=i/2
            distance=1+(j/20)
            intensity=binned[i][j]
            for k in range(0,intensity):
                newarr2.append([angle,distance])

np.savetxt('leuval 1.txt',newarr2)