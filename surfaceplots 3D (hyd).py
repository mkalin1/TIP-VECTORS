from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os.path 
from os import path
import math
import glob
import pandas as pd
from matplotlib.ticker import FuncFormatter
from sklearn.preprocessing import normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LogNorm
files=glob.glob('*-*.txt')
for filename in files[:5]:
    full = np.loadtxt(filename)[:, 20:180]
    resname=filename.split(" ")[0]
    #anglemean=filename.split(' ')[3]
    #anglevar=filename.split(' ')[4]
    #distmean=filename.split(' ')[5]
    #distvar=filename[:-4].split(' ')[6]


    counts=np.sum(full)
    new=full/counts
    
    #print(np.sum(new))
    #print(new)
    #xmax,xmin=full.max(),full.min()
    #full=(full-xmin)/(xmax-xmin)
    
    df=pd.DataFrame(data=new,index=np.array(range(0,360)),columns=np.array(range(0,160)))
    
    
    dft=df.transpose()       
    
    ylist=list()
    xlist=list() 
    for i in range(20,180):
        if i==0:
            ylist.append(i)
        if i>0:
            ylist.append(i/20)
    
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = Axes3D(fig, azim = -128, elev = 43)
    x = dft.columns
    y = dft.index
    X,Y = np.meshgrid(x,y)
    Z = dft
    ax.plot_surface(X, Y, Z, rstride = 1, cstride = 1, norm = LogNorm(), cmap = cm.jet)
    
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(X, Y, Z)
    plt.show()
    #print(xlist)
    #resplot.set(xticklabels=xlist)
    #resplot.set(xlabel="Angle (°), "+"Mean: "+anglemean+" Var: "+anglevar,ylabel="Distance (Å), "+"Mean: "+distmean+" Var: "+distvar,title="Total Cases: "+str(counts)+" Hydrophobicity: "+hyd) 
   # resname=filename.split(" ")[0]+" "+filename.split(" ")[1]
    #plt.savefig(resname+'.png',format='PNG') 
    #plt.close()
