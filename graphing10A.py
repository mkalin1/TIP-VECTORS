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


files=glob.glob('*-*.txt')
for filename in files:
    full = np.loadtxt(filename)
    hyd=filename.split(' ')[1]
    anglemean=filename.split(' ')[2]
    anglevar=filename.split(' ')[3]
    distmean=filename.split(' ')[4]
    distvar=filename.split(' ')[5]

    df=pd.DataFrame(data=full,index=np.array(range(0,90)),columns=np.array(range(0,20)))
    counts=np.sum(full)
    
    dft=df.transpose()         
    ylist=list()
    xlist=list()
    for i in range(0,20):
        if i==0:
            ylist.append(i)
        if i>0:
            ylist.append(i/2)
    
    
    resplot=sns.heatmap(dft,yticklabels=ylist)
    xlist=map(int,resplot.get_xticks()*2)
    
    print(xlist)
    resplot.set(xticklabels=xlist)
    resplot.set(xlabel="Angle (°), "+"Mean: "+anglemean+" Var: "+anglevar,ylabel="Distance (Å), "+"Mean: "+distmean+" Var: "+distvar,title="Total: "+str(counts)+" Hydrophobicity: "+hyd) 
          
    plt.savefig(filename[:-27]+'.png',format='PNG') 
    plt.close()
