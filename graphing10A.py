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


files=glob.glob('*-*.txt')
for filename in files:
    full = np.loadtxt(filename)
    hyd=filename.split(' ')[2]
    anglemean=filename.split(' ')[3]
    anglevar=filename.split(' ')[4]
    distmean=filename.split(' ')[5]
    distvar=filename[:-4].split(' ')[6]


    counts=np.sum(full)
    new=full/counts
    
    print(np.sum(new))
    #print(new)
    
    df=pd.DataFrame(data=new,index=np.array(range(0,90)),columns=np.array(range(0,20)))
    
    
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
    
    #print(xlist)
    resplot.set(xticklabels=xlist)
    resplot.set(xlabel="Angle (°), "+"Mean: "+anglemean+" Var: "+anglevar,ylabel="Distance (Å), "+"Mean: "+distmean+" Var: "+distvar,title="Total Cases: "+str(counts)+" Hydrophobicity: "+hyd) 
          
    plt.savefig(filename[:-30]+'.png',format='PNG') 
    plt.close()
