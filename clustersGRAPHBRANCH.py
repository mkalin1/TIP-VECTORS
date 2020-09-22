from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from operator import add
from itertools import combinations
import vg
from matplotlib import pyplot as plt
import seaborn as sns
import os.path 
from os import path
import math
import timeit
import glob
import pandas as pd


files=glob.glob('*-*.txt')
for filename in files:
    full = np.loadtxt(filename)
    
    df=pd.DataFrame(data=full,index=np.array(range(0,180)),columns=np.array(range(0,24)))
    counts=np.sum(full)
    
    dft=df.transpose()         
    ylist=list()
    for i in range(0,24):
        if i==0:
            ylist.append(i)
        if i>0:
            ylist.append(i/2)
    
    resplot=sns.heatmap(dft,yticklabels=ylist)
    resplot.set(xlabel="Angle",ylabel="Distance",title="Total: "+str(counts)) 
          
    plt.savefig(filename[:-4]+'.png',format='PNG') 
    plt.close()
