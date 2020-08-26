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
    dft=df.transpose()         
    sns.heatmap(dft)                      
    plt.savefig(filename[:-4]+' 2 DEGREE BINS'+'.png',format='PNG') 
    plt.close()
