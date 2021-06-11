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

'''
full=np.loadtxt('GLU-LYS COMBINED 1858461.0.txt')

counts=np.sum(full)
new=full/counts
np.savetxt('GLU-LYS.txt',new,fmt='%1.9f')
'''


#print(np.sum(new))
#print(new)
#xmax,xmin=full.max(),full.min()
#full=(full-xmin)/(xmax-xmin)


full1 = np.loadtxt('GLU-GLU.txt')
full2=np.loadtxt('GLU-LYS.txt')
distance=np.sum(np.sqrt((full1-full2)**2))

print(distance)



