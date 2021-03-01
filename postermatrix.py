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
import scipy
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import rgb2hex, colorConverter
from scipy.cluster.hierarchy import set_link_color_palette
import scipy.cluster.hierarchy as sch
import sklearn
from matplotlib import cm 
import seaborn as sns

df=pd.read_csv('ascending.csv',index_col=[0])
#df1=np.triu(df) 
dft=df.transpose()
labs=dft.columns
labs2=dft.index
df1=np.triu(df) 
df1=df1.transpose()   
'''
fig, ax = plt.subplots()
ax.set_xticklabels(labs)
print(labs)

plt.imshow(df1, cmap='hot', interpolation='None')
plt.show()
'''

plt.figure(figsize=(46, 32))
resplot=sns.heatmap(df1,cmap='hot_r',annot=False,cbar=False,fmt='', square=True,linewidths=0.5)

num_elements = len(labs)
xticklist = []
xlablist=[]
yticklist=[]
ylablist=[]
for item in range (0,num_elements):
    xticklist.append(item+0.5)
    xlablist.append(labs[item])

for item in range (0,num_elements):
    yticklist.append(item+0.5)
    ylablist.append(labs2[item])




plt.xticks(ticks=xticklist,labels=xlablist, rotation=90,fontsize=10)
plt.yticks(ticks=yticklist,labels=ylablist, rotation=0,fontsize=10)

#plt.show()
plt.savefig('t4.png',format='PNG',transparent=True) 

#plt.show()
#np.savetxt('AAA.csv',dft,'%s', ',')

'''
hmap=sns.heatmap(df1)
hmap.figure.savefig("Correlation_Heatmap_with_Seaborn.png",format='png',dpi=171)


fig = plt.figure()
ax1 = fig.add_subplot(111)
cmap = cm.get_cmap('gray_r', 10)
ax1.imshow(df1, interpolation="nearest", cmap=cmap)
ax1.grid(True)
plt.show()
'''