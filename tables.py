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
from collections import defaultdict


files=glob.glob('*-*.txt')

def tablename(files):
    table=dict()
    for filename in files:
        
        full = np.loadtxt(filename)
        


        hyd=filename.split(' ')[2]
        anglemean=filename.split(' ')[3]
        anglevar=filename.split(' ')[4]
        distmean=filename.split(' ')[5]
        distvar=filename[:-4].split(' ')[6]
        
        try:
            table[filename[0: filename.find(' ')]][hyd]=[anglemean, anglevar, distmean, distvar]
        except KeyError:
            table[filename[0: filename.find(' ')]]=defaultdict(list)
            table[filename[0: filename.find(' ')]][hyd]=[anglemean, anglevar, distmean, distvar]
        
        #table[filename[0: filename.find(' ')]]={hyd:[anglemean, anglevar, distmean, distvar]}
       


    
    return table




#print(tablename(files))
df=tablename(files)


for i in df.keys():
 
    for j in df[i].keys():
        for k in df[i][j]:
            
            reform = {(i, j): k for i, j in df.items() for j, k in df[i].items()}


#print(reform['LEU-LEU','hyd-0.858'][0])

x=pd.DataFrame.from_dict(reform)
xt=x.transpose()

xt.columns=['Angle Mean','Angle Variance','Distance Mean','Distance Variance']

xt.to_csv('FK.csv',index = True)

#tabledf=pd.DataFrame.from_dict(tablename(files))
#tabledf=tabledf.transpose() 
#print(tabledf)



