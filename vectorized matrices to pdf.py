from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import numpy 
import matplotlib.pyplot as plt
from sklearn.neighbors import radius_neighbors_graph
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csgraph
from sklearn.cluster import SpectralClustering
import glob
from scipy.stats import gaussian_kde
import pandas as pd
import os
files=glob.glob('*-*.txt')
df=pd.DataFrame()

def FileList(files,df):
    for filename in files:
        resname=filename.split(" ")[0]
        print(resname)
        #newarr=[]
        fdict={}
        
        binned = numpy.loadtxt(filename,dtype=int)[:, 20:180]
        counts=numpy.sum(binned)
        new=binned/counts
        x=[]
        y=[]
        z=[]
        for i in range(0,360):
            for j in range(0,160):
                angle=i/2
                distance=1+(j/20)
                intensity=new[i][j]
                fdict[(angle,distance)]=intensity
        
                #newarr2.append([(angle,distance,intensity)])
        frame = pd.DataFrame.from_dict(fdict, orient='index',columns=[resname])
        
        df = pd.concat([df,frame],axis=1)
        
    return df
    
            
df=FileList(files,df)

df.to_csv('ALL PAIR TYPES.csv',index = True)



    
    #numpy.savetxt(resname+' vectorized.txt', newarr2, fmt='%1.2f') 

