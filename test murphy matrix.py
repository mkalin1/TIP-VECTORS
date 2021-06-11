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

#_, p_b = scipy.stats.ttest_ind(df_a.dropna(axis=0), df_b.dropna(axis=0))
#_, p_c = scipy.stats.ttest_ind(df_a.dropna(axis=0), df_c.dropna(axis=0))
#pd.DataFrame([p_b, p_c], columns = df_a.columns, index = ['df_b', 'df_c'])

fdict2={'ARG-ARG':list(), 'ARG-ASN':list(), 'ARG-ASP':list(), 'ARG-CYS':list(), 'ARG-GLN':list(), 'ARG-GLU':list(), 'ARG-HIS':list(), 'ARG-ILE':list(), 'ARG-LEU':list(), 'ARG-LYS':list(), 'ARG-MET':list(), 'ARG-PHE':list(), 'ARG-PRO':list(), 'ARG-SER':list(), 'ARG-THR':list(), 'ARG-TRP':list(), 'ARG-TYR':list(), 'ARG-VAL':list(), 'ASN-ASN':list(), 'ASN-ASP':list(), 'ASN-CYS':list(), 'ASN-GLN':list(), 'ASN-GLU':list(), 'ASN-HIS':list(), 'ASN-ILE':list(), 'ASN-LEU':list(), 'ASN-LYS':list(), 'ASN-MET':list(), 'ASN-PHE':list(), 'ASN-PRO':list(), 'ASN-SER':list(), 'ASN-THR':list(), 'ASN-TRP':list(), 'ASN-TYR':list(), 'ASN-VAL':list(), 'ASP-ASP':list(), 'ASP-CYS':list(), 'ASP-GLN':list(), 'ASP-GLU':list(), 'ASP-HIS':list(), 'ASP-ILE':list(), 'ASP-LEU':list(), 'ASP-LYS':list(), 'ASP-MET':list(), 'ASP-PHE':list(), 'ASP-PRO':list(), 'ASP-SER':list(), 'ASP-THR':list(), 'ASP-TRP':list(), 'ASP-TYR':list(), 'ASP-VAL':list(), 'CYS-CYS':list(), 'CYS-GLN':list(), 'CYS-GLU':list(), 'CYS-HIS':list(), 'CYS-ILE':list(), 'CYS-LEU':list(), 'CYS-LYS':list(), 'CYS-MET':list(), 'CYS-PHE':list(), 'CYS-PRO':list(), 'CYS-SER':list(), 'CYS-THR':list(), 'CYS-TRP':list(), 'CYS-TYR':list(), 'CYS-VAL':list(), 'GLN-GLN':list(), 'GLN-GLU':list(), 'GLN-HIS':list(), 'GLN-ILE':list(), 'GLN-LEU':list(), 'GLN-LYS':list(), 'GLN-MET':list(), 'GLN-PHE':list(), 'GLN-PRO':list(), 'GLN-SER':list(), 'GLN-THR':list(), 'GLN-TRP':list(), 'GLN-TYR':list(), 'GLN-VAL':list(), 'GLU-GLU':list(), 'GLU-HIS':list(), 'GLU-ILE':list(), 'GLU-LEU':list(), 'GLU-LYS':list(), 'GLU-MET':list(), 'GLU-PHE':list(), 'GLU-PRO':list(), 'GLU-SER':list(), 'GLU-THR':list(), 'GLU-TRP':list(), 'GLU-TYR':list(), 'GLU-VAL':list(), 'HIS-HIS':list(), 'HIS-ILE':list(), 'HIS-LEU':list(), 'HIS-LYS':list(), 'HIS-MET':list(), 'HIS-PHE':list(), 'HIS-PRO':list(), 'HIS-SER':list(), 'HIS-THR':list(), 'HIS-TRP':list(), 'HIS-TYR':list(), 'HIS-VAL':list(), 'ILE-ILE':list(), 'ILE-LEU':list(), 'ILE-LYS':list(), 'ILE-MET':list(), 'ILE-PHE':list(), 'ILE-PRO':list(), 'ILE-SER':list(), 'ILE-THR':list(), 'ILE-TRP':list(), 'ILE-TYR':list(), 'ILE-VAL':list(), 'LEU-LEU':list(), 'LEU-LYS':list(), 'LEU-MET':list(), 'LEU-PHE':list(), 'LEU-PRO':list(), 'LEU-SER':list(), 'LEU-THR':list(), 'LEU-TRP':list(), 'LEU-TYR':list(), 'LEU-VAL':list(), 'LYS-LYS':list(), 'LYS-MET':list(), 'LYS-PHE':list(), 'LYS-PRO':list(), 'LYS-SER':list(), 'LYS-THR':list(), 'LYS-TRP':list(), 'LYS-TYR':list(), 'LYS-VAL':list(), 'MET-MET':list(), 'MET-PHE':list(), 'MET-PRO':list(), 'MET-SER':list(), 'MET-THR':list(), 'MET-TRP':list(), 'MET-TYR':list(), 'MET-VAL':list(), 'PHE-PHE':list(), 'PHE-PRO':list(), 'PHE-SER':list(), 'PHE-THR':list(), 'PHE-TRP':list(), 'PHE-TYR':list(), 'PHE-VAL':list(), 'PRO-PRO':list(), 'PRO-SER':list(), 'PRO-THR':list(), 'PRO-TRP':list(), 'PRO-TYR':list(), 'PRO-VAL':list(), 'SER-SER':list(), 'SER-THR':list(), 'SER-TRP':list(), 'SER-TYR':list(), 'SER-VAL':list(), 'THR-THR':list(), 'THR-TRP':list(), 'THR-TYR':list(), 'THR-VAL':list(), 'TRP-TRP':list(), 'TRP-TYR':list(), 'TRP-VAL':list(), 'TYR-TYR':list(), 'TYR-VAL':list(), 'VAL-VAL':list()}
nameslist=[]
for i in fdict2.keys():
    for j in range(1,11):
        nameslist.append(str(i)+' '+str(j))

files=glob.glob('*-*.txt')
names=[]


#df = pd.DataFrame(columns=nameslist, index=nameslist)
#df = pd.DataFrame(columns=nameslist, index=nameslist)
#print(df)


for filename in files:
    names.append(filename)
#print(names)
sqmat=[]
'''

'''
df=pd.read_csv('distancematrixMAX.csv',index_col=[0])



df1=np.triu(df)                               # upper triangular


'''
for i in range(2,30):
    A=sklearn.cluster.AgglomerativeClustering(n_clusters=i, affinity='euclidean', memory=None, connectivity=None, compute_full_tree='auto', linkage='ward', distance_threshold=None).fit_predict(df1)
    chscore=sklearn.metrics.calinski_harabasz_score(df, A)

    print(i,chscore)
'''
#A=sklearn.cluster.AgglomerativeClustering(n_clusters=14, affinity='euclidean', memory=None, connectivity=None, compute_full_tree='auto', linkage='ward', distance_threshold=None).fit_predict(df)
#chscore=sklearn.metrics.calinski_harabasz_score(df, A)


#print(A)



df=(df.assign(sum=df.sum(axis=0)).sort_values(by='sum', ascending=False).iloc[:, :-1])
#print(df)
fuck=[]
x=-1
z={}
for i in df.index:
    x=x+1
    y=-1
    for j in df.columns:
        y=y+1
        if df[i][j]==0:
            
            z[y]=j

z1=[]

for k,numk in enumerate(z.items()):
    #print(k,numk[0],numk[1])
    
    if 'VAL' in numk[1] or 'ILE' in numk[1] or 'MET' in numk[1] or 'THR' in numk[1] or 'PHE' in numk[1]:
        print(k,numk)
        continue
    else:
       
        z1.append(numk[0])

print(z1)
df2=df[df.columns[[z1]]]
#print(df2)
#print(df2)
#for i in z1.values():
    #newdf.append(df[i],ignore_index=True)
#df2.to_csv('AA.csv',index = True)

            #continue

#print(df.iloc[0][0])

#df=(df.iloc[169][150]).reindex

            
df2.to_csv('A.csv',index = True)

#for i in df.index:
    
    #print(df.loc[df[i]==0])


#data = df.sort_values(by=list(df.index),axis=1)

#fuck.to_csv('AAA.csv',index = True)
#np.savetxt('wtf',fuck,fmt='%f')