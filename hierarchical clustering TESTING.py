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


from collections import defaultdict
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.colors import rgb2hex, colorConverter
from scipy.cluster.hierarchy import set_link_color_palette
import pandas as pd
import scipy.cluster.hierarchy as sch


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
for i in names:
    
    matrixA=i.split(" ")[0]+' '+i.split(" ")[1]
    df1 = np.loadtxt(i)
    counts1=np.max(df1)
    new1=df1/counts1
    #print(new1)
    for j in names:
        
        matrixB=j.split(" ")[0]+' '+j.split(" ")[1]
        
        df2 = np.loadtxt(j)
        counts2=np.sum(df2)
        new2=df2/counts2
        distance=math.sqrt(np.sum((new1-new2)**2))
        df.at[matrixA,matrixB]=distance
'''
df=pd.read_csv('distancematrix.csv',index_col=[0])

df.head(1)
'''
dendrogram = sch.dendrogram(sch.linkage(df, method  = "ward"))
plt.title('Dendrogram')
plt.xlabel('Customers')
plt.ylabel('Euclidean distances')
plt.show()
'''
#c_dist = pdist(df) # computing the distance

#c_link = linkage(df,  metric='correlation', method='complete')# computing the linkage
#c_link=sch.linkage(df,method='ward', metric='euclidean', optimal_ordering=True)
#B=dendrogram(c_link,labels=list(df.columns))


#plt.show()

#df.to_csv('distancematrix.csv',index = True)
        
        
class Clusters(dict):
    def _repr_html_(self):
        html = '<table style="border: 0;">'
        for c in self:
            hx = rgb2hex(colorConverter.to_rgb(c))
            html += '<tr style="border: 0;">' \
            '<td style="background-color: {0}; ' \
                       'border: 0;">' \
            '<code style="background-color: {0};">'.format(hx)
            html += c + '</code></td>'
            html += '<td style="border: 0"><code>' 
            html += repr(self[c]) + '</code>'
            html += '</td></tr>'

        html += '</table>'

        return html
def get_cluster_classes(den, label='ivl'):
    cluster_idxs = defaultdict(list)
    for c, pi in zip(den['color_list'], den['icoord']):
        for leg in pi[1:3]:
            i = (leg - 5.0) / 10.0
            if abs(i - int(i)) < 1e-5:
                cluster_idxs[c].append(int(i))
    
    cluster_classes = Clusters()
    for c, l in cluster_idxs.items():
        i_l = [den[label][i] for i in l]
        cluster_classes[c] = i_l
    
    return cluster_classes
def get_clust_graph(df, numclust, transpose=False, dataname=None, save=False, xticksize=8):
    if transpose==True:
        aml=df.transpose()
        xl="x-axis"
    else:
        aml=df
        xl="y-axis"
    data_dist = pdist(aml.transpose()) # computing the distance
    data_link = linkage(data_dist,  metric='correlation', method='complete')#method="complete") # computing the linkage
    B=dendrogram(data_link,labels=list(aml.columns),p=numclust, truncate_mode="lastp",get_leaves=True, count_sort='ascending', show_contracted=True)
    #myInd = [i for i, c in zip(B['ivl'], B['color_list']) if c=='g']
    get_cluster_classes(B)
    ax=plt.gca()
    ax.tick_params(axis='x', which='major', labelsize=xticksize)
    ax.tick_params(axis='y', which='major', labelsize=15)
    plt.xlabel(xl)
    #plt.set_size_inches(18.5, 10.5)
    plt.ylabel('Distance')
    plt.suptitle(xl+" clustering for "+dataname, fontweight='bold', fontsize=16);
    if save:
        plt.savefig(str(df.index.name)+str(numclust)+"tr_"+str(transpose)+"dn_"+str(dataname)+save+'.png')
    else:
        print("Not saving")
    return get_cluster_classes(B)
    
def give_cluster_assigns(df, numclust, transpose=True):
    if transpose==True:
        data_dist = pdist(df.transpose())
        data_link = linkage(data_dist,  metric='correlation', method='complete')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.columns)
    else:
        data_dist = pdist(df)
        data_link = linkage(data_dist,  metric='correlation', method='complete')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.index)
    for i in range(1,numclust+1):
        print("Cluster ",str(i),": ( N =",len(cluster_assigns[cluster_assigns==i].index),")", ", ".join(list(cluster_assigns[cluster_assigns==i].index)))
def give_cluster_assigns(df, numclust, transpose=True):
    if transpose==True:
        data_dist = pdist(df.transpose())
        data_link = linkage(data_dist,  metric='correlation', method='complete')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.columns)
    else:
        data_dist = pdist(df)
        data_link = linkage(data_dist,  metric='correlation', method='complete')
        cluster_assigns=pd.Series(sch.fcluster(data_link, numclust, criterion='maxclust', monocrit=None), index=df.index)
    for i in range(1,numclust+1):
        print("Cluster ",str(i),": ( N =",len(cluster_assigns[cluster_assigns==i].index),")", ", ".join(list(cluster_assigns[cluster_assigns==i].index)))
    
get_clust_graph(df, 20, transpose=True,dataname="Residue Pairs", save="HIERARCHICAL CLUSTERING", xticksize=9)
give_cluster_assigns(df,20,transpose=True)
plt.show()