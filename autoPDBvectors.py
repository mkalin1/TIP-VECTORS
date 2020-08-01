from packman import molecule
import vg
import numpy as np
import os.path 
from os import path
import math
from sklearn import preprocessing
import sklearn 
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

from scipy import stats


with open('namestring.txt', 'r') as f:              #txt of all pdb file names to download
    namestring=f.read().split(",")
    


leng=len(namestring)
names=[]
for i in range(0,leng):
    names.append(namestring[i])
'''
for i in range(0,leng):
    if path.exists(namestring[i]+'.pdb'):
        continue 
    else:
        print(namestring[i])
        molecule.download_structure(namestring[i],(namestring[i]+'.pdb'))
'''


#molecule.download_structure('4KXV','4KXV.pdb')
#mol=molecule.load_structure('4KXV.pdb')


for k in names:                                              
    mol = molecule.load_structure(k+'.pdb')
    restricted_res = ['ALA', 'GLY']
    ang_dist = ()
    tl = []
    name = []
    Data = {}
    try:
        for i in mol[0].get_residues():
            for j in mol[0].get_residues():          

                if(i.get_id() < j.get_id() and i.get_name() not in restricted_res and j.get_name() not in restricted_res and 12 >= np.linalg.norm(i.get_tip().get_location()-j.get_tip().get_location())):

                    res1_vec = i.get_tip().get_location() - i.get_calpha().get_location() 
                    res2_vec = j.get_tip().get_location() - j.get_calpha().get_location()  
                    #name.append(i.get_name() + '-' + j.get_name())
                    #name1=i.get_name() + '-' + j.get_name()
                    # angle=np.arccos(np.dot(res1_vec,res2_vec)/(np.linalg.norm(res1_vec)*np.linalg.norm(res2_vec)))
                    angle = vg.signed_angle(res1_vec, res2_vec, look=vg.basis.z)
                    if angle < 0:
                        angle = angle+360
                    
                    res1 = i.get_tip().get_location()
                    #print(res1)
                    res2 = j.get_tip().get_location()
                    #print(res2)
                    
                    distmin = np.linalg.norm(res1-res2)

                    # if distmin > 13:
                    #print(i.get_id(), j.get_id(), distmin)
                    key = sorted([i.get_name(), j.get_name()])

                    ang_dist = (angle, distmin)
                    tl.append(ang_dist)

                    try:
                        Data[key[0]+'-'+key[1]].append(ang_dist)
                    except:
                        Data[key[0]+'-'+key[1]] = []
                        Data[key[0]+'-'+key[1]].append(ang_dist)
                    # print(ang_dist)
        print(k)
        from matplotlib import pyplot as plt
        import seaborn as sns
        import pandas as pd

        from scipy import stats

        addition = np.zeros(shape=(360, 36))
        for i in Data.keys():
            x = [j[0] for j in Data[i]]
            y = [j[1] for j in Data[i]]
            mtx = np.histogram2d(x, y, bins=(360, 36), range=[[0, 360], [0, 12]])

           

            if path.exists(i):
                full = np.loadtxt(i)
                addition = np.add(mtx[0], full)
            else:
                addition = mtx[0]

            '''
            for m in range(0,len(yint)):
                if any(n>100 for n in yint[m]):
                    counter=counter+1
            '''
            np.savetxt(i, addition, fmt='%i')
            # np.savetxt(i,mtx[0],fmt='%i')

            # exit()
    except:
        pass        
'''
with open('OUTLIERS.txt', 'w') as f:
  f.write('%i' % counter)

'''

#for j in Data.keys():
    #full = np.loadtxt(j+'tip')                    #enable
    #print(full)
    #y=()
    #y1=[]
    #x=[]
    #for numi,i in enumerate(full):
    #df1=sklearn.preprocessing.normalize(full)
    #print(df)
    #df=pd.DataFrame(data=df1,index=np.array(range(0,72)),columns=np.array(range(0,15)))
'''
    full1=np.array(full)
    f=0
    farray=[]
    for i in range(0,360):
        for j in range(0,30):
            full[i][j]
            f=f+full[i][j]
        farray.append(f/30)
        f=0
'''
    #print(farray)
    #full1=np.delete(full1, np.s_[14:30], axis=1)
    #full1 = np.delete(full1, [14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30], axis=1)
    #full1.resize(360,30)
    #df=pd.DataFrame(data=full,index=np.array(range(0,360)),columns=np.array(range(0,30)))               #enable this
    #df = pd.DataFrame(full, columns=["Angle", "Distance"])
    #print(df)
    
    #dft=df.transpose()                                                 #enable this
    #print(dft)
'''   
      
    
    sns.heatmap(farray)
    plt.savefig('averagedist'+j+'.png',format='PNG')
    plt.close()
'''
    #plt.imshow(df, cmap='hot', interpolation='nearest')
    #plt.show()
    #sns.heatmap(dft)                       #enable this
    
'''
    trials=str(len(df))
    xstd=str(round(np.std(df["Angle"]),2))
    ystd=str(round(np.std(df["Distance"]),2))
    xmean=str(round(df["Angle"].mean(),2))
    ymean=str(round(df["Distance"].mean(),2))
    xvar=str(round(np.var(df["Angle"]),2))
    yvar=str(round(np.var(df["Distance"]),2))
'''

    #g=sns.JointGrid(x="Angle", y="Distance", data=df)#, kind="scatter",marginal_kws=False)
    #g = g.plot_joint(sns.scatterplot)
    #g.fig.suptitle('TOTAL: ' +trials + '\n XMEAN: ' + xmean + ', XVAR: ' + xvar + ', XSTD: ' +xstd + '\n YMEAN: '+ ymean+ ', YVAR: ' + yvar +', YSTD: ' +ystd +'\n')


    
    #plt.savefig('tipNEW'+j+'.png',format='PNG')            #enable
    #plt.close()                                                    #enable
#dot vs angle instead of distance

    