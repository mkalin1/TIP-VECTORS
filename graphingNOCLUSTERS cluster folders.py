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

namesdict={1: ['ARG-GLU', 'ASP-LYS'], 2: ['ARG-ASP', 'LYS-LYS', 'LYS-TYR'], 3: ['GLU-LYS'], 4: ['PRO-PRO'], 5: ['ARG-PRO', 'GLN-SER', 'GLN-THR', 'GLU-PRO', 'GLU-SER', 'GLU-THR', 'HIS-SER', 'HIS-THR', 'MET-PRO', 'SER-TYR', 'TYR-VAL'], 6: ['ARG-ILE', 'ARG-LEU', 'ARG-MET', 
'ARG-PHE', 'ARG-SER', 'ARG-THR', 'ASN-HIS', 'ASN-TRP', 'ASN-TYR', 'ASP-TRP', 'GLN-MET', 'GLN-PHE', 'GLU-MET', 'GLU-PHE', 'HIS-MET', 'HIS-PHE', 'ILE-LYS', 'LEU-LYS', 'LYS-MET', 'LYS-PHE', 'LYS-PRO', 'LYS-SER', 'LYS-THR', 'SER-TRP', 'THR-TRP', 'THR-TYR'], 7: ['CYS-PHE', 'ILE-MET', 'ILE-PHE', 'ILE-TRP', 'LEU-MET', 'MET-VAL', 'PHE-VAL', 'PRO-TRP'], 8: ['ARG-TRP', 'ASN-SER', 'ASN-THR', 'ASP-SER', 'ASP-THR', 'GLN-PRO', 'GLN-TRP', 'GLU-TRP', 'HIS-TRP', 'ILE-TYR', 'LEU-TYR', 'MET-TYR', 'PHE-PRO', 'PHE-TYR', 'SER-SER', 'SER-THR', 'TRP-TYR', 'TRP-VAL'], 9: ['CYS-HIS', 'CYS-ILE', 'CYS-LEU', 'CYS-MET', 'CYS-VAL', 'ILE-VAL', 'LEU-VAL'], 10: ['ASN-CYS', 'ASN-VAL', 'ASP-CYS', 'ASP-VAL', 'CYS-SER', 'CYS-THR', 'ILE-SER', 'ILE-THR', 'LEU-SER', 'LEU-THR', 'SER-VAL'], 11: ['ASN-ILE', 'ASN-LEU', 'ASP-ILE', 'ASP-LEU', 'CYS-GLN', 'CYS-GLU', 'GLU-VAL'], 12: ['ARG-CYS', 'ARG-VAL', 'ASN-MET', 'ASN-PHE', 'ASP-MET', 'ASP-PHE', 'CYS-LYS', 'CYS-TYR', 'GLN-ILE', 'GLN-LEU', 'GLN-VAL', 'GLU-ILE', 'GLU-LEU', 'HIS-ILE', 'HIS-LEU', 'HIS-VAL', 'LYS-VAL', 'MET-SER', 'MET-THR', 'PHE-SER', 'PHE-THR'], 13: ['ARG-LYS', 'ASN-LYS', 'GLN-GLU', 'GLN-LYS', 'GLU-GLU'], 14: ['ARG-TYR', 'ASN-ASN', 'ASN-ASP', 'ASP-ASP', 'ASP-HIS', 'HIS-HIS', 'HIS-TYR', 'LYS-TRP'], 15: ['ARG-ARG', 'ARG-ASN', 'ARG-GLN', 'ARG-HIS', 'ASN-GLN', 'ASN-GLU', 
'ASP-GLN', 'ASP-GLU', 'ASP-TYR', 'GLN-GLN', 'GLN-HIS', 'GLN-TYR', 'GLU-HIS', 'GLU-TYR', 'HIS-LYS'], 16: ['CYS-TRP', 'ILE-ILE', 'ILE-LEU', 'LEU-LEU', 'MET-TRP', 'PHE-PHE', 'PHE-TRP', 'TRP-TRP'], 17: ['THR-VAL', 'VAL-VAL'], 18: ['CYS-PRO', 'ILE-PRO', 'LEU-PRO', 'PRO-SER', 'PRO-THR', 'PRO-VAL', 'THR-THR'], 19: ['LEU-PHE', 'LEU-TRP', 'MET-MET', 'MET-PHE'], 20: ['ASN-PRO', 'ASP-PRO', 'HIS-PRO', 'PRO-TYR', 'TYR-TYR'], 21: ['CYS-CYS']}
try:
    for key in namesdict.keys():
        os.mkdir(os.path.join('./', str(key)))
except:
    pass
files=glob.glob('*-*.txt')
for filename in files:
    full = np.loadtxt(filename)
    resname=filename.split(" ")[0]
    for key in namesdict.keys():
        for value in namesdict[key]:
            if value in resname:
                save_results_to = "/Users/Michael/Desktop/graph stuff/prontoCOMBINEDGRAPHS/"+str(key)+'/'+str(resname)+'.png'
                #print(save_results_to)
                #print(key,value,resname)
                counts=np.sum(full)
                new=full/counts
    




    #print(np.sum(new))
    #print(new)
    #xmax,xmin=full.max(),full.min()
    #full=(full-xmin)/(xmax-xmin)
    
                df=pd.DataFrame(data=new,index=np.array(range(0,90)),columns=np.array(range(0,20)))
    
    
                dft=df.transpose()       
                
                ylist=list()
                xlist=list() 
                for i in range(0,20):
                    if i==0:
                        ylist.append(i)
                    if i>0:
                        ylist.append(i/2)
                
                fig, ax = plt.subplots()
                
                resplot=sns.heatmap(dft,yticklabels=ylist)
                xlist=map(int,resplot.get_xticks()*2)

                #print(xlist)
                resplot.set(xticklabels=xlist)
                resplot.set(xlabel="Angle (°)",ylabel="Distance (Å)",title="Total Cases: "+str(counts)) 
                resname=filename.split(" ")[0]
                
                plt.savefig(save_results_to) 
                plt.close()


'''
for key in namesdict.keys():
    os.mkdir(os.path.join('./', str(key)))           ### create directory from dictionary key


def make_dirs_from_dict(d, current_dir='./'):
    for key, val in d.items():
        os.mkdir(os.path.join('./', key))
        if type(val) == dict:
            make_dirs_from_dict(val, os.path.join(current_dir, key))
'''