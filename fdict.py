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
import kmeans1d
from collections import defaultdict
def pull_clusters(filename, cutoff_val, chain_id,fdict2,nhyd):
    
    hydrophobicity={
    'ALA':0.20,
    'CYS':4.10,
    'ASP':-3.1,
    'GLU':-1.8,
    'PHE':4.40,
    'GLY':0.00,
    'HIS':0.50,
    'ILE':4.80,
    'LYS':-3.1,
    'LEU':5.70,
    'MET':4.20,
    'ASN':-0.5,
    'PRO':-2.2,
    'GLN':-2.8,
    'ARG':1.40,
    'SER':-0.5,
    'THR':-1.9,
    'VAL':4.70,
    'TRP':1.00,
    'TYR':3.20
    }


    
    HYD = ["GLY", "ALA", "VAL", "LEU", "ISO", "PRO", "PHE", "MET", "TRP"] # Edit this as per your requirement
    mol=molecule.load_structure(filename)
    coordObj_atoms=[i for i in mol[0].get_calpha()]
    coordObj=[i.get_location() for i in mol[0].get_calpha()]
    distance_mat = squareform(pdist(coordObj))
    n = len(distance_mat)
    interacting_residue_pairs = np.where( distance_mat <= cutoff_val )
    
    store = dict()
    storetip=dict()
    storecalpha=dict()
    storevec=dict()
    store_center = dict()
    cluster_types = dict()
    retDict = dict()
    
    for i, j in zip(interacting_residue_pairs[0], interacting_residue_pairs[1]):
        try:
            store[i]
            storetip[coordObj_atoms[i]]
            storecalpha[coordObj_atoms[i]]
            storevec[coordObj_atoms[i]]
        except KeyError:
            store[i] = list()
            storetip[coordObj_atoms[i]]=list()
            storecalpha[coordObj_atoms[i]]=list()
            storevec[coordObj_atoms[i]]=list()
            cluster_types[i] = list()
        
        
        store[i].append(coordObj_atoms[j].get_parent().get_name())                          #dict of cluster w/ resnames
        if i>=j:
            storevec[coordObj_atoms[i]].append([(coordObj_atoms[j].get_location()-coordObj_atoms[j].get_parent().get_tip().get_location()),coordObj_atoms[j].get_parent()])      #stores calpha-tip vctors with mol OBJECT
            #storevec[i].append([(coordObj_atoms[j].get_location()-coordObj_atoms[j].get_parent().get_tip().get_location()),coordObj_atoms[j].get_parent().get_name()]) #stores calpha-tip vctors with resname
            #storetip[i].append([coordObj_atoms[j].get_parent().get_tip().get_location(),coordObj_atoms[j].get_parent().get_name()])            #dict of cluster w/tip coords
            storetip[coordObj_atoms[i]].append([coordObj_atoms[j].get_parent().get_tip().get_location(),coordObj_atoms[j]])            #dict of cluster w/tip coords and names
        
        #storecalpha[i].append([coordObj_atoms[j].get_location(),coordObj_atoms[j].get_parent().get_name()])           #dict of cluster w/calpha coords
        #storetip[i].append([coordObj_atoms[j].get_parent().get_tip().get_location()])
        #storecalpha[i].append([coordObj_atoms[j].get_location()])
        
        #store[i].append(np.array([coordObj[j][0], coordObj[j][1], coordObj[j][2]]))

        if coordObj_atoms[j].get_parent().get_name() in HYD:
            cluster_types[i].append("HB")
        else:
            cluster_types[i].append("HP")
        
    
    
    ##########################################################################
    
    #for i in storetip.keys():                       
        #print(mol[0].get_residues()[i].get_tip().get_location())          s
    
    rescount=[]
    for i in store.keys():                       
        rescount.append(len(store[i]))                               # number of atoms in cluster
        
    
    ###########################################################################
    clusterhyd=dict()
    clusterhyd = {k: [hydrophobicity.get(v, v) for v in v] for k, v in store.items()}              #hydrophobicity of cluster (list)
    for key,val in clusterhyd.items():                                                              #sum hyd of cluster list
        clusterhyd[key]=sum(clusterhyd[key])
    ##########################################################################
    hydperc=[]  
    for key,val in clusterhyd.items():
        hydperc.append(clusterhyd[key]/(5.70*rescount[key]))                                 # % hydrophobicity
    

   
    hyddict=dict(zip(storetip.keys(),hydperc))
    
    
   
 #############################################################################################################             creating dict of [key]:[[angle],[distmin],[(respair)]]
    angles=dict()
   
    for i in storevec.keys():
        if i.get_parent().get_name()!='GLY' and i.get_parent().get_name()!='ALA':                   ########## deal with this later
            try:
                angles[i]
                
            except KeyError:
                
                angles[i]=list()
    
   
    fdict=dict()
    
  
#################   #################   #################   #################   #################   #################   

    for i in angles.keys():
        for j in range(0,len(storevec[i])):  
            if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA':                 
                
                distmin=np.linalg.norm(i.get_parent().get_tip().get_location()-storetip[i][j][0])

                if i.get_parent()==storevec[i][j][1]:
                    angle=0.0
                else:
                    angle=vg.signed_angle(i.get_location()-i.get_parent().get_tip().get_location(),storevec[i][j][0], look=vg.basis.z)
            
                if angle < 0:
                    angle = angle+360

                hydclust=hyddict[i]
                
                angles[i].append([(angle),(distmin),(hydclust),(i.get_parent().get_name(),storevec[i][j][1].get_name())])
                
                
#################   #################   #################   #################   #################   #################                   
                        
    this=[]
    fdict=dict()
        
    for i in angles.keys():
        for j in range(0,len(angles[i])):
            key=sorted(angles[i][j][3])
            this=angles[i][j]    
            #print(this)
            try:
                fdict[ key[0]+'-'+key[1] ].append(this)
            except:
                fdict[ key[0]+'-'+key[1] ]=[]
                fdict[ key[0]+'-'+key[1] ].append(this)
    
   
    
    anglesdists=[]
    print(k)
    for i in fdict.keys():
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                try:
                    nhyd[i].append(fdict[i][j][2])
                except:
                    nhyd[i]=list()
                    nhyd[i].append(fdict[i][j][2])

    
    for i in fdict.keys():
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                fdict2[i].append(fdict[i][j])
               


    #############################################################################################################

    for key, val in store.items():
        try:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
        except:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))] = list()
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
            
    return retDict,fdict2,nhyd




###########################################################################################################################################################



with open('newnewnamestring.txt', 'r') as f:              #txt of all pdb file names to download
    namestring=f.read().split(",")
    


leng=len(namestring)
names=[]
for i in range(0,leng):
    names.append(namestring[i])


fdict2={'ARG-ARG':list(), 'ARG-ASN':list(), 'ARG-ASP':list(), 'ARG-CYS':list(), 'ARG-GLN':list(), 'ARG-GLU':list(), 'ARG-HIS':list(), 'ARG-ILE':list(), 'ARG-LEU':list(), 'ARG-LYS':list(), 'ARG-MET':list(), 'ARG-PHE':list(), 'ARG-PRO':list(), 'ARG-SER':list(), 'ARG-THR':list(), 'ARG-TRP':list(), 'ARG-TYR':list(), 'ARG-VAL':list(), 'ASN-ASN':list(), 'ASN-ASP':list(), 'ASN-CYS':list(), 'ASN-GLN':list(), 'ASN-GLU':list(), 'ASN-HIS':list(), 'ASN-ILE':list(), 'ASN-LEU':list(), 'ASN-LYS':list(), 'ASN-MET':list(), 'ASN-PHE':list(), 'ASN-PRO':list(), 'ASN-SER':list(), 'ASN-THR':list(), 'ASN-TRP':list(), 'ASN-TYR':list(), 'ASN-VAL':list(), 'ASP-ASP':list(), 'ASP-CYS':list(), 'ASP-GLN':list(), 'ASP-GLU':list(), 'ASP-HIS':list(), 'ASP-ILE':list(), 'ASP-LEU':list(), 'ASP-LYS':list(), 'ASP-MET':list(), 'ASP-PHE':list(), 'ASP-PRO':list(), 'ASP-SER':list(), 'ASP-THR':list(), 'ASP-TRP':list(), 'ASP-TYR':list(), 'ASP-VAL':list(), 'CYS-CYS':list(), 'CYS-GLN':list(), 'CYS-GLU':list(), 'CYS-HIS':list(), 'CYS-ILE':list(), 'CYS-LEU':list(), 'CYS-LYS':list(), 'CYS-MET':list(), 'CYS-PHE':list(), 'CYS-PRO':list(), 'CYS-SER':list(), 'CYS-THR':list(), 'CYS-TRP':list(), 'CYS-TYR':list(), 'CYS-VAL':list(), 'GLN-GLN':list(), 'GLN-GLU':list(), 'GLN-HIS':list(), 'GLN-ILE':list(), 'GLN-LEU':list(), 'GLN-LYS':list(), 'GLN-MET':list(), 'GLN-PHE':list(), 'GLN-PRO':list(), 'GLN-SER':list(), 'GLN-THR':list(), 'GLN-TRP':list(), 'GLN-TYR':list(), 'GLN-VAL':list(), 'GLU-GLU':list(), 'GLU-HIS':list(), 'GLU-ILE':list(), 'GLU-LEU':list(), 'GLU-LYS':list(), 'GLU-MET':list(), 'GLU-PHE':list(), 'GLU-PRO':list(), 'GLU-SER':list(), 'GLU-THR':list(), 'GLU-TRP':list(), 'GLU-TYR':list(), 'GLU-VAL':list(), 'HIS-HIS':list(), 'HIS-ILE':list(), 'HIS-LEU':list(), 'HIS-LYS':list(), 'HIS-MET':list(), 'HIS-PHE':list(), 'HIS-PRO':list(), 'HIS-SER':list(), 'HIS-THR':list(), 'HIS-TRP':list(), 'HIS-TYR':list(), 'HIS-VAL':list(), 'ILE-ILE':list(), 'ILE-LEU':list(), 'ILE-LYS':list(), 'ILE-MET':list(), 'ILE-PHE':list(), 'ILE-PRO':list(), 'ILE-SER':list(), 'ILE-THR':list(), 'ILE-TRP':list(), 'ILE-TYR':list(), 'ILE-VAL':list(), 'LEU-LEU':list(), 'LEU-LYS':list(), 'LEU-MET':list(), 'LEU-PHE':list(), 'LEU-PRO':list(), 'LEU-SER':list(), 'LEU-THR':list(), 'LEU-TRP':list(), 'LEU-TYR':list(), 'LEU-VAL':list(), 'LYS-LYS':list(), 'LYS-MET':list(), 'LYS-PHE':list(), 'LYS-PRO':list(), 'LYS-SER':list(), 'LYS-THR':list(), 'LYS-TRP':list(), 'LYS-TYR':list(), 'LYS-VAL':list(), 'MET-MET':list(), 'MET-PHE':list(), 'MET-PRO':list(), 'MET-SER':list(), 'MET-THR':list(), 'MET-TRP':list(), 'MET-TYR':list(), 'MET-VAL':list(), 'PHE-PHE':list(), 'PHE-PRO':list(), 'PHE-SER':list(), 'PHE-THR':list(), 'PHE-TRP':list(), 'PHE-TYR':list(), 'PHE-VAL':list(), 'PRO-PRO':list(), 'PRO-SER':list(), 'PRO-THR':list(), 'PRO-TRP':list(), 'PRO-TYR':list(), 'PRO-VAL':list(), 'SER-SER':list(), 'SER-THR':list(), 'SER-TRP':list(), 'SER-TYR':list(), 'SER-VAL':list(), 'THR-THR':list(), 'THR-TRP':list(), 'THR-TYR':list(), 'THR-VAL':list(), 'TRP-TRP':list(), 'TRP-TYR':list(), 'TRP-VAL':list(), 'TYR-TYR':list(), 'TYR-VAL':list(), 'VAL-VAL':list()}

nhyd=dict()

for i,k in enumerate(names[:5000]):

    print(i)
    #out,fdict2,nhyd = pull_clusters(k+'.pdb', 10.0, "A",fdict2,nhyd)
    try:   
        out,fdict2,nhyd = pull_clusters(k+'.pdb', 12.0, "A",fdict2,nhyd)  #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues
    except:
        continue

def hyd(nhyd,fdict2):
    k=10
    newdict=dict()
    newdict2=dict()
    for i in fdict2.keys():
        newdict2[i]=defaultdict(list)
        newdict[i]=defaultdict(list)
        clusters, centroids = kmeans1d.cluster(nhyd[i], k)
        for j,z in zip(clusters,fdict2[i]):
            try:
                newdict2[i][j].append([z[0],z[1],z[2],j,i])
            except KeyError:
                
                newdict2[i][j]=list()
                newdict2[i][j].append([z[0],z[1],z[2],j,i])
                #newdict2[i].append({j:[z[0],z[1],z[2],j]})
        for n in centroids:
            try:
                newdict[i][n]
            except KeyError:
                newdict[i][n]=list()
                #newdict[i][n].append(list())
        #for j in range(0,newdict[i]):
            #print(newdict[i][j])
        for b,v in enumerate(newdict[i]):
            #print(newdict2[i][b])
            #if b==newdict2[i]
            #fuck=list(newdict2[i][b].values())
            #print(b,v)
        
            newdict[i][v].append(newdict2[i][b])
            
    return newdict
    
def writefile(fulldict):
    
    addition = np.zeros(shape=(360, 24))
    
    keydict=dict()
    #keydict2=dict()
    

    for i in fulldict.keys():
        for j in fulldict[i].keys():
            #print(i,j,fulldict[i][j][0][0][0],fulldict[i][j][0][0][2],fulldict[i][j][0][0][4])
            for numk,k in enumerate(fulldict[i][j][0]):
                #print(i,j,fulldict[i][j][0][numk][0],fulldict[i][j][0][numk][2],fulldict[i][j][0][numk][4])
                
                try:
                    keydict[str(i)+ ' '+str(j)+'.txt']
                except:

                    keydict[str(i)+ ' '+str(j)+'.txt']=[list(),list()]
                
                x=fulldict[i][j][0][numk][0]
                y=fulldict[i][j][0][numk][1]
                
                keydict[str(i)+ ' '+str(j)+'.txt'][0].append(x)
                keydict[str(i)+ ' '+str(j)+'.txt'][1].append(y)
    for key,val in keydict.items():                   #### keydict2 for no repititions within a hydcluster
        try:
            mtx=np.histogram2d(val[0],val[1],bins=(180, 24),range=[[0,360],[0,12]])
    
            
        
            np.savetxt(key,mtx[0],fmt='%i')
        except:
            print(val[0],val[1])
            continue
        

  
 
    return True
 
                
writefile(hyd(nhyd,fdict2))


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.
    
    
'''
            z=[2]
            try:
                keydict[i+ ' '+j+'.txt']
            except:

                keydict[i+ ' '+j+'.txt']=[list(),list()]
            
            x=j[0]
            y=j[1]
            
            keydict[i+ ' '+z+'.txt'][0].append(x)
            keydict[i+ ' '+z+'.txt'][1].append(y)
          
    for key,val in keydict.items():                   #### keydict2 for no repititions within a hydcluster
        try:
            mtx=np.histogram2d(val[0],val[1],bins=(180, 24),range=[[0,360],[0,12]])
    
            
        
            np.savetxt(key,mtx[0],fmt='%i')
        except:
            print(val[0],val[1])
            continue
    for key,val in keydict2.items():
        check=dict()
        for n,i in enumerate(val[0]):
            try:
                check[i]
                continue
            except:
                check[i]="NA"
                try:
                    keydict2[key][0].append(i)
                    keydict2[key][1].append(val[1][n])
                except:
                    keydict2[key]=[list(),list()]
                    keydict2[key][0].append(i)
                    keydict2[key][1].append(val[1][n])
'''
            