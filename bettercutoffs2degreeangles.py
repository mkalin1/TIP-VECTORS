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
    
    hydrophobicity={'ALA':0.20,'CYS':4.10,'ASP':-3.1,'GLU':-1.8,'PHE':4.40,'GLY':0.00,'HIS':0.50,'ILE':4.80,'LYS':-3.1,                      #dict containing hydrophobicity value of each residue
    'LEU':5.70,'MET':4.20,'ASN':-0.5,'PRO':-2.2,'GLN':-2.8,'ARG':1.40,'SER':-0.5,'THR':-1.9,'VAL':4.70,'TRP':1.00,'TYR':3.20}

    mol=molecule.load_structure(filename)
    
    coordObj_atoms=[i.get_parent().get_tip() for i in mol[0].get_calpha()]
    #coordObj_atoms=[i for i in mol[0].get_calpha()]
    #print(coordObj_atoms)
    coordObj=[i.get_parent().get_tip().get_location() for i in mol[0].get_calpha()]
    #coordObj=[i.get_location() for i in mol[0].get_calpha()]
    #print(coordObj)
    distance_mat = squareform(pdist(coordObj))
    n = len(distance_mat)
    interacting_residue_pairs = np.where( distance_mat <= cutoff_val )
    
    store = dict()
    storevec=dict()
    cluster_types = dict()
    retDict = dict()
    
    for i, j in zip(interacting_residue_pairs[0], interacting_residue_pairs[1]):
        try:
            store[i]
            storevec[coordObj_atoms[i]]
        except KeyError:
            store[i] = list()
            storevec[coordObj_atoms[i]]=list()
        

        store[i].append(coordObj_atoms[j].get_parent().get_name())                          #dict of cluster w/ resnames
        if i>=j:
            storevec[coordObj_atoms[i]].append([(coordObj_atoms[j].get_parent().get_calpha().get_location()-coordObj_atoms[j].get_location()),coordObj_atoms[j].get_parent(),coordObj_atoms[j].get_location()])      #stores key=central res (i): values(calpha-tip vctors, molobj (j), tiplocation(j))
       
        #print(coordObj_atoms[j].get_parent().get_calpha().get_location(),coordObj_atoms[j].get_parent().get_name())
        #print(coordObj_atoms[j].get_location(),coordObj_atoms[j].get_parent().get_name())
    
    ##########################################################################
    
    
    rescount=[]
    for i in store.keys():                       
        rescount.append(len(store[i]))                               # number of atoms in cluster (used for hydrophobicity% calculation)
        

    clusterhyd=dict()
    clusterhyd = {k: [hydrophobicity.get(v, v) for v in v] for k, v in store.items()}                        #hydrophobicity of cluster (list)
    for key,val in clusterhyd.items():                                                                        #sum hyd of cluster list
        clusterhyd[key]=sum(clusterhyd[key])
    
    hydperc=[]  
    for key,val in clusterhyd.items():
        hydperc.append(clusterhyd[key]/(5.70*rescount[key]))                                             # % hydrophobicity
    

   
    hyddict=dict(zip(storevec.keys(),hydperc))
    
    
   
 #############################################################################################################            
    angles=dict()
    

    for i in storevec.keys():                                                                           
        if i.get_parent().get_name()!='GLY' and i.get_parent().get_name()!='ALA':                   # create new dict (angles) which removes clusters in which gly/ala is central atom
            try:
                angles[i]
                
            except KeyError:
                
                angles[i]=list()
    
   
    
    
    for i in angles.keys():   
        #print(i)                                                                                                  
        for j in range(0,len(storevec[i])):  
            
            if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA':                                     # skips over gly/ala within a cluster
                
                distmin=np.linalg.norm(i.get_location()-storevec[i][j][2])
                
                vec1=i.get_parent().get_calpha().get_location()-i.get_location()
                vec2=storevec[i][j][0]
                if i.get_parent()==storevec[i][j][1]:                                                                           #no angle calculation if comparing center to itself
                    angle=0.0
                else:
                    #angle=vg.signed_angle(vec1,vec2, look=vg.basis.neg_z)
                    angle=vg.angle(vec1,vec2,look=None,units="deg")                              #angle calculation between i central res and j cluster res
                    
                #if angle < 0:
                    #angle = angle+360

                hydclust=hyddict[i]                                                                                                #get cluster information for each i central atom
                
                angles[i].append([(angle),(distmin),(hydclust),(i.get_parent().get_name(),storevec[i][j][1].get_name())])           #stores angle,distmin,hydrophobicity of cluster, and respair in dict angles
                
    #print(angles)
 #############################################################################################################            
             
    this=[]
    fdict=dict()
        
    for i in angles.keys():                                                                 #create sorted dictionary with respairs as keys
        for j in range(0,len(angles[i])):
            key=sorted(angles[i][j][3])
            this=angles[i][j]    
            #print(this)
            try:
                fdict[ key[0]+'-'+key[1] ].append(this)
            except:
                fdict[ key[0]+'-'+key[1] ]=[]
                fdict[ key[0]+'-'+key[1] ].append(this)
    
   
    
    print(k)
    for i in fdict.keys():                                                                  #append nhyd with cluster hydrophobicity data
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                try:
                    nhyd[i].append(fdict[i][j][2])
                except:
                    nhyd[i]=list()
                    nhyd[i].append(fdict[i][j][2])                                             

    
    for i in fdict.keys():                                                                      # append fdict2 with fdict data
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                fdict2[i].append(fdict[i][j])
  
            
    return fdict2,nhyd

###########################################################################################################################################################



with open('newnewnamestring.txt', 'r') as f:              #txt of all pdb file names to download
    namestring=f.read().split(",")
    


leng=len(namestring)     
names=[]
for i in range(0,leng):
    names.append(namestring[i])                          

#fdict2 is a inter-PDBfile dictionary with [respairs] as keys that contains all [angle,distance,hydrophbicity] data 
fdict2={'ARG-ARG':list(), 'ARG-ASN':list(), 'ARG-ASP':list(), 'ARG-CYS':list(), 'ARG-GLN':list(), 'ARG-GLU':list(), 'ARG-HIS':list(), 'ARG-ILE':list(), 'ARG-LEU':list(), 'ARG-LYS':list(), 'ARG-MET':list(), 'ARG-PHE':list(), 'ARG-PRO':list(), 'ARG-SER':list(), 'ARG-THR':list(), 'ARG-TRP':list(), 'ARG-TYR':list(), 'ARG-VAL':list(), 'ASN-ASN':list(), 'ASN-ASP':list(), 'ASN-CYS':list(), 'ASN-GLN':list(), 'ASN-GLU':list(), 'ASN-HIS':list(), 'ASN-ILE':list(), 'ASN-LEU':list(), 'ASN-LYS':list(), 'ASN-MET':list(), 'ASN-PHE':list(), 'ASN-PRO':list(), 'ASN-SER':list(), 'ASN-THR':list(), 'ASN-TRP':list(), 'ASN-TYR':list(), 'ASN-VAL':list(), 'ASP-ASP':list(), 'ASP-CYS':list(), 'ASP-GLN':list(), 'ASP-GLU':list(), 'ASP-HIS':list(), 'ASP-ILE':list(), 'ASP-LEU':list(), 'ASP-LYS':list(), 'ASP-MET':list(), 'ASP-PHE':list(), 'ASP-PRO':list(), 'ASP-SER':list(), 'ASP-THR':list(), 'ASP-TRP':list(), 'ASP-TYR':list(), 'ASP-VAL':list(), 'CYS-CYS':list(), 'CYS-GLN':list(), 'CYS-GLU':list(), 'CYS-HIS':list(), 'CYS-ILE':list(), 'CYS-LEU':list(), 'CYS-LYS':list(), 'CYS-MET':list(), 'CYS-PHE':list(), 'CYS-PRO':list(), 'CYS-SER':list(), 'CYS-THR':list(), 'CYS-TRP':list(), 'CYS-TYR':list(), 'CYS-VAL':list(), 'GLN-GLN':list(), 'GLN-GLU':list(), 'GLN-HIS':list(), 'GLN-ILE':list(), 'GLN-LEU':list(), 'GLN-LYS':list(), 'GLN-MET':list(), 'GLN-PHE':list(), 'GLN-PRO':list(), 'GLN-SER':list(), 'GLN-THR':list(), 'GLN-TRP':list(), 'GLN-TYR':list(), 'GLN-VAL':list(), 'GLU-GLU':list(), 'GLU-HIS':list(), 'GLU-ILE':list(), 'GLU-LEU':list(), 'GLU-LYS':list(), 'GLU-MET':list(), 'GLU-PHE':list(), 'GLU-PRO':list(), 'GLU-SER':list(), 'GLU-THR':list(), 'GLU-TRP':list(), 'GLU-TYR':list(), 'GLU-VAL':list(), 'HIS-HIS':list(), 'HIS-ILE':list(), 'HIS-LEU':list(), 'HIS-LYS':list(), 'HIS-MET':list(), 'HIS-PHE':list(), 'HIS-PRO':list(), 'HIS-SER':list(), 'HIS-THR':list(), 'HIS-TRP':list(), 'HIS-TYR':list(), 'HIS-VAL':list(), 'ILE-ILE':list(), 'ILE-LEU':list(), 'ILE-LYS':list(), 'ILE-MET':list(), 'ILE-PHE':list(), 'ILE-PRO':list(), 'ILE-SER':list(), 'ILE-THR':list(), 'ILE-TRP':list(), 'ILE-TYR':list(), 'ILE-VAL':list(), 'LEU-LEU':list(), 'LEU-LYS':list(), 'LEU-MET':list(), 'LEU-PHE':list(), 'LEU-PRO':list(), 'LEU-SER':list(), 'LEU-THR':list(), 'LEU-TRP':list(), 'LEU-TYR':list(), 'LEU-VAL':list(), 'LYS-LYS':list(), 'LYS-MET':list(), 'LYS-PHE':list(), 'LYS-PRO':list(), 'LYS-SER':list(), 'LYS-THR':list(), 'LYS-TRP':list(), 'LYS-TYR':list(), 'LYS-VAL':list(), 'MET-MET':list(), 'MET-PHE':list(), 'MET-PRO':list(), 'MET-SER':list(), 'MET-THR':list(), 'MET-TRP':list(), 'MET-TYR':list(), 'MET-VAL':list(), 'PHE-PHE':list(), 'PHE-PRO':list(), 'PHE-SER':list(), 'PHE-THR':list(), 'PHE-TRP':list(), 'PHE-TYR':list(), 'PHE-VAL':list(), 'PRO-PRO':list(), 'PRO-SER':list(), 'PRO-THR':list(), 'PRO-TRP':list(), 'PRO-TYR':list(), 'PRO-VAL':list(), 'SER-SER':list(), 'SER-THR':list(), 'SER-TRP':list(), 'SER-TYR':list(), 'SER-VAL':list(), 'THR-THR':list(), 'THR-TRP':list(), 'THR-TYR':list(), 'THR-VAL':list(), 'TRP-TRP':list(), 'TRP-TYR':list(), 'TRP-VAL':list(), 'TYR-TYR':list(), 'TYR-VAL':list(), 'VAL-VAL':list()}

nhyd=dict()

for i,k in enumerate(names[:3000]):

    print(i)
    #out,fdict2,nhyd = pull_clusters(k+'.pdb', 10.0, "A",fdict2,nhyd)
    try:   
    #fdict2,nhyd = pull_clusters(k+'.pdb', 12.0, "A",fdict2,nhyd)
        fdict2,nhyd = pull_clusters(k+'.pdb', 10.0, "A",fdict2,nhyd)  #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues
    except:
        continue





    
def hyd(nhyd,fdict2):                                              #hyd creates kmeans clusters from fdict and sorts the data into nested dictionary then returns it. newdict=dict with structure respair:centroids: ang/dist data                                                              #### plot hydrophobicity to confirm # of kclusters                                                                                                                                                                        #centroids: ang/dist data
    newdict=dict()
    newdict2=dict()
    for i in fdict2.keys():
        if i=="TRP-TYR"  or i== "ARG-TYR"  or i== "ASN-PHE"  or i== "CYS-GLU"  or i== "TRP-TRP"  or i== "HIS-PHE"  or i== "ASN-HIS"  or i== "HIS-MET"  or i== "ILE-MET"  or i== "LEU-MET"  or i== "CYS-TRP"  or i== "ASN-CYS"  or i== "ASP-CYS"  or i== "ASN-MET"  or i== "MET-TRP"  or i== "MET-TYR"  or i== "ARG-MET"  or i== "ARG-PHE"  or i== "ARG-TRP"  or i== "ARG-SER"  or i== "ARG-CYS"  or i== "SER-TYR"  or i== "ARG-LEU"  or i== "MET-MET"  or i== "MET-SER"  or i== "MET-PRO"  or i== "CYS-LEU"  or i== "SER-TRP"  or i== "HIS-HIS"  or i== "THR-TYR"  or i== "ARG-ILE"  or i== "CYS-SER"  or i== "CYS-LYS"  or i== "ARG-HIS"  or i== "ARG-ARG"  or i== "ASN-TYR"  or i== "CYS-MET"  or i== "PHE-SER"  or i== "PHE-PHE"  or i== "PRO-TRP"  or i== "ASN-TRP"  or i== "CYS-ILE"  or i== "CYS-PHE"  or i== "CYS-HIS"  or i== "CYS-PRO":
            k=1  
        if i== "PHE-TYR"  or i== "PHE-VAL"  or i== "GLN-TRP"  or i== "PRO-VAL"  or i== "GLU-PRO"  or i== "LYS-TYR"  or i== "TYR-TYR"  or i== "MET-THR"  or i== "ASP-PRO"  or i== "GLU-TYR"  or i== "LYS-PHE"  or i== "TYR-VAL"  or i== "LYS-VAL"  or i== "CYS-TYR"  or i== "CYS-THR"  or i== "ASP-ILE"  or i== "CYS-VAL"  or i== "ASP-SER"  or i== "ARG-GLN"  or i== "PHE-THR"  or i== "ASP-PHE"  or i== "ASP-MET"  or i== "ASP-THR"  or i== "HIS-TRP"  or i== "MET-PHE"  or i== "GLU-SER"  or i== "GLN-PHE"  or i== "LYS-MET"  or i== "GLN-TYR"  or i== "ARG-ASN"  or i== "GLN-ILE"  or i== "GLN-GLN"  or i== "ARG-GLU"  or i== "ARG-ASP"  or i== "GLN-LEU"  or i== "GLN-GLU"  or i== "GLN-VAL"  or i== "GLU-LYS"  or i== "GLU-GLU"  or i== "GLN-MET"  or i== "ARG-THR"  or i== "GLN-HIS"  or i== "ARG-LYS"  or i== "GLN-LYS"  or i== "CYS-CYS"  or i== "ARG-VAL"  or i== "GLU-LEU"  or i== "ARG-PRO"  or i== "GLU-ILE"  or i== "GLU-HIS"  or i== "GLU-TRP"  or i== "HIS-LYS"  or i== "ASP-LEU"  or i== "PRO-THR"  or i== "HIS-TYR"  or i== "LYS-THR"  or i== "ASP-LYS"  or i== "ASN-GLU"  or i== "ASN-LYS"  or i== "ILE-ILE"  or i== "HIS-PRO"  or i== "LEU-LYS"  or i== "ASP-VAL"  or i== "SER-THR"  or i== "LEU-LEU"  or i== "HIS-THR"  or i== "HIS-SER"  or i== "LEU-PHE"  or i== "ILE-VAL"  or i== "GLU-VAL"  or i== "SER-SER"  or i== "LYS-SER"  or i== "ILE-THR"  or i== "LEU-THR"  or i== "ILE-PHE"  or i== "LEU-TYR"  or i== "LEU-VAL"  or i== "HIS-LEU"  or i== "LEU-TRP"  or i== "ILE-PRO"  or i== "ILE-SER"  or i== "ILE-TRP"  or i== "ILE-LYS"  or i== "HIS-ILE"  or i== "ILE-LEU"  or i== "LYS-LYS"  or i== "ILE-TYR"  or i== "LEU-SER"  or i== "ASN-GLN"  or i== "ASP-GLU"  or i== "PHE-TRP"  or i== "MET-VAL"  or i== "ASP-TRP"  or i== "THR-TRP"  or i== "THR-THR"  or i== "SER-VAL"  or i== "ASN-PRO"  or i== "ASN-ASP"  or i== "PHE-PRO"  or i== "THR-VAL"  or i== "TRP-VAL"  or i== "ASN-SER"  or i== "ASN-THR"  or i== "ASN-LEU"  or i== "ASP-ASP"  or i== "ASP-HIS"  or i== "ASN-VAL"  or i== "ASN-ILE"  or i== "VAL-VAL"  or i== "LEU-PRO"  or i== "ASN-ASN"  or i== "ASP-TYR"  or i== "HIS-VAL"  or i== "ASP-GLN":
            k=2
        if i== "GLN-THR"  or i== "GLN-SER"  or i== "GLU-PHE"  or i== "LYS-TRP"  or i== "GLN-PRO"  or i== "PRO-PRO"  or i== "GLU-MET"  or i== "PRO-TYR"  or i== "CYS-GLN"  or i== "LYS-PRO"  or i== "PRO-SER"  or i== "GLU-THR":
            k=3
        newdict2[i]=defaultdict(list)
        newdict[i]=defaultdict(list)
        clusters, centroids = kmeans1d.cluster(nhyd[i], k)
        for j,z in zip(clusters,fdict2[i]):
            try:
                newdict2[i][j].append([z[0],z[1]])
            except KeyError:
                newdict2[i][j]=list()
                newdict2[i][j].append([z[0],z[1]])
        for n in centroids:
            try:
                newdict[i][n]
            except KeyError:
                newdict[i][n]=list()
        for b,v in enumerate(newdict[i]):
            newdict[i][v].append(newdict2[i][b])
            
    return newdict

'''
def stats(fdict2):
    for i in fdict2.keys():
        resangles=list()
        resdistance=list()
        reshyd=list()
        for j in range(0,len(fdict2[i])):
            resangles.append(fdict2[i][j][0])
            resdistance.append(fdict2[i][j][1])
            reshyd.append(fdict2[i][j][2])
    #print(i, " ", resangles, " ", resdistance, " ", reshyd)
'''



def writefile(fulldict):                                                                                    #writefile takes in dictionary from hyd and creates np angle/distance matrix file for each respair and centroid
    anglesdict=dict()
    distdict=dict()
    keydict=dict()
    #keydict2=dict()
    for i in fulldict.keys():
        for j in fulldict[i].keys():
            #print(i,j,fulldict[i][j][0][0][0],fulldict[i][j][0][0][2],fulldict[i][j][0][0][4])
            resangles=list()
            resdistance=list()
            distdict[i]=defaultdict(list)
            anglesdict[i]=defaultdict(list)
            
            for l in range(0,len(fulldict[i][j][0])):
                
                resangles.append(fulldict[i][j][0][l][0])
                resdistance.append(fulldict[i][j][0][l][1])
            try:
                anglesdict[i][j].append([np.mean(resangles),np.var(resangles)])
                distdict[i][j].append([np.mean(resdistance),np.var(resdistance)])
            except KeyError:
                anglesdict[i][j] = list()
                distdict[i][j]=list()
                anglesdict[i][j].append([np.mean(resangles),np.var(resangles)])
                distdict[i][j].append([np.mean(resdistance),np.var(resdistance)])
          

            for k in range(0,len(fulldict[i][j][0])):
                
                
                #print(i,j,fulldict[i][j][0][numk][0],fulldict[i][j][0][numk][2],fulldict[i][j][0][numk][4])
                try:
                    keydict[str(i)+ ' '+str(round(j,2))+' '+str(round(anglesdict[i][j][0][0],2))+' ' +str(round(anglesdict[i][j][0][1],2))+' '+str(round(distdict[i][j][0][0],2))+' '+str(round(distdict[i][j][0][1],2))+'.txt']
                except:

                    keydict[str(i)+ ' '+str(round(j,2))+' '+str(round(anglesdict[i][j][0][0],2))+' ' +str(round(anglesdict[i][j][0][1],2))+' '+str(round(distdict[i][j][0][0],2))+' '+str(round(distdict[i][j][0][1],2))+'.txt']=[list(),list()]
                
 
                x=fulldict[i][j][0][k][0]
                y=fulldict[i][j][0][k][1]
                
                keydict[str(i)+ ' '+str(round(j,2))+' '+str(round(anglesdict[i][j][0][0],2))+' ' +str(round(anglesdict[i][j][0][1],2))+' '+str(round(distdict[i][j][0][0],2))+' '+str(round(distdict[i][j][0][1],2))+'.txt'][0].append(x)
                keydict[str(i)+ ' '+str(round(j,2))+' '+str(round(anglesdict[i][j][0][0],2))+' ' +str(round(anglesdict[i][j][0][1],2))+' '+str(round(distdict[i][j][0][0],2))+' '+str(round(distdict[i][j][0][1],2))+'.txt'][1].append(y)

          
            #str(i)+ ' '+str(round(j,2))+' '+str(round(anglesdict[i][j][0],2))+' ' +str(round(anglesdict[i][j][1],2))+' '+str(round(distdict[i][j][0],2))+' '+str(round(distdict[i][j][1],2))+'.txt'

    for key,val in keydict.items():                                                               #### create and save matrix
        
        try:
            mtx=np.histogram2d(val[0],val[1],bins=(90, 20),range=[[0,180],[0,10]])
    
            
        
            np.savetxt(key,mtx[0],fmt='%i')
        except:
            print(val[0],val[1])
            continue
        

  
 
    return True


writefile(hyd(nhyd,fdict2))


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.
    
    
            