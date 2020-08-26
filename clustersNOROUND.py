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

def pull_clusters(filename, cutoff_val, chain_id,fulldict):
    
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
    

    for i in fdict.keys():
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                if 0.9<fdict[i][j][2]<1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    
                    fulldict['ninty'].append(anglesdists)
                if 0.8<fdict[i][j][2]<0.9:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['eighty'].append(anglesdists)
                if 0.7<fdict[i][j][2]<0.8:  
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['seventy'].append(anglesdists)
                if 0.6<fdict[i][j][2]<0.7:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['sixty'].append(anglesdists)
                if 0.5<fdict[i][j][2]<0.6:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['fifty'].append(anglesdists)
                if 0.4<fdict[i][j][2]<0.5:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['forty'].append(anglesdists)
                    
                if 0.3<fdict[i][j][2]<0.4:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['thirty'].append(anglesdists)
                if 0.2<fdict[i][j][2]<0.3:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['twenty'].append(anglesdists)
                if 0.1<fdict[i][j][2]<0.2:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['ten'].append(anglesdists)
                if 0<fdict[i][j][2]<0.1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['zeros'].append(anglesdists)
                if -0.1<fdict[i][j][2]<0:
                    
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['nzeros'].append(anglesdists)
                if -0.2<fdict[i][j][2]<-0.1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['nten'].append(anglesdists)
                if -0.3<fdict[i][j][2]<-0.2:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['ntwenty'].append(anglesdists)
                if -0.4<fdict[i][j][2]<-0.3:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fulldict['nthirty'].append(anglesdists)

    
    
    

    
    
  
        #for i in angles.keys():
            #print(i.get_parent().get_id())
        #print(angles)
    
    print(k)
    #############################################################################################################

    for key, val in store.items():
        try:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
        except:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))] = list()
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
            
    return retDict,fulldict

def writefile(fulldict):
    
    addition = np.zeros(shape=(360, 24))
    
    keydict=dict()
    #keydict2=dict()
    

    for i in fulldict.keys():
        for j in fulldict[i]:
            z=j[2]
            try:
                keydict[i+ ' '+z+'.txt']
            except:

                keydict[i+ ' '+z+'.txt']=[list(),list()]
            
            x=j[0]
            y=j[1]
            
            keydict[i+ ' '+z+'.txt'][0].append(x)
            keydict[i+ ' '+z+'.txt'][1].append(y)
    '''
    for key,val in keydict.items():
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
            
   
    for key,val in keydict.items():                   #### keydict2 for no repititions within a hydcluster
        try:
            mtx=np.histogram2d(val[0],val[1],bins=(180, 24),range=[[0,360],[0,12]])
    
            
        
            np.savetxt(key,mtx[0],fmt='%i')
        except:
            print(val[0],val[1])
            continue
    
    return True


with open('newnewnamestring.txt', 'r') as f:              #txt of all pdb file names to download
    namestring=f.read().split(",")
    


leng=len(namestring)
names=[]
for i in range(0,leng):
    names.append(namestring[i])

fulldict={'ninty':list(),'eighty':list(),'seventy':list(),'sixty':list(),'fifty':list(),'forty':list(),'thirty':list(),'twenty':list(),'ten':list(),'zeros':list(),'nzeros':list(),'nten':list(),'ntwenty':list(),'nthirty':list()}
#fulldict={'0.9':list(),'0.8':list(),'0.7':list(),'0.6':list(),'0.5':list(),'0.4':list(),'0.3':list(),'0.2':list(),'0.1':list(),'-0.0':list(),'-0.1':list(),'-0.2':list(),'-0.3':list(),'0.4':list()}

for i,k in enumerate(names):

    print(i)
    #out,fulldict = pull_clusters(k+'.pdb', 12.0, "A", fulldict)
    try:   
        out,fulldict = pull_clusters(k+'.pdb', 7.0, "A", fulldict)  #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues
    except:
        continue

writefile(fulldict)


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.