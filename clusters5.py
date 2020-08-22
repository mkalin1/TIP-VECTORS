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

def pull_clusters(filename, cutoff_val, chain_id):
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
    
    '''
    for i in storevec.keys():
        print(i.get_parent().get_name(),hyddict[i])
    '''
    fdict=dict()
    
    for i in angles.keys():
        for j in range(0,len(storevec[i])):  
            if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA':     #################               
                
                #print(k.get_parent().get_id())
                    #if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA' and i.get_name!='GLY' and i.get_name!='ALA':                           #use mol object instead of name???
                distmin=np.linalg.norm(i.get_parent().get_tip().get_location()-storetip[i][j][0])
                
                

                if i.get_parent()==storevec[i][j][1]:
                    angle=0.0
                else:
                    angle=vg.signed_angle(i.get_location()-i.get_parent().get_tip().get_location(),storevec[i][j][0], look=vg.basis.z)
                
                
                #print(i.get_parent().get_tip().get_location(),storevec[i][j][0])

                if angle < 0:
                    angle = angle+360

                hydclust=hyddict[i]
                
                angles[i].append([(angle),(distmin),(hydclust),(i.get_parent().get_name(),storevec[i][j][1].get_name())])
                #for n in range(0,len(angles.items()))
                #get_parent().get_id()==get_parent().get_id()-n
                
                '''
                key=sorted(angles[i][j][3])
                print(angles[i][j][3])
                try:
                    fdict[ key[0]+'-'+key[1] ].append(angles[i])
                except:
                    fdict[ key[0]+'-'+key[1] ]=[]
                    fdict[ key[0]+'-'+key[1] ].append(angles[i])
                    '''
                        #print(distmin,storetip[i][j][1],storetip[i][k][1])
                        #print(vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z))
                        #print(storevec[i][j][1].get_name(),storevec[i][k][1].get_name())
                        #print(storevec[i][j+1<len(storevec[i])][0],storevec[i][j+1<len(storevec[i])][1])
                        #print(storevec[i][j][0],storevec[i][j][1])
                        #print(storevec[i][j][0])
                        #print("\n")
                        #print((vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z)))
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
    
 


    
    ninty=[]
    eighty=[]
    seventy=[]
    sixty=[]
    fifty=[]
    forty=[]
    thirty=[]
    twenty=[]
    ten=[]
    zeros=[]
    nzeros=[]
    nten=[]
    ntwenty=[]
    nthirty=[]
    anglesdists=[]
    


    
    for i in fdict.keys():
        for j in range(0,len(fdict[i])):
            if fdict[i][j][0]==0.0:
                continue
            else:
                if 0.9<fdict[i][j][2]<1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    ninty.append(anglesdists)
                if 0.8<fdict[i][j][2]<0.9:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    eighty.append(anglesdists)
                if 0.7<fdict[i][j][2]<0.8:  
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    seventy.append(anglesdists)
                if 0.6<fdict[i][j][2]<0.7:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    sixty.append(anglesdists)
                if 0.5<fdict[i][j][2]<0.6:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    fifty.append(anglesdists)
                if 0.4<fdict[i][j][2]<0.5:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    forty.append(anglesdists)
                    
                if 0.3<fdict[i][j][2]<0.4:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    thirty.append(anglesdists)
                if 0.2<fdict[i][j][2]<0.3:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    twenty.append(anglesdists)
                if 0.1<fdict[i][j][2]<0.2:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    ten.append(anglesdists)
                if 0<fdict[i][j][2]<0.1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    zeros.append(anglesdists)
                if -0.1<fdict[i][j][2]<0:
                    
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    nzeros.append(anglesdists)
                if -0.2<fdict[i][j][2]<-0.1:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    nten.append(anglesdists)
                if -0.3<fdict[i][j][2]<-0.2:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    ntwenty.append(anglesdists)
                if -0.4<fdict[i][j][2]<-0.3:
                    anglesdists=(fdict[i][j][0],fdict[i][j][1],i)
                    nthirty.append(anglesdists)
    
    

    fulldict={'ninty':ninty,'eighty':eighty,'seventy':seventy,'sixty':sixty,'fifty':fifty,'forty':forty,'thirty':thirty,'twenty':twenty,'ten':ten,'zeros':zeros,'nzeros':nzeros,'nten':nten,'ntwenty':ntwenty,'nthirty':nthirty}
    
    
    addition = np.zeros(shape=(360, 24))

    for i in fulldict.keys():
        for j in fulldict[i]:
            x=[j[0]]
            y=[j[1]]
            z=j[2]
            
            mtx=np.histogram2d(x,y,bins=(360, 24),range=[[0,360],[0,12]])
            

            if path.exists(i + ' ' + z):
                full = np.loadtxt(i + ' ' + z)
                addition=np.add(mtx[0],full)
            else:
                addition=mtx[0]
            
        
            np.savetxt(i + ' ' + z,addition,fmt='%i')
    
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
            
    return retDict



with open('namestring.txt', 'r') as f:              #txt of all pdb file names to download
    namestring=f.read().split(",")
    


leng=len(namestring)
names=[]
for i in range(0,leng):
    names.append(namestring[i])
    
for k in names:   
    out = pull_clusters(k+'.pdb', 12.0, "A") #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.



