from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from operator import add
from itertools import combinations
import vg
from matplotlib import pyplot as plt
import seaborn as sns

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
        if i.get_parent().get_name()!='GLY' and i.get_parent().get_name()!='ALA':
            try:
                angles[i]
            except KeyError:
                angles[i]=list()
    
    for i in angles.keys():
        for j in range(0,len(storevec[i])):     
            if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA':
                
                #print(k.get_parent().get_id())
                    #if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA' and i.get_name!='GLY' and i.get_name!='ALA':                                            #use mol object instead of name???
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
                
                        #print(distmin,storetip[i][j][1],storetip[i][k][1])
                        #print(vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z))
                        #print(storevec[i][j][1].get_name(),storevec[i][k][1].get_name())
                        #print(storevec[i][j+1<len(storevec[i])][0],storevec[i][j+1<len(storevec[i])][1])
                        #print(storevec[i][j][0],storevec[i][j][1])
                        #print(storevec[i][j][0])
                        #print("\n")
                        #print((vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z)))
    
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
    
    anglesdists=[]
    fortydict=dict()
    
    for i in angles.keys():
        for j in range(0,len(angles[i])):
            if angles[i][j][0]==angles[i][j][1]:
                continue
            else:
                
                
                if 0.9<angles[i][j][2]<1:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    ninty.append(anglesdists)
                if 0.8<angles[i][j][2]<0.9:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    eighty.append(anglesdists)
                if 0.7<angles[i][j][2]<0.8:  
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    seventy.append(anglesdists)
                if 0.6<angles[i][j][2]<0.7:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    sixty.append(anglesdists)
                if 0.5<angles[i][j][2]<0.6:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    fifty.append(anglesdists)
                if 0.4<angles[i][j][2]<0.5:
                    #print(angles[i][j])
                    
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    key=sorted(angles[i][j][3])
                    
                    try:
                        fortydict[ key[0]+'-'+key[1] ].append(anglesdists)
                    except:
                        fortydict[ key[0]+'-'+key[1] ]=[]
                        fortydict[ key[0]+'-'+key[1] ].append(anglesdists)
                    
                if 0.3<angles[i][j][2]<0.4:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    thirty.append(anglesdists)
                if 0.2<angles[i][j][2]<0.3:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    

                if 0.1<angles[i][j][2]<0.2:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    ten.append(anglesdists)
                if 0<angles[i][j][2]<0.1:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    zeros.append(anglesdists)
                if -0.1<angles[i][j][2]<0:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    nzeros.append(anglesdists)
                if -0.2<angles[i][j][2]<-0.1:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    nten.append(anglesdists)
                if -0.3<angles[i][j][2]<-0.2:
                    anglesdists=(angles[i][j][0],angles[i][j][1])
                    ntwenty.append(anglesdists)
    print(fortydict)
    #print(twentydict)
    #sns.heatmap(twenty)
   # plt.show()
    #df = sns.load_dataset('iris')
    #sns.regplot(x=df["sepal_length"], y=df["sepal_width"])
    #for i in angles.keys():
        #print(i.get_parent().get_id())
    #print(angles)
    '''   
    for i in angles.keys():
        print(angles[i])
        print("\n")
    '''

    #############################################################################################################

    for key, val in store.items():
        try:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
        except:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))] = list()
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
            
    return retDict
    


out = pull_clusters("5lxw.pdb", 7.0, "A") #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.



