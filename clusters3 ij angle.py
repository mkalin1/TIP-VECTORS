from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from operator import add
from itertools import combinations
import vg

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
            storetip[i]
            storecalpha[i]
            storevec[i]
        except KeyError:
            store[i] = list()
            storetip[i]=list()
            storecalpha[i]=list()
            storevec[i]=list()
            cluster_types[i] = list()

        store[i].append(coordObj_atoms[j].get_parent().get_name())                          #dict of cluster w/ resnames
        storevec[i].append([(coordObj_atoms[j].get_location()-coordObj_atoms[j].get_parent().get_tip().get_location()),coordObj_atoms[j].get_parent()])      #stores calpha-tip vctors with mol OBJECT
        #storevec[i].append([(coordObj_atoms[j].get_location()-coordObj_atoms[j].get_parent().get_tip().get_location()),coordObj_atoms[j].get_parent().get_name()]) #stores calpha-tip vctors with resname
        #storetip[i].append([coordObj_atoms[j].get_parent().get_tip().get_location(),coordObj_atoms[j].get_parent().get_name()])            #dict of cluster w/tip coords
        storetip[i].append([coordObj_atoms[j].get_parent().get_tip().get_location(),coordObj_atoms[j].get_parent().get_name()])            #dict of cluster w/tip coords and names

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
        
    
    hyddict=dict(zip(hydperc,store.values()))
    
    
       ###################################################################################################             creating dict of [key]:[[angle],[distmin],[(respair)]]
    angles=dict()
    for i in storevec.keys():
        try:
            angles[i]
        except KeyError:
            angles[i]=list()
        for j in range(0,len(storevec[i])):
            if storevec[i][j][1].get_name()!='GLY' and storevec[i][j][1].get_name()!='ALA':                                            #use mol object instead of name???
                for k in range(j,len(storevec[i])):
                    if storevec[i][k][1].get_name()!='GLY' and storevec[i][k][1].get_name()!='ALA':
                        if k>j:
                            distmin=np.linalg.norm(storetip[i][j][0]-storetip[i][k][0])
                            angle=vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z)
                            if angle < 0:
                                angle = angle+360
                            hydclust=hydperc[i]
                            angles[i].append([(angle),(distmin),(hydclust),(storevec[i][j][1].get_name(),storevec[i][k][1].get_name())])
                            #print(distmin,storetip[i][j][1],storetip[i][k][1])
                            #print(vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z))
                            #print(storevec[i][j][1].get_name(),storevec[i][k][1].get_name())
                #print(storevec[i][j+1<len(storevec[i])][0],storevec[i][j+1<len(storevec[i])][1])
                    #print(storevec[i][j][0],storevec[i][j][1])
                #print(storevec[i][j][0])
        #print("\n")
                    #print((vg.signed_angle(storevec[i][j][0],storevec[i][k][0], look=vg.basis.z)))

    for i in angles.keys():
        print(angles[i])

   ##############################################################################################################
    #hydtip=dict(zip(hydperc,angles.values()))
    #print(hydperc)




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



