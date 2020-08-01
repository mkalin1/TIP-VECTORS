from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from operator import add



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
    store_center = dict()
    cluster_types = dict()
    retDict = dict()

    

    for i, j in zip(interacting_residue_pairs[0], interacting_residue_pairs[1]):
        try:
            store[i]
            storetip[i]
            storecalpha[i]
        except KeyError:
            store[i] = list()
            storetip[i]=list()
            storecalpha[i]=list()
            cluster_types[i] = list()



        store[i].append(coordObj_atoms[j].get_parent().get_name())                          #dict of cluster w/ resnames
        storetip[i].append(coordObj_atoms[j].get_parent().get_tip().get_location())            #dict of cluster w/tip coords
        storecalpha[i].append(coordObj_atoms[j].get_location())           #dict of cluster w/calpha coords
        
        
        
        #store[i].append(np.array([coordObj[j][0], coordObj[j][1], coordObj[j][2]]))

        if coordObj_atoms[j].get_parent().get_name() in HYD:
            cluster_types[i].append("HB")
        else:
            cluster_types[i].append("HP")
    
    print(storetip)
    vec=[]
    c=[]
    for key,value in storetip.items():
        #c.append([[i[0]-j[0]] for i,j in zip(storecalpha[key],storetip[key])])
        c=np.subtract(storecalpha[key],storetip[key])
        #print(c)
    
        
    
        
    ##########################################################################
    
    #for i in storetip.keys():                       
        #print(mol[0].get_residues()[i].get_tip().get_location())          s
                                    
        
    rescount=[]
    
    for i in store.keys():                       
        rescount.append(len(store[i]))                               # number of atoms in cluster
    
    ###########################################################################
   
    clusterhyd=dict()

    clusterhyd = {k: [hydrophobicity.get(v, v) for v in v] for k, v in store.items()}              #hydrophobicity of cluster (list)


    for key,val in clusterhyd.items():                             #sum hyd of cluster list
        clusterhyd[key]=sum(clusterhyd[key])
    
    
    hydperc=[]  
    ##########################################################################

    for key,val in clusterhyd.items():
        hydperc.append(clusterhyd[key]/(5.70*rescount[key]))                                 # % hydrophobicity
    
    
    hyddict=dict(zip(hydperc,store.values()))
   
    ##########################################################################

    for key, val in store.items():
        try:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
        except:
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))] = list()
            retDict[str(cluster_types[key].count("HB"))+"/"+str(cluster_types[key].count("HP"))].append(val)
            
    return retDict
    


out = pull_clusters("5lxw.pdb", 7.0, "A") #Here are all your clusters with ids -number of hydrophobic residue/number of hydrophilic residues


# This is for one PDB id. You can collect this for many pdb ids and merge the clusters.



