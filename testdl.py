import os
import re
from packman import molecule
'''
f=open("Input_Step2.txt","r")
for line in f:
    ID=line[13:18].split("\t")[1]
    
    runchain=line.split("_")[1]
    print(ID,runchain)

def downloads(filename):
    names=[]
    IDChain=()
    with open("Input_Step2.txt","r") as f:              #txt of all pdb file names to download
        for line in f:
            line = line.strip('\n')
            ID=line[13:18].split("\t")[1]
            runchain=line.split("_")[1]
            #print(ID,runchain)
            IDChain=(ID,runchain)
            names.append(IDChain)
    return names
'''
names=[]
def downloads(filename):
    IDChain=()
    with open(filename,"r") as f:              #txt of all pdb file names to download
        for line in f:
            line = line.strip('\n')
            ID=line[13:18].split("\t")[1]
            runchain=line.split("_")[1]
            #print(ID,runchain)
            IDChain=(ID,runchain)
            names.append(IDChain)
    return names

downloads("Input_Step2.txt")

for i,k in enumerate(names[0:2]):
    print(i,k[0])
    molecule.download_structure(k[0],('STRUCTURE.pdb'))

    