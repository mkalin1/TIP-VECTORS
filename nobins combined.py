import numpy as np
import glob
import re
files=glob.glob('*-*.txt')
resdict={}
for filename in files:
    try:
        resdict[filename.split(" ")[0]].append(str(filename))
        
    except:
        resdict[filename.split(" ")[0]]=list()
        resdict[filename.split(" ")[0]].append(str(filename))

angdict=dict()
distdict=dict()
newdict=dict()
for key in resdict.keys():
    angles=[]
    dists=[]
    for i in range(0,10):
        

        with open(resdict[key][i],'r') as f:
            txt = f.read()
        nums = re.findall(r'\[([^][]+)\]', txt)
        data=np.loadtxt(nums,delimiter=",")
        
        angles=data[0]
        dists=data[1]
        print(np.median(dists))
        for i,j in zip(angles,dists):
            t=(i,j)
            try:
                newdict[key].append(t)
            except:
                newdict[key]=list()
                newdict[key].append(t)


'''
for i in newdict.keys():
    print(i)
    np.savetxt(i+' nobins.txt',newdict[i])


    angdict[key]=angles
    distdict[key]=dists


    try:
        angdict[key].append(angles)
        distdict[key].append(dists)
    except:
        angdict[key]=list()
        distdict[key]=list()
        angdict[key].append(angles)
        distdict[key].append(dists)
'''

        
