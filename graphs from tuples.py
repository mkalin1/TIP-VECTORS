import numpy as np
import glob
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
#files=glob.glob('*-*.txt')
data1=np.loadtxt('CYS-CYS nobins.txt')
#df=pd.DataFrame(data=data1)


x=[]
y=[]
for line in data1:
    x.append(line[0])
    y.append(line[1])
avg=np.median(y)

print(len(y))


#data_list = [[row[0][0][0][0], row[0][0][1][0], row[0][1]] for row in data1]
#resplot=sns.heatmap(data)
#plt.show()
'''
for filename in files:
    for line in filename:
        print('yes')
'''