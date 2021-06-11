import cv2 as cv
from packman import molecule
from scipy.spatial.distance import pdist, squareform
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os.path 
from os import path
import math
import glob
import pandas as pd
from matplotlib.ticker import FuncFormatter
from sklearn.preprocessing import normalize

files=glob.glob('*-*.txt')
for filename in files[:1]:
    full = np.loadtxt(filename)
   
    counts=np.sum(full)
    new=full/counts
    f=plt.imshow(new)
    plt.show()
    
    #img = cv.imread(str(filename),0)
    kernel = np.ones((10,10),np.uint8)
    grad = cv.morphologyEx(img,cv.MORPH_GRADIENT,kernel)
    plt.imshow(grad,'gray')
    plt.show()
    #rosion.savefig(filename+'.png',format='PNG') 
    