from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import numpy 
import matplotlib.pyplot as plt
from sklearn.neighbors import radius_neighbors_graph
from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csgraph
from sklearn.cluster import SpectralClustering
import glob
from scipy.stats import gaussian_kde
import pandas as pd
files=glob.glob('CYS-CYS nobins.txt')
for filename in files:
    x, y = numpy.loadtxt(filename, unpack=True)
    
    HVE_X = x
    HVE_Y = y


   # for i in ResidueType_Entropy:
  #      for j in ResidueType_Entropy[i]:
   #         try:
  #              HVE_X.append(Hydrophobicity[i])
   #             HVE_Y.append(j)
    #        except:
   #             None


    delta_x = (numpy.max(HVE_X) - numpy.min(HVE_X))/ 10
    delta_y = (numpy.max(HVE_Y) - numpy.min(HVE_Y))/ 10
    x_min = numpy.min(HVE_X) - delta_x
    y_min = numpy.min(HVE_Y) - delta_y
    x_max = numpy.max(HVE_X) + delta_x
    y_max = numpy.max(HVE_Y) + delta_y
    xx, yy = numpy.mgrid[x_min:x_max:1000j, y_min:y_max:1000j]
    positions = numpy.vstack([xx.ravel(), yy.ravel()])
    values = numpy.vstack([HVE_X, HVE_Y])
    kernal = gaussian_kde(values)                              ### gives all peaks and all variances of peaks
    f = numpy.reshape(kernal(positions).T, xx.shape)
    plt.rcParams["font.family"] = "Arial"
    fig = plt.figure(figsize=(8,8))
    ax = fig.gca()
    #ax.set_xlim(xmin, xmax)
    #ax.set_ylim(ymin, ymax)
    cfset = ax.contour( yy, xx, f, cmap='Greys')
    ax.set_xlabel('distnace')
    ax.set_ylabel('angle')
    plt.tight_layout()
    plt.savefig(filename+'.png', dpi=1000)
    plt.close()
