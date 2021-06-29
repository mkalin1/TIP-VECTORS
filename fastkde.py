
from fastkde import fastkde
import numpy 
import pandas as pd
import glob
import pylab as PP

files=glob.glob('LEU-VAL nobins.txt')
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


    delta_x = (numpy.max(HVE_X) - numpy.min(HVE_X))/10
    delta_y = (numpy.max(HVE_Y) - numpy.min(HVE_Y))/10
    x_min = numpy.min(HVE_X) - delta_x
    y_min = numpy.min(HVE_Y) - delta_y
    x_max = numpy.max(HVE_X) + delta_x
    y_max = numpy.max(HVE_Y) + delta_y
    #xx, yy = numpy.mgrid[x_min:x_max:50j, y_min:y_max:50j]
    #positions = numpy.vstack([xx.ravel(), yy.ravel()])
    #values = numpy.vstack([HVE_X, HVE_Y])
    pOfYGivenX,axes = fastKDE.conditional(y,x)
    fig,axs = PP.subplots(1,2,figsize=(10,5))
    axs[0].plot(x,y,'k.',alpha=0.1)
    axs[0].set_title('Original (x,y) data')
    axs[1].contourf(axes[0],axes[1],pOfYGivenX,64)
#Overplot the original underlying relationship
    
    

    #Set axis limits to be the same
    


    fig.tight_layout()

    PP.show()
