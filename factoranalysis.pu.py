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
from factor_analyzer import FactorAnalyzer
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
from factor_analyzer.factor_analyzer import calculate_kmo
from factor_analyzer import FactorAnalyzer

files=glob.glob('CYS-CYS nobins.txt')
for filename in files:
    df = pd.read_csv(filename, delimiter = " ",header=None)
    print(df.info)
    
    chi_square_value,p_value=calculate_bartlett_sphericity(df)
    kmo_all,kmo_model=calculate_kmo(df)
    print(chi_square_value, p_value)
    print(kmo_model)
    fa = FactorAnalyzer()
    fa.fit(df)
    
    ev, v = fa.get_eigenvalues()
    print(fa.loadings_)