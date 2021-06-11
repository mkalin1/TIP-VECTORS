import numpy as np
A=np.loadtxt('ARG-ARG COMBINED 616972.0.txt')
counts=np.max(A)
new=A/counts

C=np.concatenate((new,np.zeros((90,70))), axis=1)



#B=np.resize(A,(90,90))


vals, vecs = np.linalg.eig(C)


for i, value in enumerate(vals):
    print("Eigenvector:", vecs[:,i], ", Eigenvalue:", value)
