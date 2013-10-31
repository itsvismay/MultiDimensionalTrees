# main.py

from ContactGeometry1 import *
import numpy as np
import matplotlib as plt

sig = 8 
FileName = 'PDBFiles/1QKU.pdb'

pdb = ReadPDB(FileName)
X = np.array(pdb[['x','y','z', 'serial']],dtype='float64')
X1 = X[:, 0:10]
C = ContactMatrix(X1,X1,sig)
R = ContactCoord(C)

# display contact matrix as heat map
#np.set_printoptions(precision=8,suppress=True)
#print(C.tolist())
#np.savetxt("contactmatrix.csv", C.tolist(),fmt='%.5f', delimiter=",")
#imgplot = plt.imshow(C)
#plt.colorbar()

