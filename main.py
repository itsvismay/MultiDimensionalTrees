# main.py

from ContactGeometry1 import *
import numpy as np

sig = 8 
FileName1 = 'd1cida1.ent'
FileName2 = 'd2rhea_.ent'

# read PDB files and extract CA atom coordinates
pdb1 = ReadPDB(FileName1)
pdb2 = ReadPDB(FileName2)
X1 = CA_coordinates(pdb1) 
X2 = CA_coordinates(pdb2) 

# compute contact matrices and contact coordinates
C1 = ContactMatrix(X1,X1,sig)
C2 = ContactMatrix(X2,X2,sig)
R1 = ContactCoord(C1)
R2 = ContactCoord(C2)


