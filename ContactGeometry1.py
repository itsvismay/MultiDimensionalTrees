# ContactGeometry1.py

import pandas as pd
from numpy import arange,argsort,diag,dot,exp,sqrt,tile,zeros
from numpy.linalg import eig

def ReadPDB(FileName):
    # reads any PDB file and returns only ATOM records
    names = ['record','serial','name','altLoc','resName','chainID','resSeq','iCode','x','y','z','occupancy','tempFactor','element','charge']
    cols = [(0,6),(6,12),(12,16),(16,17),(17,21),(21,22),(22,26),(26,30),(30,38),(38,46),(46,54),(54,60),(60,76),(76,78),(78,80)]
    pdb0 = pd.read_fwf(FileName,colspecs=cols,header=None,names=names)
    pdb1 = pdb0[pdb0['record']=='ATOM']
    pdb = pdb1[pdb0['name'] =='CA']
    print FileName
    return pdb

def ContactMatrix(Xa,Xb,sig):
	Na = Xa.shape[1]
	Nb = Xb.shape[1]
	D2 = zeros((Na,Nb))
	for i in arange(3):
		Da = tile(Xa[i,:],(Nb,1))
		Db = tile(Xb[i,:],(Na,1))
		D  = Da.T - Db
		D2 = D2 + D*D
	return exp(-0.5*(D2/sig**2)) 

def Eig(M):
	# computes sorted eigensystem
	lam,V = eig(M)
	I = argsort(-lam)
	if any(lam<0):
		print 'Eig --> warning: negative eigenvalue detecpdted.'	
	return lam[I],V[I,:][:,I] 

def ContactCoord(C):
	[lam,V] = Eig(C)
	return dot(diag(sqrt(lam)),V.T)

