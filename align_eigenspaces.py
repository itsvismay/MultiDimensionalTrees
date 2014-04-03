# align_eigenspaces.py

import numpy as np
from numpy.linalg import det,eig,svd
from ContactGeometry1 import *
from scipy.optimize import fmin_bfgs

def RMSD(Xa,Xb,Ia,Ib): 
	dX = Xa[:,Ia] - Xb[:,Ib] 
	return np.sqrt((dX*dX).mean())

def Superimpose(Xa,Xb,Ia,Ib):
	XaI = Xa[:,Ia]
	XbI = Xb[:,Ib]
	XaI_centroid = np.reshape(XaI.mean(axis=1),(3,1)) 
	XbI_centroid = np.reshape(XbI.mean(axis=1),(3,1)) 
	YaI = XaI - XaI_centroid
	YbI = XbI - XbI_centroid
	C = np.dot(YaI,YbI.T)
	(V,D,W) = svd(C)
	E = np.eye(3)
	E[2,2] = np.sign(det(C))
	U = np.dot(np.dot(W,E),V.T) 
	Ya = Xa - XaI_centroid
	Yb = Xb - XbI_centroid
	Yar = np.dot(U,Ya)
	rmsd = RMSD(Yar,Yb,Ia,Ib)
	return Yar,Yb,rmsd

def Eig(M):
	# computes sorted eigensystem
	lam,V = eig(M)
	I = np.argsort(-lam)
	print("Lenght of V")
	print(len(V))
	if any(lam<0):
		print 'Eig --> warning: negative eigenvalue detected.'
	return (lam[I],V[:,I])

def ContactMatrix(Xa,Xb,sig):
	Na = Xa.shape[1]
	Nb = Xb.shape[1]
	D2 = np.zeros((Na,Nb))
	for i in np.arange(3):
		Da = np.tile(Xa[i,:],(Nb,1))
		Db = np.tile(Xb[i,:],(Na,1))
		D  = Da.T - Db
		D2 = D2 + D*D
	return np.exp(-0.5*(D2/sig**2)) 

def ContactCoord(C):
	[lam,V] = Eig(C)
	R = np.dot(np.diag(np.sqrt(lam)),V.T)
	return R

def Rot(o):
	v = np.sqrt(np.dot(o,o.T))
	n = o/v
	c = np.cos(v/2.0)
	s = np.sin(v/2.0)
	q0 = c 
	q1 = s*n[0]
	q2 = s*n[1]
	q3 = s*n[2]
	rot = np.zeros((3,3))
	
	rot[0,0] = q0**2.0 + q1**2.0 - q2**2.0 - q3**2.0
	rot[1,0] = 2.0*q1*q2 - 2.0*q0*q3
	rot[2,0] = 2.0*q1*q3 + 2.0*q0*q2

	rot[0,1] = 2.0*q1*q2 + 2.0*q0*q3
	rot[1,1] = q0**2.0 - q1**2.0 + q2**2.0 - q3**2.0
	rot[2,1] = 2.0*q2*q3 - 2.0*q0*q1

	rot[0,2] = 2.0*q1*q3 - 2.0*q0*q2
	rot[1,2] = 2.0*q2*q3 + 2.0*q0*q1
	rot[2,2] = q0**2.0 - q1**2.0 - q2**2.0 + q3**2.0

	return rot

def move(x,X):
	o = x[0:3]
	p = x[3:]
	Xr = np.dot(Rot(o),X) + reshape(p,(3,1)) 
	return Xr

def ContactScore(x,Xa,Xb,Va,Vb,sig):
	#(n1, n2, n3, px, py, pz)
	o = x[0:3]
	p = x[3:]
	Xar = np.dot(Rot(o),Xa) + np.reshape(p,(3,1)) 
	Cab = ContactMatrix(Xar,Xb,sig)
	scr = -(np.dot(Va.T,np.dot(Cab,Vb))).sum()
	return scr

def PrincipleAxesTransform(X):
	# center structure
	X_centroid = np.reshape(X.mean(axis=1),(3,1))
	Y = X - X_centroid 
	# rotate structure
	C = np.dot(Y,Y.T)
	(lam,V) = Eig(C)
	return np.dot(V.T,Y)

sig = 8 
tol = 0.001 
itermax = 100
FileName = 'PDBFiles/d1cida1.ent'
pdb = ReadPDB(FileName)
Xa = np.array(pdb[['x','y','z']],dtype='float64')

FileName = 'PDBFiles/d2rhea_.ent'
pdb2 = ReadPDB(FileName)
Xb = np.array(pdb[['x','y','z']],dtype='float64')
Xa = PrincipleAxesTransform(Xa[:10].T)
Xb = PrincipleAxesTransform(Xb[:10].T)
Ca = ContactMatrix(Xa,Xa,sig)
Cb = ContactMatrix(Xb,Xb,sig)
tempa= Eig(Ca)
la = tempa[0]
Va = tempa[1]
# #################################
tempb = Eig(Cb)
lb = tempb[0]
Vb = tempb[1]
x0 = np.array([0.001,0,0,0,0,0])

#####Following is the score calculation function
xopt,scr,A,B,a,b,c = fmin_bfgs(ContactScore,x0,args=(Xa,Xb,Va,Vb,sig),full_output=1,disp=1)
#xopt is the optimal n,p / scr is the optimal score
o = xopt[0:3]
p = xopt[3:]
v = np.sqrt(np.dot(o,o.T))
n = o/v
fid = open('rot_trans.txt','w')
s = '%f %f %f\n' % (n[0],n[1],n[2])
fid.write(s)
print 'axis of rotation'
print s
theta = (180.0/np.pi)*v
s = '%f\n' % theta 
fid.write(s)
print 'angle of rotation'
print s
s = '%f %f %f\n' % (p[0],p[1],p[2])
fid.write(s)
print 'displacement'
print s
fid.close()

# print 'contact score before'
# print ContactScore(x0,Xa,Xb,Va,Vb,sig)
# print 'contact score after'
# print ContactScore(xopt,Xa,Xb,Va,Vb,sig)
