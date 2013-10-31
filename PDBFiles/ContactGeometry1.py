# ContactGeometry1.py

from numpy import arange,argmax,argsort,array,cos,diag,dot,eye,exp,insert,max,real,reshape,sign,sin,sqrt,sum,tile,zeros
from numpy.linalg import det,eig,svd
from scipy.optimize import fmin_ncg,fmin_bfgs
import pandas as pd

def ReadPDB(FileName):
    # reads PDB files and returns only ATOM records
    names = ['record','serial','name','altLoc','resName','chainID','resSeq','iCode','x','y','z','occupancy','tempFactor','element','charge']
    cols = [(0,6),(6,12),(12,16),(16,17),(17,21),(21,22),(22,26),(26,30),(30,38),(38,46),(46,54),(54,60),(60,76),(76,78),(78,80)]
    pdb0 = pd.read_fwf(FileName,colspecs=cols,header=None,names=names)
    pdb = pdb0[pdb0['record']=='ATOM']
    print FileName
    print pdb[0:10]
    print
    return pdb

def CA_coordinates(pdb):
	# reads coordinates of CA atoms
	return array(pdb[pdb.name=='CA'][['x','y','z']].values.T,dtype='float64')

def ContactMatrix(Xa,Xb,sig):
	# computes either a protein or a protein-protein contact matrix
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
		print 'Eig --> warning: negative eigenvalue detected.'	
	return lam[I],V[I,:][:,I] 

def ContactCoord(C):
	# computes contact coordinates of a protein
	[lam,V] = Eig(C)
	return dot(diag(sqrt(lam)),V.T)

def RMSD(Xa,Xb,Ia,Ib): 
	# computes RMSD of paired coordinates
	dX = Xa[:,Ia] - Xb[:,Ib] 
	return sqrt((dX*dX).mean())

def Superimpose(Xa,Xb,Ia,Ib):
	# superimposes two proteins and computes RMSD

	# center structures 
	XaI = Xa[:,Ia]
	XbI = Xb[:,Ib]
	XaI_centroid = reshape(XaI.mean(axis=1),(3,1)) 
	XbI_centroid = reshape(XbI.mean(axis=1),(3,1)) 
	YaI = XaI - XaI_centroid
	YbI = XbI - XbI_centroid
	# rotate first structure 
	C = dot(YaI,YbI.T)
	(V,D,W) = svd(C)
	E = eye(3)
	E[2,2] = sign(det(C))
	U = dot(dot(W,E),V.T) 
	Ya = Xa - XaI_centroid
	Yb = Xb - XbI_centroid
	Yar = dot(U,Ya)
	rmsd = RMSD(Yar,Yb,Ia,Ib)
	return Yar,Yb,rmsd

def DP(Score,gap): 
	# Dynamic Programming
	#     Score    - cost matrix
	#     gap      - gap penalty
	#     Ia,Ib    - optimal alignment
	#     maxScore - optimal alignment score

	# initialization 
	M,N = Score.shape
	S   = zeros((M+1,N+1))
	Idx = zeros((M+1,N+1))
	for j in arange(1,N+1):
		S[0,j]   = S[0,j-1] - gap
		Idx[0,j] = 0
	for i in arange(1,M+1): 
		S[i,0]   = S[i-1,0] - gap
		Idx[i,0] = 1;
	# forward sweep
	for i in arange(1,M+1):
		for j in arange(1,N+1):
			gap_a = S[i,j-1] - gap
			gap_b = S[i-1,j] - gap
			match = S[i-1,j-1] + Score[i-1,j-1]
			S[i,j]   = max([gap_a,gap_b,match])
			Idx[i,j] = argmax([gap_a,gap_b,match])
	maxScore = S[M,N]
	# traceback 
	Ia = array([],dtype='int64')
	Ib = array([],dtype='int64')
	i = M
	j = N
	while (i>0) or (j>0):
		if Idx[i,j]==2:
			Ia = insert(Ia,0,i-1)
			Ib = insert(Ib,0,j-1)
			i = i-1
			j = j-1
		elif Idx[i,j]==0:
			j = j-1
		elif Idx[i,j]==1:
			i = i-1
	return Ia,Ib,maxScore

def AlignEigenvalues(la,lb,gap):
	# aligns eigenspaces by aligning corresponding eigenvalues
	Na = len(la)
	Nb = len(lb)
	lA = tile(la,(Nb,1))
	lB = tile(lb,(Na,1))
	Score = -2*(abs(lA.T-lB)/(lA.T+lB))
	Ia,Ib,maxScore = DP(Score,gap)
	return Ia,Ib,maxScore 

def Rot(o):
	v = sqrt(dot(o,o.T))
	n = o/v
	c = cos(v/2.0)
	s = sin(v/2.0)
	q0 = c 
	q1 = s*n[0]
	q2 = s*n[1]
	q3 = s*n[2]
	rot = zeros((3,3))
	
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

def DRot(o):
	v1 = o[0]
	v2 = o[1]
	v3 = o[2]
	v = sqrt(dot(o,o.T))
	n = o/v
	c = cos(v/2.0)
	s = sin(v/2.0)
	a = c*v - 2.0*s
	q0 = c 
	q1 = s*n[0]
	q2 = s*n[1]
	q3 = s*n[2]

	drot_dq = zeros((4,3,3))
	drot_dq[0,:,:] = 2.0*array([[q0,q3,-q2],[-q3,q0,q1],[q2,-q1,q0]])
	drot_dq[1,:,:] = 2.0*array([[q1,q2,q3],[q2,-q1,q0],[q3,-q0,-q1]])
	drot_dq[2,:,:] = 2.0*array([[-q2,q1,-q0],[q1,q2,q3],[q0,q3,-q2]])
	drot_dq[3,:,:] = 2.0*array([[-q3,q0,q1],[-q0,-q3,q2],[q1,q2,q3]])

	dqdo = zeros((4,3))
	dqdo[0,:] = (-v**2.0*s)*array([v1,v2,v3])
	dqdo[1,:] = array([2.0*v**2.0*s+v1**2.0*a, v1*v2*a, v1*v3*a])
	dqdo[2,:] = array([v1*v2*a,2.0*v**2.0*s+v2**2.0*a,v2*v3*a])
	dqdo[3,:] = array([v1*v3*a, v2*v3*a, 2.0*v**2.0*s+v3**2.0*a])
	dqdo = (1.0/(2.0*v**3))*dqdo

	M1 = drot_dq.reshape((4,9)).T
	M2 = dot(M1,dqdo).T
	drot = reshape(M2,(3,3,3))

	return drot

def EigSpaceScore(x,Xa,Xb,Vab,sig):
	o = x[0:3]
	p = x[3:]
	Xar = dot(Rot(o),Xa) + reshape(p,(3,1)) 
	Na,Nb = Vab.shape
	dX = zeros((3,Na,Nb))
	D2 = zeros((Na,Nb))
	for i in arange(3):
		dX[i,:,:] = tile(Xar[i,:],(Nb,1)).T - tile(Xb[i,:],(Na,1))
		D2 = D2 + dX[i,:,:]*dX[i,:,:]
	Cab = exp(-0.5*(D2/sig**2))
	CV = Cab*Vab
	scr = sum(CV)
	return -scr

def EigSpaceGrad(x,Xa,Xb,Vab,sig):
	o = x[0:3]
	p = x[3:]
	Xar = dot(Rot(o),Xa) + reshape(p,(3,1)) 
	Na,Nb = Vab.shape
	dX = zeros((3,Na,Nb))
	D2 = zeros((Na,Nb))
	for i in arange(3):
		dX[i,:,:] = tile(Xar[i,:],(Nb,1)).T - tile(Xb[i,:],(Na,1))
		D2 = D2 + dX[i,:,:]*dX[i,:,:]
	Cab = exp(-0.5*(D2/sig**2))
	CV = Cab*Vab
	dfdX = zeros(Xar.shape)
	for i in arange(3):
		dfdX[i,:] = -(1.0/sig**2)*sum(dX[i,:,:]*CV,axis=1).T
	drot = DRot(o)
	grad = zeros(6)
	for i in arange(3):
	   grad[i] = sum(dot(drot[i,:,:],Xa)*dfdX)
	grad[3:6] = sum(dfdX,axis=1)
	return -grad

def PrincipleAxesTransform(X):
	# center structure
	X_centroid = reshape(X.mean(axis=1),(3,1))
	Y = X - X_centroid 
	# rotate structure
	C = dot(Y,Y.T)
	(lam,V) = Eig(C)
	return dot(V.T,Y)

def AlignEigenspaces(Xa,Xb,sig,tol,itermax):
	#Xa = PrincipleAxesTransform(Xa)
	#Xb = PrincipleAxesTransform(Xb)
	Ca = ContactMatrix(Xa,Xa,sig)
	Cb = ContactMatrix(Xb,Xb,sig)
	Cab = ContactMatrix(Xa,Xb,sig)
	[la,Va] = Eig(Ca)
	[lb,Vb] = Eig(Cb)
	Ia,Ib,maxScore = AlignEigenvalues(la,lb,0.01)
	x = array([0.001,0,0,0,0,0])
	iter = 0
	scr = 0 
	dscr = tol
	while (dscr >= tol) and (iter <= itermax):
		Vab = dot(Va[:,Ia],Vb[:,Ib].T)
		results = fmin_bfgs(EigSpaceScore,x,EigSpaceGrad,args=(Xa,Xb,Vab,sig),full_output=1,disp=0)
		#results = fmin_ncg(EigSpaceScore,x,EigSpaceGrad,args=(Xa,Xb,Vab,sig),full_output=1,disp=0)
		xopt  = results[0]
		scr1  = results[1]
		Xaopt = dot(Rot(xopt[0:3]),Xa) + reshape(xopt[3:],(3,1))
		Cab = ContactMatrix(Xaopt,Xb,sig)
		Score = dot(Va.T,dot(Cab,Vb))
		Ia,Ib,scr2 = DP(Score,0)
		dscr = abs(scr - scr2)
		x = xopt
		scr = scr2
		iter = iter + 1
		print
		print 'AlignEigenspaces --> iter = %d' % iter
		print 'AlignEigenspaces --> superposition score = %f' % -scr1
		print 'AlignEigenspaces --> DP score = %f' % scr2
		print 'AlignEigenspaces --> number modes aligned = %d' % len(Ia)
		print
	return Ia,Ib,Xaopt,xopt,scr

