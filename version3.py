from kdtree import KDTree
import sys, time, random, csv, math
from align_eigenspaces import *
from ContactGeometry1 import *
import numpy as np
import matplotlib as plt
from scipy.optimize import fmin_bfgs


def dot_product(Va, Vb):
    #No need to do transpose here. its done earlier
    d=0
    for i in range(len(Va)):
        d+=Va[i]*Vb[i]
    return d
	
def newContactScoreFxn(x, Xa, Xb, Va, Vb, sig):
	#(n1, n2, n3, px, py, pz)
	o = x[0:3]
	p = x[3:]
	Xar = np.dot(Rot(o), Xa)+ np.reshape(p, (3,1))
	Cab = ContactMatrix(Xar,Xb,sig)
	scr = -(np.dot(Va.T,np.dot(Cab,Vb))).sum()
	return scr
FileName = 'PDBFiles/d1cida1.ent'
pdb = ReadPDB(FileName)
FileName = 'PDBFiles/d2rhea_.ent'
pdb2 = ReadPDB(FileName)

rawProteinA = CA_coordinates(pdb2)
rawProteinB = CA_coordinates(pdb2)

sig = 8 
tol = 0.001 
itermax = 100


#full list of atoms in protein a
Xafull = np.array(pdb[['x','y','z']],dtype='float64')
Xa1 = Xafull[:15]#only take first 15 atoms (for now)


#full list of atoms in protein b
Xbfull = np.array(pdb2[['x','y','z']],dtype='float64')
Xb1 =Xbfull[:15]#only take first 15 atoms (for now)

Xa = PrincipleAxesTransform(rawProteinA)#cut all but 15 atoms
Xb = PrincipleAxesTransform(rawProteinB)
Ca = ContactMatrix(Xa,Xa,sig)
Cb = ContactMatrix(Xb,Xb,sig)
tempa= Eig(Ca)
la = tempa[0]
Va = tempa[1]
#####################
tempb = Eig(Cb)
lb = tempb[0]
Vb = tempb[1]

#--------------------------------------------------
x0 = np.array([0.001,0,0,0,0,0])
xopt,scr,A,B,a,b,c = fmin_bfgs(ContactScore,x0,args=(Xa,Xb,Va,Vb,sig),full_output=1,disp=1)
#--------------------------------------------------

CA_Coord_ProteinA = rawProteinA.T
CA_Coord_ProteinB = rawProteinB.T
limitedListCA_ProtA = CA_Coord_ProteinA
limitedListCA_ProtB = CA_Coord_ProteinB

#listXa is the list of atoms in protein a
listXa=[]
for i in range(len(limitedListCA_ProtA)):
	tupp = tuple([round(limitedListCA_ProtA[i][0], 3)])+tuple([round(limitedListCA_ProtA[i][1],3)])+tuple([round(limitedListCA_ProtA[i][2], 3)])+tuple([Va[i].T])
	#eigenvector Va.T is appended to each node
	listXa.append(tupp)

#listXb is the list of atoms in protein b
listXb=[]
for i in range(len(limitedListCA_ProtB)):
	tupp = tuple([round(limitedListCA_ProtB[i][0], 3)])+tuple([round(limitedListCA_ProtB[i][1],3)])+tuple([round(limitedListCA_ProtB[i][2], 3)])+tuple([Vb[i]])
	#eigenvector Vb is appended to each node
	listXb.append(tupp)

data1 = listXa
data2 = listXb
Tree1 = KDTree.construct_from_data(data1)
Tree2 = KDTree.construct_from_data(data2)
score = 0
#print("####################################")
#Times for KD Tree approach
startT = time.time()
for i in range(len(data1)):
	#finds the atoms within radius 30 of query pt
	score += Tree2.queryrange(query_point=data1[i], r = 50)
solveTime = time.time() - startT
print(solveTime)



#Time for non-tree approach
#startT = time.time()
#total =0
#for i in range(len(data1)):
#	#all points within dist r of the nn points on Tree2!
#	rangepoints =Tree2.queryrange(query_point = data1[i], r = 50)
#	#print(rangepoints)
#	for k in range(len(rangepoints)):
#		sd = rangepoints[k][-1] #Cij value of contact matrix
#		C = math.exp(-0.5*(sd/64))
#		#No Longer Necessary, calc happen recursively
#		#only prints non-zero values of C
#		if(C>.01): 
#			#print(C) 
#			total = total+C*dot_product(data1[i][-1], rangepoints[k][0][-1])		
#solveTime = time.time() - startT
#print(solveTime)