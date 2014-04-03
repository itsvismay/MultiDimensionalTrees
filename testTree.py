from kdtree import KDTree
import sys, time, random, csv, math
from align_eigenspaces import *
from ContactGeometry1 import *
import numpy as np
import matplotlib as plt

#Not Used anymore
# def readInputPoints(filename):
# 	floats =[]
# 	with open(filename, "rb") as f:
# 		r = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
# 		next(r)
# 		for row in r:
# 			floats.append([item for number, item in enumerate(row) if item])
# 	return floats

sig = 8 
tol = 0.001 
itermax = 100
FileName = 'PDBFiles/d1cida1.ent'
pdb = ReadPDB(FileName)
Xa1 = np.array(pdb[['x','y','z']],dtype='float64')

FileName = 'PDBFiles/d2rhea_.ent'
pdb2 = ReadPDB(FileName)
Xb1 = np.array(pdb[['x','y','z']],dtype='float64')
Xa = PrincipleAxesTransform(Xa1[:15].T)
Xb = PrincipleAxesTransform(Xb1[:15].T)
Ca = ContactMatrix(Xa,Xa,sig)
Cb = ContactMatrix(Xb,Xb,sig)
tempa= Eig(Ca)
la = tempa[0]
Va = tempa[1]
#####################
tempb = Eig(Cb)
lb = tempb[0]
Vb = tempb[1]
print(Va)
#--------------------------------------------------
x0 = np.array([0.001,0,0,0,0,0])
xopt,scr,A,B,a,b,c = fmin_bfgs(ContactScore,x0,args=(Xa,Xb,Va,Vb,sig),full_output=1,disp=1)
#--------------------------------------------------

#pdb1 = ReadPDB('PDBFiles/1QKU.pdb')
#X = np.array(pdb1[['x','y','z']],dtype='float64')
#X1 = X[:len(X)]
#print(X1)
listXa=[]
for i in range(len(Xa1)):
	tupp = tuple([round(X1[i][0], 3)])+tuple([round(X1[i][1],3)])+tuple([round(X1[i][2], 3)])+tuple([i])
	listXa.append(tupp)

listXb=[]
for i in range(len(Xb1)):
	tupp = tuple([round(X1[i][0], 3)])+tuple([round(X1[i][1],3)])+tuple([round(X1[i][2], 3)])+tuple([i])
	listXb.append(tupp)

data1 = listXa
data2 = listXb
Tree1 = KDTree.construct_from_data(data1)
Tree2 = KDTree.construct_from_data(data2)
"""
querypoint = (42.5150,16.8110,16.5530, 1) 
nnpoints = Tree1.querynn(querypoint, t=4) #t nearest neighbors of query point
#print(nnpoints)
#total = 0

for i in range(len(nnpoints)):
	#all points within dist r of the nn points on Tree2!
	rangepoints =Tree2.queryrange(query_point = nnpoints[i], r = 4)
	#print(rangepoints)
	for k in range(len(rangepoints)):
		C = rangepoints[k][-1] #Cij value of contact matrix

		#No Longer Necessary, calc happen recursively
		#only prints non-zero values of C
		#if(C!=0): 
			#print(C) 
			#atom1
			#print(nnpoints[i][-1])
			#total = total+V*Distance*V
			#atom2
			#print(rangepoints[k][0][-1])
		
		# Vrangepoints[k][-1])
	#print(rangepoints)
"""