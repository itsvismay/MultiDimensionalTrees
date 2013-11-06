from kdtree import KDTree
import sys, time, random, csv, math

from ContactGeometry1 import *
import numpy as np
import matplotlib as plt

def readInputPoints(filename):
	floats =[]
	with open(filename, "rb") as f:
		r = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
		next(r)
		for row in r:
			floats.append([item for number, item in enumerate(row) if item])
	return floats

FileName = 'PDBFiles/1QKU.pdb'

pdb = ReadPDB(FileName)
X = np.array(pdb[['x','y','z']],dtype='float64')
X1 = X[:10]
listX=[]
for i in range(len(X1)):
	print(i)
	tupp = tuple([round(X1[i][0], 3)])+tuple([round(X1[i][1],3)])+tuple([round(X1[i][2], 3)])+tuple([i])
	listX.append(tupp)
data1 = listX
data2 = listX
Tree1 = KDTree.construct_from_data(data1)
Tree2 = KDTree.construct_from_data(data2)

querypoint = (42.5150,16.8110,16.5530, 1)
nnpoints = Tree1.querynn(querypoint, t=4)
#print(nnpoints)
total = 0
for i in range(len(nnpoints)):
	rangepoints =Tree2.queryrange(query_point = nnpoints[i], r = 4)
	#print(rangepoints)
	for k in range(len(rangepoints)):
		C = rangepoints[k][-1]
		if(C!=0):
			print(C)
			#atom1
			#print(nnpoints[i][-1])
			#total = total+V*Distance*V
			#atom2
			#print(rangepoints[k][0][-1])
		
		# Vrangepoints[k][-1])
	#print(rangepoints)