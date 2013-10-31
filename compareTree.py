from kdtree import KDTree
import sys, time, random, csv, math
def randPoints(amt):
	data = []
	for i in range(amt):
		x = random.randint(-1000,1000)
		y = random.randint(-1000,1000)
		z = random.randint(-1000,1000)
		#c = random.randint(-100000,100000)
		#d = random.randint(-100000,100000)
		#e = random.randint(-100000,100000)
		#f = random.randint(-100000,100000)
		#r = random.randint(-100000,100000)
		#g = random.randint(-100000,100000)
		data.append((x, y, z))
	return data
def square_distance(pointA, pointB):
    # squared euclidean distance
    distance = 0
    dimensions = len(pointA) # assumes both points have the same dimensions
    for dimension in range(dimensions):
        distance += (pointA[dimension] - pointB[dimension])**2
    return distance
def average(list):
	tot = 0
	for i in list:
		tot +=i
	return tot/len(list)

def bruteForceNN(pointsList, point):
	minimum = square_distance(point, pointsList[0])
	point1 = pointsList[0]
	for i in range(len(pointsList)):
		d = square_distance(point, pointsList[i])
		if(d<minimum):
			minimum = d
			point1 = pointsList[i]
	return(minimum, point1)

def bruteForceRange(pointsList, point, range):
	results=[]
	for p in pointsList:
		d = square_distance(p, point)
		if(d<=range):
			results.append(p)
	return results

def readInputPoints(filename):
	floats =[]
	with open(filename, "rb") as f:
		r = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
		next(r)
		for row in r:
			floats.append([item for number, item in enumerate(row) if item])
	return floats

point = (6,10,0)
#data = readInputPoints("contactmatrix.csv")
kdTreeConstruct = []
kdTreeSearch = []
bruteForceSearch =[]
for i in range(25):
	data = randPoints(10000)
	startCon = time.time()
	tree = KDTree.construct_from_data(data)
	constructionTime = time.time() - startCon
	start = time.time()
	#nearest = tree.querynn(query_point=point, t=3)
	nearest = tree.queryrange(query_point = point, r = 250)
	searchTime = time.time() -start
	kdTreeConstruct.append(constructionTime)
	kdTreeSearch.append(searchTime)
	st = time.time()
	t = bruteForceRange(data, point, 250*250)
	end = time.time()-st
	bruteForceSearch.append(end)
print(average(kdTreeConstruct))
print(average(kdTreeSearch))
print(average(bruteForceSearch))