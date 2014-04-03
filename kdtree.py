#!/usr/bin/python
# encoding: utf-8

import math
__version__ = "1r11.1.2013"
__all__ = ["KDTree"]

def square_distance(pointA, pointB):
    # squared euclidean distance
    distance = 0
    dimensions = len(pointA)-1 # assumes both points have the same dimensions
    for dimension in range(dimensions):
        distance += (pointA[dimension] - pointB[dimension])**2
    return distance

def dot_product(Va, Vb):
    #No need to do transpose here. its done earlier
    d=0
    for i in range(len(Va)):
        d+=Va[i]*Vb[i]
    return d
class KDTreeNode():
    def __init__(self, point, left, right):
        self.point = point
        self.left = left
        self.right = right
    
    def is_leaf(self):
        return (self.left == None and self.right == None)

class KDTreeRange():
    """Internal structure used in range search"""
    def __init__(self, query_point, r):
        self.query_point = query_point
        self.range = r
        #self.found_points = [] #Toggle for full C matrix calculation
        self.ContactScore = 0 #contact score between query point and other points
    def add(self, point):
        sd = square_distance(point, self.query_point) #square distance between 2 atoms
        if (sd<=self.range):
            #Contact Score calculation happens here
            Cij = math.exp(-0.5*(sd/64))
            self.ContactScore += Cij*dot_product(self.query_point[-1], point[-1])
            #self.found_points.append(((point),sd)) #Toggle for full C matrix calculation
        return
    def getRange(self):
        #return self.found_points #Toggle for full C matrix calculation
        return self.ContactScore #Toggle for recursive calculation


class KDTreeNeighbours():
    """ Internal structure used in nearest-neighbours search.
    """
    def __init__(self, query_point, t):
        self.query_point = query_point
        self.t = t # neighbours wanted
        self.largest_distance = 0 # squared
        self.current_best = [] #structure [[(point), distance], [(point), distance]...]

    #The largest radius to search within. 
    #Until t points are found, it is infinite*
    def calculate_largest(self):
        if self.t >= len(self.current_best):
            #infinite *
            self.largest_distance = self.current_best[-1][1]
        else:
            #largest radius  = current_best's distance to the query_point
            self.largest_distance = self.current_best[self.t-1][1]

    def add(self, point):
        sd = square_distance(point, self.query_point)

        # run through current_best, try to find appropriate place
        #enumerate documentation: http://docs.python.org/2.3/whatsnew/section-enumerate.html
        for i, e in enumerate(self.current_best):
            if i == self.t:
                return # enough neighbours, this one is farther, let's forget it
            if e[1] > sd:
                self.current_best.insert(i, [point, sd])
                self.calculate_largest()
                return
        # append it to the end otherwise
        self.current_best.append([point, sd])
        self.calculate_largest()

    #returns the t best points
    def get_best(self):
        return [element[0] for element in self.current_best[:self.t]]
        
class KDTree():
    """ KDTree implementation.
    
        Example usage:
        
            from kdtree import KDTree
            
            data = <load data> # iterable of points (which are also iterable, same length)
            point = <the point of which neighbours we're looking for>
            
            tree = KDTree.construct_from_data(data)
            nearest = tree.query(point, t=4) # find nearest 4 points
    """
    
    def __init__(self, data):
        def build_kdtree(point_list, depth):
            # code based on wikipedia article: http://en.wikipedia.org/wiki/Kd-tree
            if not point_list:
                return None

            #As one moves down the tree, 
            #one cycles through the axes used to select the splitting planes. 
            #http://en.wikipedia.org/wiki/K-d_tree#Construction
            axis = depth % (len(point_list[0])-1) # assumes all points have the same dimension

            #Sort points list based on selected axis,
            #and choose median as pivot point.
            #Currently this is nlog(n) avg time
            #TODO: better selection method, linear-time selection, distribution
            #TODO: use quickselect to find the K^th value in O(n) time
            point_list.sort(key=lambda point: point[axis])
            median = len(point_list)/2 # choose median

            # create node and recursively construct subtrees
            node = KDTreeNode(point=point_list[median],
                              left=build_kdtree(point_list[0:median], depth+1),
                              right=build_kdtree(point_list[median+1:], depth+1))
            return node
        
        self.root_node = build_kdtree(data, depth=0)
    
    #call this method to initialize construction of Tree
    @staticmethod
    def construct_from_data(data):
        tree = KDTree(data)
        return tree
    def queryrange(self, query_point, r):
        statistics = {'nodes_visited': 0, 'far_search': 0, 'leafs_reached': 0}
        def r_search(node, query_point, r, depth, points_in_range):
            if node == None:
                return
            if node.is_leaf():
                points_in_range.add(node.point)
                return
            axis = depth % (len(query_point)-1)
            near_subtree = None
            far_subtree = None

            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left

            r_search(near_subtree, query_point, r, depth+1, points_in_range)
            points_in_range.add(node.point)
            if (node.point[axis] - query_point[axis])**2 < points_in_range.range:
                r_search(far_subtree, query_point, r, depth+1, points_in_range)

            return

        # if there's no tree, there's no range
        if self.root_node != None:
            r = r*r #we used distance sqd to make it optimized
            points = KDTreeRange(query_point, r)
            r_search(self.root_node, query_point, r, depth=0, points_in_range=points)
            result = points.getRange()
        else:
            result = []
        return result


    #nearest neighbor search
    def querynn(self, query_point, t=1):
        statistics = {'nodes_visited': 0, 'far_search': 0, 'leafs_reached': 0}

        
        #can be done using a priority Q: http://web.engr.oregonstate.edu/~tgd/classes/534/slides/part3.pdf
        #future reference
        def nn_search(node, query_point, t, depth, best_neighbours):
            if node == None:
                return
            
            #statistics['nodes_visited'] += 1
            
            # if we have reached a leaf, let's add to current best neighbours,
            # (if it's better than the worst one or if there is not enough neighbours)
            if node.is_leaf():
                #statistics['leafs_reached'] += 1
                best_neighbours.add(node.point)
                return
            
            # this node is no leaft
            
            # select dimension for comparison (based on current depth)
            axis = depth % (len(query_point)-1)
            
            # figure out which subtree to search
            near_subtree = None # near subtree
            far_subtree = None # far subtree (perhaps we'll have to traverse it as well)
            
            # compare query_point and point of current node in selected dimension
            # and figure out which subtree is farther than the other
            if query_point[axis] < node.point[axis]:
                near_subtree = node.left
                far_subtree = node.right
            else:
                near_subtree = node.right
                far_subtree = node.left

            # recursively search through the tree until a leaf is found
            nn_search(near_subtree, query_point, t, depth+1, best_neighbours)

            # while unwinding the recursion, check if the current node
            # is closer to query point than the current best,
            # also, until t points have been found, search radius is infinity
            best_neighbours.add(node.point)
            
            #check whether there could be any points on the other side of the
            #splitting plane that are closer to the query point than the current best
            #True if the radius of the best_distance crosses over the splitting plane.
            if (node.point[axis] - query_point[axis])**2 < best_neighbours.largest_distance:
                nn_search(far_subtree, query_point, t, depth+1, best_neighbours)
            
            return
        
        # if there's no tree, there's no neighbors
        if self.root_node != None:
            neighbours = KDTreeNeighbours(query_point, t)
            nn_search(self.root_node, query_point, t, depth=0, best_neighbours=neighbours)
            result = neighbours.get_best()
        else:
            result = []
        
        #print statistics
        return result
