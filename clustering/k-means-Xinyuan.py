#!/usr/bin/env python3
"""
Author: Xinyuan Min
Student number: 950829573070
Implementation of the k-means clustering algorithm
"""
# import statements
import matplotlib.pyplot as plt
import numpy as np


def csv_parser(lines):
    """Return numpy.matrix of points coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]

    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: matrix with each row as data points, colunmns as coordinates of the
    point(without label).
    """
    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(list(map(float, items[1:]))) #first item is the label
        except ValueError: #must be the header
            continue
    return np.mat(data_points) # return a matrix

def euc_distance(point1, point2):
    """Return a numpy.float which is the Euclidian distance between point1 and point2

    point1, point2: numpy array with coordinates as [x, y, z...] of a single point

    this function is to calculate the Euclidian distance between point1 and point2
    """
    distance = np.sqrt(np.sum(np.power(point1-point2, 2)))
    return distance

def init_centroids(dataset, k):
    """Return a numpy matrix of centorids coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]

    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    k: a positive int, the number of clusters, normally k is larger than 1

    this function is to randomly select k points as initial centroids.each row of return
    contains the information of a centroid， the columns store the coordinates
    on different dimensions of a centroid as [x, y, z...]
    """
    datanr, dim = dataset.shape
    # initiate centroids as a matrix in shape (k, dimension of data)
    centroids = np.zeros((k, dim))
    for i in range(k):
        rand_int = np.random.randint(0,datanr)
        # randomly select k datapoints as initial centroid
        centroids[i,:] = dataset[rand_int,:]
    return centroids


def assign_centroid(dataset, centroids, k):
    """Return an array of assignment result, element in array is a int between 0 to k-1

    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    centroids：a numpy matrix of a set of centroids， each row is the information of
    a centroid， the columns store the coordinates on different dimensionsof a
    centroid as [x, y, z...] 
    k: a positive int, the number of clusters, normally k is larger than 1

    given the centroids, recorde the cluster being assigned of each data point as
    an int ranges between 0 to (k-1)
    """
    datanr, dim = dataset.shape
    # store the distances and belonging luster for each datapoint
    recorder = np.zeros((datanr, 1))
    for i in range(datanr): 
        data = dataset[i,:] # for each data point
        min_distance = np.inf
        cluster_index = None
        for j in range(k): # for each centriod
            centroid = centroids[j,:]
            distance = euc_distance(data, centroid)
            if distance < min_distance:
                # update min_distance & cluster index if closer distance appears
                min_distance = distance
                cluster_index = j
        recorder[i,:] =cluster_index
    return recorder

        
def calculate_centroids(dataset, centroids, k):
    """Return a matrix of coordinates of centorids as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    centroids：a numpy matrix of a set of centroids， each row is the information of
    a centroid， the columns store the coordinates on different dimensionsof a
    centroid as [x, y, z...] 
    k: a positive int, the number of clusters, normally k is larger than 1

    Given current centorids, assgin points to each clusters and recalculate centroids,
    the coordinates of new centroids are the mean coordinates of points in the cluster 
    """
    datanr, dim = dataset.shape
    new_centroids = np.zeros((k,dim))
    empty = False
    recorder = assign_centroid(dataset, centroids, k)             
    for j in range(k):               
        # get indices of all data points in cluster k
        indices = np.nonzero(recorder[:, 0] == j)[0]
        if len(indices) == 0: # if cluster k is empty
            # replace the centroid with random datapoints 
            centroids[j] = dataset[np.random.randint(0,datanr),:]
            return calculate_centroids(dataset, centroids, k)
        else:
            points_in_cluster = dataset[indices]
            new_centroids[j,:] = np.mean(points_in_cluster[:], axis = 0)
    return new_centroids


def k_means(dataset, k):
    """Return matrix of final centroids coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]],
    an array of assignment result,for which element is an int between 0 and k-1
    a positve int indicates number of iterations, 
    a list of int to reflect the size of each cluster
    
    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    k: a positive int, the number of clusters, normally k is larger than 1

    implement k-means on dataset for specified k. Calculate centroids and assign
    clusers iteratively utill the location of centroids converge
    """
    centroids = init_centroids(dataset, k)
    centriod_changed = True
    iteration = 0
    while centriod_changed:
        updated_centroids = calculate_centroids(dataset, centroids, k)
        iteration += 1
        if (updated_centroids != centroids).any():
            # if centroids changed, go back to while loop and update centroids
            centriod_changed = True
            centroids = updated_centroids # update centroids
        else:
            centriod_changed = False
    assign_result = assign_centroid(dataset, centroids, k)
    cluster_size = []
    for j in range(k):
        indices = np.nonzero(assign_result[:, 0] == j)[0]
        cluster_size.append(len(indices))
    return centroids, assign_result, iteration, cluster_size

        
class DimensionError(Exception): 
    """self-defined exception for wrong dimensions of input data
    """
    def __init__(self,ErrorInfo):
        self.errorinfo = ErrorInfo
 
    
def visulization(dataset, k):
    """plot a scatter to present points of different cluster in different colors

    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    k: a positive int, the number of clusters, normally k is larger than 1
    """
    datanr, dim = dataset.shape
    if dim != 2: 
        raise DimensionError('dataset must be of 2 dimension!')
    assign_result = k_means(dataset, k)[1]
    for i in range(datanr): # for each data points
         # get the cluster index the points belongs to
        cluster_index = assign_result[i,0]
        # give different color for points from different clusters 
        plt.scatter(dataset[i,0], dataset[i,1], c='C'+str(int(cluster_index)))
    plt.show()
        

def wgss_bgss(dataset, k):
    """Return a 2-decimal percent to represent the WGSS-BGSS ration of cluster result

    dataset: a numpy matrix of coordinates of serval data points, each point is
    represented as [x, y, z...]
    k: a positive int, the number of clusters, normally k is larger than 1

    WGSS: within-cluster sum of squares. BGSS: between-cluster sum of squares
    WGSS-BGSS ratio = WGSS/BGSS, lower ration indicates "better clusers"
    """
    wgss = 0
    k_result = k_means(dataset, k)
    assign_result = k_result[1]
    centroids = k_result[0]
    for j in range(k):
        # get indices of all data points in cluster k
        indices = np.nonzero(assign_result[:, 0] == j)[0]
        cluster_nr = len(indices) # number of points in this cluster
        points_in_cluster = dataset[indices]
        a = np.sum(np.power(points_in_cluster-centroids[j], 2))/cluster_nr
        wgss += a
        
    bgss = 0
    # calculate the centroid of the whole dataset
    overall_centroid = np.zeros((1,dataset.shape[1]))
    overall_centroid[0,:] = np.mean(dataset,axis = 0)
    for j in range(k):
        indices = np.nonzero(assign_result[:, 0] == j)[0]
        cluster_nr = len(indices)
        bgss += cluster_nr*euc_distance(centroids[j], overall_centroid)**2
    return  '{:.2%}'.format(wgss/bgss)  

    
       
    
    
    

if __name__ == "__main__":
    file1 = open('2dtest.csv')
    lines1 = file1.readlines()
    dataSet = csv_parser(lines1)

    file_nd = open('LargeSet_1.csv')
    lines_nd = file_nd.readlines()
    dataSet_nd = csv_parser(lines_nd)

    file_big = open('LargeSet_2.csv')
    lines_big = file_big.readlines()
    dataSet_big = csv_parser(lines_big)


    # Question 1
    for i in range(5):
        print('iteration needed for converge:',k_means(dataSet, 3)[2],'\n')
  
    # Question 2
    for k in range(2,7):
        for i in range(5):
            print('k:',k,\
                  'iteration needed for converge:',k_means(dataSet_nd, k)[2],\
                  'cluster size:',k_means(dataSet_nd, k)[3])

    # Question 3
    visulization(dataSet_big, 2)

    # Question 4
    for k in range(2,7):
        for i in range(10):
            print('k:',k,'WGSS/BGSS ratio:',wgss_bgss(dataSet_nd, k))
        print('\n')
    print(k_means(dataSet_nd, 5)[1][:,0])

 


  
          
    
