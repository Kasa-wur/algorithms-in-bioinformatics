#!/usr/bin/env python
"""
Author:
Student number:
Implementation of the k-means clustering algorithm

Hints:
- write a function to obtain Euclidean distance between two points.
- write a function to initialize centroids by randomly selecting points 
    from the initial set of points. You can use the random.sample() method
- write a function to find the closest centroid to each point and assign 
    the points to the clusters.     
- write a function to calculate the centroids given the clusters
- write a function to implement k-means
- write a function to calculate WGSS given the clusters and the points
- write a function to calculate BGSS given the clusters and the points
"""

# import statements
import random   


def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(map(float, items[1:]))#first item is the label
        except ValueError: #must be the header
            continue
    return data_points


if __name__ == "__main__":
    file = open('2dtest.csv')
    lines1 = file.readlines()
    print(csv_parser(lines1))
    

    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.
    
   
