#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Xinyuan Min
Student nr: 950829573070
Script to impelement Eulerian cycle/ path algorithm
"""

################ determine Eulerian / has Eulerian path############

def all_vertices(graph):
    """Returns a set of string/ int，each string/ integer
    corresponds to a vertex in the graph

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    to get every distinct vertex in a graph
    """
    vertices = set()
    for key, value in graph.items():
        vertices.add(key)
        for node in value:
            vertices.add(node)
    return vertices

def count_in_out(graph):
    """Return a dictionary with string/ int as key, each
    key corresponds to a vertex. a list containing 2 integer
    as values. The first int is the indegree and second
    is the outdegree of the vertex

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    to count the indegree and outdegree of each vertex in graph
    """
    count_in_out = {}
    # get all vertices in the graph 
    vertices = all_vertices(graph)
    for vertex in vertices:
        count_in_out[vertex] = [0, 0]
    for node in graph:
        # the length of connected nodes is outdegree
        outdegree = len(graph[node])
        count_in_out[node][1] = outdegree
        for connected_node in graph[node]:
            # the number of a node been taken as values is indegree
            count_in_out[connected_node][0] += 1
    return count_in_out
        
    
def is_eulerian(graph):
    """Return True if the graph is Eulerian, False otherwise

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    a graph is Eulerian if indegree = outdegree holds for
    all vertices
    """
    count = count_in_out(graph)
    eulerian = True
    for vertex in count:
        indegree = count[vertex][0]
        outdegree = count[vertex][1]
        if indegree != outdegree:
            eulerian = False
            break
    return eulerian


def has_eulerian_path(graph):
    """Return True if graph has Eulerian path, False otherwise

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    a connected graph has a Eulerian path of and only if it
    contains at most 2 semi-balanced vertices and all other
    vertices are balanced
    """
    count = count_in_out(graph)
    semi_balanced = 0
    eulerian_path = False
    for vertex in count:
        indegree = count[vertex][0]
        outdegree = count[vertex][1]
        if abs(indegree - outdegree) == 1:
            semi_balanced += 1
    if semi_balanced <= 2:
        eulerian_path = True
    return eulerian_path


        
######################## find Eulerian cycle ##########################

def delete_dictionary(graph, node):
    """Return dictionary similar to graph except for 1 value been
    deleted
    
    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices
    node：a single-character string or int，the first value
    of node will be deleted

    to delete the fist value（connected node）of node in graph，
    because the used edge cannot be taken again
    """
    if len(graph[node]) == 1:
        # clear node and its value in graph
        graph.pop(node)
    else:
        # delete the first value of node
        next_node = graph[node][0]
        graph[node].remove(next_node)
    return graph

def find_subcycle(graph, node, start_node, used_edges):
    """First return a list of tuples to denote the sub-cycle，
    each tuple represents one edge. e.g.（'A','B') means edge AB
    the second return value is a dictionary represent the
    remaining graph excluding the sub-cycle

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices
    node：string or int， the node currently traverse through
    start_node: string or int, the start node of the sub-cycle
    used_edge:list, containing used edges 
    
    Find a sub-cycle given a start node, as well as store the
    traversed edges and return the remaining graph
    """
    # choose the first connectd node in value list
    next_node = graph[node][0]
    choosen_edge = (node, next_node)
    # append the choosen_edges of each recursion
    used_edges.append(choosen_edge)
    # get graph after the chosen edge been deleted
    graph = delete_dictionary(graph, node)
    # keep on traversing until going back to start node 
    if next_node != start_node:
        # recursion，next time start traversing from nextnode
        return find_subcycle(graph, next_node, start_node, used_edges)
    else:
        return used_edges, graph


def find_all_cycle(graph, paths, used_edges):
    """Return a list of list, each inner list stores the path
    of a sub-cycle

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices
    paths: an empty list, waits for recording sub-paths
    used_edges: an empty list, waits for storing edges of sub-cycle

    This function recursively calls find_subcycle function to
    get the overall paths for every subcycles in graph
    """
    if graph != {}:
        # take the key in the first position of the dictionary
        # as start vertex
        start_node = list(graph.keys())[0]
        # initiate used_edges list
        used_edges = []
        search_sub_cycle = find_subcycle(graph, start_node,
                                         start_node, used_edges)
        # append the path of a sub-cycle in paths list 
        paths.append(search_sub_cycle[0])
        # remaining graph after one recursion 
        graph = search_sub_cycle[1]
        return find_all_cycle(graph, paths, used_edges)
    else:
        return paths

def find_eulerian_cycle(graph):
    """Return a list of list, each inner list stores the path
    of a sub-cycle

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    this function is used to pass the two empty list as parameters
    to the find_all_cycle function
    """
    paths = []
    used_edges = []
    path = find_all_cycle(graph, paths, used_edges)
    return path

def integrate_2_paths(path1, path2, cross_node):
    """Return a list of strings, the list is the integrated
    path1 and path2

    path1：list of string， first sub-path in paths list
    path2：list of string， second sub-path in paths list
    cross_point: single- character string/ integer
    
    to integrate path1 and path2 into one single path
    """
    for edge in path1:
        # if the connected vertex of a edge is the cross node
        if edge[-1] == cross_node:
            insert_index = path1.index(edge) + 1
    combined_path = path1[:insert_index] + path2 + path1[insert_index:]
    return combined_path
        
def combine(graph,paths):
    """Return list of tuples to represent the Eulerian cycle of graph,
    each tuple is a edge, traverse edges in order forms a Eulerian cycle

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices
    paths: list of list, return value of find_eulerian_cycle(graph)

    this function integrate the first and second sub-paths recursively,
    the final integrated list is the Eulerian cycle of graph
    """
    if is_eulerian(graph) == True:
        if len(paths) >= 2:
            path1 = paths[0]
            path2 = paths[1]
            # corss_node of fisrt and second sub-path is the first
            # character of first edge in the second sub-path 
            cross_node = path2[0][0]        
            combined_path = integrate_2_paths(path1, path2, cross_node)
            paths.pop(0)
            paths.pop(0)
            paths.insert(0,combined_path)
            return combine(graph, paths)
        else:
            return paths[0]
    # if input graph is not Eulerian, raise an error     
    else:
        raise InputError("invalid eulerian graph")


    
######################### find Eulerian path ##########################

def find_start_end(graph):
    """Return two string， the first return value denotes the start
    vertex in graph， the second end vertex
    
    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices

    to get the start and vertex of a graph
    """
    count = count_in_out(graph)
    for vertex, in_out in count.items():
        if in_out[1] - in_out[0] == 1:
            start = vertex
        if in_out[0] - in_out[1] == 1:
            end = vertex
    return start, end

def find_eulerian_path(graph):
    """Return a list of tuples，the list represents the eulerian
    path of graph， each tuple means an edge

    graph：dictionary, string/ int as key to denote the vertex，
    list as value to represent the conncect vertices
    """
    start = find_start_end(graph)[0]
    end = find_start_end(graph)[1]
    # add an extra edge between start and end vertices
    graph[end] = [start]
    # find eulerian cycle of the ajusted graph 
    cycle = find_eulerian_cycle(graph)
    combined_cycle = combine(graph,cycle)
    # get the index of added edge 
    index = combined_cycle.index((end, start))
    # break and recombine the cycle into eulerian path
    path = combined_cycle[index+1:] + combined_cycle[:index]
    return path 
    
    
###################### create graph from spectrum #####################
    
def create_graph(spectrum):
    """Return a dictionary, represents the graph for the given spectrum

    spectrum: list of strings, each string is a kmer for DNA sequence

    this function creates a graph for a given spectrum
    """
    two_mer_set = set()
    for three_mer in spectrum:
        # add distinct 2-mer to two_mer_set
        two_mer_set.add(three_mer[:2])
        two_mer_set.add(three_mer[-2:])
    graph = {}
    for two_mer in two_mer_set:
        graph[two_mer] = []
        for three_mer in spectrum:
            if three_mer[:2] == two_mer:
                graph[two_mer].append(three_mer[-2:])         
    return graph
        
                  
        

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}

    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {1:[2], 2:[3], 3:[7,1],\
        4:[5,10],5:[6],6:[7],7:[8,9,11],\
        8:[4],9:[3,4],\
        10:[9],11:[12],12:[7]}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']


   # Question 1
    if is_eulerian(graph_822) == True:
        print('graph 822 is Eulerian')
    else:
        print('graph822 is not Eulerian')

    # Question 2
    if has_eulerian_path(graph_822) == True:
        print('graph_822 has a Eulerian path')
    else:
        print('graph_822 does not have a Eulerian path')

    # Question 3
    cycle_822 = find_eulerian_cycle(graph_822)
    
    print('\neulerian cycle for graph_822:\n',\
          combine(graph_822,cycle_822))
   
    # Question 5
    print('\ngraph from fig.8.20:')
    for k, v in create_graph(s).items():
        print(k,v)
    
   # Question 6
    graph_s = create_graph(s)
    if is_eulerian(graph_s) == True:
        print('\ngraph s is Eulerian')
    else:
        print('\ngraph s is not Eulerian')
        
    if has_eulerian_path(graph_s) == True:
        print('graph_s has a Eulerian path')
    else:
        print('graph_s does not have a Eulerian path')
    
    # Question 7
    path_s = find_eulerian_path(graph_s)
    print('\npath of s:\n',path_s)

    seq = path_s[0][0]
    for edge in path_s[1:-1]:
        seq += edge[0][1]
    seq += path_s[-1][-1]
    print('\nDNA sequence:', seq)
        
    # Qusetion 8
    paths = find_eulerian_cycle(bigger_graph)
    print('\nEulerian cycle of bigger graph:\n',
          combine(bigger_graph,paths))





