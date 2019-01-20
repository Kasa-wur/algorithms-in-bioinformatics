#!/usr/bin/env python

"""
Author: Xinyuan Min
Student nr: 950829573070
Script to:
"""
# Import statements


# Implement your functions here
def is_eulerian(graph):
    # use set() to get every distinct node in graph
    vertices = set()
    for key_value in graph.items():
        vertices.add(key_value[0])
        for value in key_value[1]:
            vertices.add(value)

    count_in_out = {}
    for vertex in vertices:
        # the first element of value of a vertex indicates its indegree
        # similarly, the second element represents a vertex's outdegree
        count_in_out[vertex] = [0, 0]

    # the length of value(list) of each node(key) is the outdegree of a vertex
    for node in graph:
        outdegree = len(graph[node])
        count_in_out[node][1] = outdegree
        # indegree is the number of a node occurs as a value of other nodes
        for connected_node in graph[node]:
            count_in_out[connected_node][0] += 1

    # a connected graph is Eulerian if all vertices is balanced
    eulerian = True
    for vertex in count_in_out:
        indegree = count_in_out[vertex][0]
        outdegree = count_in_out[vertex][1]
        if indegree != outdegree:
            eulerian = False
            break

    return eulerian, count_in_out



def has_eulerian_path(graph):
    count_in_out = is_eulerian(graph)[1]
    semi_balanced = 0
    eulerian_path = False
    for vertex in count_in_out:
        indegree = count_in_out[vertex][0]
        outdegree = count_in_out[vertex][1]
        if abs(indegree - outdegree) == 1:
            semi_balanced += 1
    if semi_balanced <= 2:
        eulerian_path = True

    if eulerian_path == True:
        print('This graph has a Eulerian path')
    else:
        print('This graph does not have a Eulerian path')

    return eulerian_path


def get_edges(graph):
    """Return a list of strings which contains all the edges in the graph

    :param a:
    :param node:
    :return:
    """
    edges = []
    for vertex in graph.keys():
        connected_nodes = graph[vertex]
        for node in connected_nodes:
            edges.append(str(vertex + node))

    return edges

def get_next_edge(graph, node):
    candidate_nodes = graph[node]
    return node + candidate_nodes[0]

unused_edges = get_edges(graph_822)
used_edges = []


def find_eulerian_cycle(graph, unused_edges, node, used_edges):
    next_node = graph[node][0]
    choosen_edge = get_next_edge(graph, node)
    used_edges.append(choosen_edge)
    unused_edges.remove(choosen_edge)
    if next_node != 'A':
        return find_eulerian_cycle(graph, unused_edges, next_node, used_edges)
    else:
        return used_edges, unused_edges


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
    if is_eulerian(graph_822)[0] == True:
        print('this graph is Eulerian')
    else:
        print('This graph is not Eulerian')

    # Question 2
    has_eulerian_path(graph_822)

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
