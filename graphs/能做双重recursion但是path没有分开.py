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



def get_next_edge(graph, node):
    candidate_nodes = graph[node]
    return node + candidate_nodes[0]


def find_subcycle(graph, node, start_node, used_edges, used_nodes):
    next_node = graph[node][0]
    used_nodes.append(node)
    choosen_edge = get_next_edge(graph, node)
    used_edges.append(choosen_edge)   
    if len(graph[node]) == 1:
        graph.pop(node)
    else:
        graph[node].remove(next_node)
    if next_node != start_node:
        return find_subcycle(graph, next_node, start_node, used_edges, used_nodes)
    else:
        return used_edges, graph, used_nodes
        print(used_edges)



def find_eulerian_cycle(graph, paths,used_nodes):
    if graph != {}:
        start_node = list(graph.keys())[0]
        search_sub_cycle = find_subcycle(graph, start_node,start_node, used_edges, used_nodes)
        paths.append(search_sub_cycle[0])
        print(paths)
        used_nodes = search_sub_cycle[2]
        graph = search_sub_cycle[1]
        return find_eulerian_cycle(graph, paths, used_nodes)
    else:
        
        return paths

    


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


    used_edges = []
    used_nodes = []
    paths = []
    print(find_eulerian_cycle(graph_822, paths, used_nodes))
    #print(find_subcycle(graph_822, 'A', 'A', used_edges, used_nodes))





    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)
