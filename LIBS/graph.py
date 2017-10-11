#!/usr/bin/python

"""
a library built to define graphs.

Here I use the adjacent list implementation to represent weighted directional graphs
(+) class Vertex defining a Vertex
(+) class Graph to built up a Graph (Vertices + Edges)
"""

class Vertex:
    def __init__(self,key):
        self.id = key
        self.connectedTo = {}

    def addNeighbor(self,nbr,weight=0):
        self.connectedTo[nbr] = weight

    def __str__(self):
        return str(self.id) + ' connectedTo: ' + str([x.id for x in self.connectedTo])

    def getConnections(self):
        return self.connectedTo.keys()

    def getId(self):
        return self.id

    def getWeight(self,nbr):
        return self.connectedTo[nbr]


class Graph:
    def __init__(self):
        self.vertList = {}    # vertList: dictionary
        self.numVertices = 0  # number of vertices

    def addVertex(self,key):
        self.numVertices = self.numVertices + 1
        newVertex = Vertex(key)         # initialize new vertex
        self.vertList[key] = newVertex  # add an entry (key) to dictionary
        return newVertex

    def getVertex(self,n):
        if n in self.vertList:
            return self.vertList[n]
        else:
            return None

    def __contains__(self,n):
        return n in self.vertList

    def addEdge(self,f,t,cost=0):
        if f not in self.vertList: # if vertex1 of edge not yet in vertList --> add it
            nv = self.addVertex(f)
        if t not in self.vertList: # if vertex2 if edge bit yet in vertList --> add it
            nv = self.addVertex(t)

        # add directional edge from f --> t
        self.vertList[f].addNeighbor(self.vertList[t], cost)

    def getVertices(self):
        return self.vertList.keys()

    def __iter__(self):
        return iter(self.vertList.values())
