import sys
import unittest
import numpy as np
from collections import defaultdict

def getOverlap(left, right):
    overlap = ""
    for c in right:
        temp  = overlap + c
        if left.find(temp) >= 0:
            overlap = temp
        else:
            if left.endswith(overlap) and right.startswith(overlap):
                return overlap
    return


def buildGraph(F, t):
    graph = defaultdict(list)
    weights = {}
    overlaps ={}

    for i in range(len(F)):
        for j in range(len(F)):
            if i != j:
                overlap = getOverlap(F[i], F[j])
                if overlap and len(overlap) >= t:
                    graph[i].append(j)
                    weights[i, j] = len(overlap)
                    overlaps[i, j] = overlap

    return graph, weights, overlaps


def topologicalSortRec(graph, v, visited, stack):
    visited[v] = True

    for i in graph[v]:
        if not visited[i]:
            topologicalSortRec(graph, i, visited, stack)

    stack.insert(0, v)


def topologialSort(graph, V):
    visited = [False]*V
    stack = []

    for i in range(V):
        if not visited[i]:
            topologicalSortRec(graph, i, visited, stack)

    print(stack)


def main():
    t = int(sys.argv[1])
    F = sys.argv[2:]
    F = np.array(F)
    V = len(F)

    graph, weights, overlaps = buildGraph(F, t)
    print(overlaps)

    topologialSort(graph, V)

main()
