  # AO PREENCHER ESSE CABEÇALHO COM O MEU NOME E O MEU NÚMERO USP,
  # DECLARO QUE SOU O ÚNICO AUTOR E RESPONSÁVEL POR ESSE PROGRAMA.
  # TODAS AS PARTES ORIGINAIS DESSE EXERCÍCIO-PROGRAMA (EP) FORAM
  # DESENVOLVIDAS E IMPLEMENTADAS POR MIM SEGUINDO AS INSTRUÇÕES
  # DESSE EP E QUE PORTANTO NÃO CONSTITUEM DESONESTIDADE ACADÊMICA
  # OU PLÁGIO.
  # DECLARO TAMBÉM QUE SOU RESPONSÁVEL POR TODAS AS CÓPIAS DESSE
  # PROGRAMA E QUE EU NÃO DISTRIBUI OU FACILITEI A SUA DISTRIBUIÇÃO.
  # ESTOU CIENTE QUE OS CASOS DE PLÁGIO E DESONESTIDADE ACADÊMICA
  # SERÃO SEVERAMENTE PUNIDOS.
  # ENTENDO QUE EPS SEM DECLARAÇÃO DE RESPONSABILIDADE NÃO SERÃO
  # CORRIGIDOS E PODERÃO CAUSAR PUNIÇÃO DO AUTOR POR DESONESTIDADE
  # ACADÊMICA.
  #
  # Nome: Isabela Blucher
  # Número USP: 9298170
  # Curso: Biologia Computacional - MAC0465         Prof. Alair
  # Exercício-Programa 3
  #
  #  OBS: o algoritmo de ordenação topológica tem como fonte: https://www.geeksforgeeks.org/topological-sorting/

# -*- coding: utf-8 -*-

import sys
import unittest
import numpy as np
from collections import defaultdict

# devolve o overlap entre dois fragmentos
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


# constrói o grafo de overlap dos fragmentos
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


# função recursiva para a ordenação topológica
def topologicalSortRec(graph, v, visited, stack):
    visited[v] = True

    for i in graph[v]:
        if not visited[i]:
            topologicalSortRec(graph, i, visited, stack)

    stack.insert(0, v)


# chamada principal para a ordenação topológica
def topologialSort(graph, V):
    visited = [False]*V
    stack = []

    for i in range(V):
        if not visited[i]:
            topologicalSortRec(graph, i, visited, stack)

    return stack


# imprime os vértices do grafo
def printVertex(F):
    print('Vértices:')
    for i in range(len(F)):
        print(i, F[i])
    print('\n')


# imprime as arestas do grafo
def printEdges(weights):
    print('Arestas:')
    for k, v in weights.items():
        a, b = k
        print('{} -> {} ^ {}'.format(a, b, v))
    print('\n')


# imprime o layout dos fragmentos
def printFragmentLayout(F, order, positions, super_len):
    for i in order:
        cnt = 0
        s = ''
        while cnt != super_len:
            if cnt == positions[i]:
                s += F[i]
                cnt += len(F[i])
            else:
                s += '-'
                cnt += 1
        print(s, i)


# Exemplo de chamada do código:
# python ep3.py 3 AGTATTGGCAATC AATCGATG ATGCAAACCT CCTTTTGG TTGGCAATCACT
def main():
    t = int(sys.argv[1])
    F = sys.argv[2:]
    F = np.array(F)
    V = len(F)

    graph, weights, overlaps = buildGraph(F, t)

    order = topologialSort(graph, V)

    printVertex(F)
    printEdges(weights)

    superseq = F[order[0]]
    init_position = []
    init_position.insert(order[0], 0)
    pos = 0
    for s, t in zip(order, order[1:]):
        w = weights[s, t]
        superseq += F[t][w:]

        pos += len(F[s]) - w
        init_position.insert(t, pos)

    print('Supersequência comum e Layout dos fragmentos:')
    print(superseq)
    printFragmentLayout(F, order, init_position, len(superseq))
    print(superseq)

main()
