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
  # Exercício-Programa 2
  #
  # - O algoritmo "similarity"  é uma generalização do algoritmo do professor
  #   no arquivo sim3.py, postado no PACA da disciplina


import sys
import numpy as np
from itertools import combinations
from itertools import product

memo = {}
k = 0
r = 0
q = 0
g = 0
alignment = []

# verfica se as sequências estão vazias para o algoritmo similarity
def check_sequence_lengths_sim(last_index):
    for i in last_index:
        if i < 0:
            return 0
    return 1

# verifica se pelo menos uma sequência não está vazia para o algoritmo backtrace
def check_sequence_lengths_back(last_index):
    if all(v == 0 for v in last_index) == False:
        return 1
    else:
        return 0

# retorna uma lista com os índices dos últimos elementos das sequências
def get_last_index(seqs):
    index = []
    for s in seqs:
        index.append(len(s))
    return index

# faz match ou mismatch dos elementos das sequências
def p(a, b):
    if a==b: return r
    else: return q

# calcula o sp-measure de k sequências
def sp_measure(seqs, delta, last_index):
    seqs_with_1 = []
    for i in range(len(delta)):
        if delta[i] == 1:
            seqs_with_1.append(i)

    comb = combinations(delta, 2)
    sp = 0
    for c in comb:
        a, b = c
        if a == 1 and b == 0 or a == 0 and b == 1:
            sp += g

    if len(seqs_with_1) > 1:
        comb_seqs = combinations(seqs_with_1, 2)
        for c in comb_seqs:
            a, b = c
            score = p(seqs[a][last_index[a]], seqs[b][last_index[b]])
            sp += score

    return sp

# computa o score do alinhamento das k sequências
def similarity(seqs, last_index):
    global memo

    if tuple(last_index) in memo: return memo[tuple(last_index)]
    if check_sequence_lengths_sim(last_index) == 0: return float('-inf')

    delta = list(product([0, 1], repeat=k))
    delta = delta[1:]
    m = float('-inf')

    for i in range(2**k - 1):
        d = delta[i]
        d = np.array(d)

        # remover 1 dos last index de quem tem 1
        last_index = last_index - d

        # pegar o máximo m que ocorreu em certo delta
        m = max(m, similarity(seqs, last_index) + sp_measure(seqs, d, last_index))
        last_index = last_index + d


    memo[tuple(last_index)] = m
    return m

# retorna um alinhamento ótimo de k sequências
def backtrace(seqs, m):
    global alignment

    last_index = get_last_index(seqs)
    score = m[tuple(last_index)]

    delta = list(product([0, 1], repeat=k))
    delta = delta[1:]

    while check_sequence_lengths_back(last_index):
        for i in range(2**k - 1):
            d = delta[i]
            d = np.array(d)

            ind = last_index - d
            align = []
            if all(v >= 0 for v in ind):
                if score == m[tuple(ind)] + sp_measure(seqs, d, ind):
                    for j in range(len(d)):
                        if d[j] == 1:
                            align.append(seqs[j][ind[j]])
                        else:
                            align.append('-')

                    alignment.append(align)
                    align = ' '.join(align)

                    last_index = last_index - d

        score = m[tuple(last_index)]

    return alignment

def main():
    global memo
    global k
    global r
    global q
    global g

    k = int(sys.argv[1])
    r, q, g = sys.argv[2:5]
    r = int(r)
    q = int(q)
    g = int(g)

    memo[(0,) * k] = 0

    seqs = []
    for i in range(k):
        s = sys.argv[5 + i]
        seqs.append([x for x in s])

    last_index = np.array(get_last_index(seqs))
    p = similarity(seqs, last_index)

    b = backtrace(seqs, memo)
    b = np.array(b)
    b = b[::-1]
    b = b.T

    for i in range(len(b)):
        print(' '.join(b[i]))

main()
