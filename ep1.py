### EP1 de MAC0465 - Alinhador local de sequências

### Nome: Isabela Blucher
### NUSP: 9298170

import numpy as np
import math

# Pontuação da similaridade do alinhamento
match= 1
mismatch= -1
gap= -2

# Incialização de variáveis globais
pref = np.zeros(100)
suff = np.zeros(100)
align_s = [None] * 100
align_t = [None] * 100

# Retorna match ou mismatch para a comparação entre caracteres de s e t
def p(a,b):
    if a == b:
        return match
    return mismatch

# Retorna as strings cortadas nas extremidades associadas ao alinhamento local ótimo
def getStringExtremeties(s, t):
    global pref
    best_score, best_i, best_j = (0, 0, 0)

    for i in range(len(s) + 1):
        for j in range(len(t) + 1):
            tmp = pref[j]
            if i == 0 or j == 0:
                pref[j] = 0
            else:
                pref[j] = max(pref[j - 1] + gap, \
                old + p(s[i - 1], t[j - 1]), \
                pref[j] + gap, \
                0)

            if pref[j] > best_score:
                best_score = pref[j];
                best_i = i
                best_j = j
            old = tmp

    while len(s) > best_i:
        s = s[:-1]
    while len(t) > best_j:
        t = t[:-1]

    return s, t

# Retorna as strings cortadas nas extremidades associadas ao alinhamento local ótimo para as strings reversas
def getReverseStringExtremeties(s, t):
    rev_s = s[::-1]
    rev_t = t[::-1]
    rev_s, rev_t =  getStringExtremeties(rev_s, rev_t)
    s = rev_s[::-1]
    t = rev_t[::-1]
    return s, t

# Retorna as substrings das sequências originais que serão usadas no alinhamento
def localSequences(s, t):
    s, t = getStringExtremeties(s, t)
    return getReverseStringExtremeties(s, t)

# Pontuação de similaridade entre sequências
def bestScore(s, t, a, b, c, d):
    global pref
    for j in range(c - 1, d + 1):
        pref[j] = (j - c + 1) * gap

    for i in range(a, b + 1):
        old = pref[c - 1]
        pref[c - 1] = (i - a + 1) * gap
        for j in range(c, d):
            tmp = pref[j]
            pref[j] = max(pref[j - 1] + gap, \
                old + p(s[i], t[j]), \
                pref[j] + gap)
            old = tmp


# Pontuação de similaridade entre sequências reversas
def bestScoreReverse(s, t, a, b, c, d):
    global suff
    for j in range(d + 1, c - 1, -1):
        suff[j] = (d - j + 1) * gap

    for i in range(b, a - 1, -1):
        old = suff[d + 1]
        suff[d + 1] = (b - i + 1) * gap
        for j in range(d, c - 1, -1):
            tmp = suff[j]
            suff[j] = max(suff[j - 1] + gap, \
                old + p(s[i - 1], t[j - 1]), \
                suff[j] + gap)
            old = tmp

# Construção do alinhamento local ótimo
# Algoritmo retirado do livro Introduction to Computational Molecular Biology (algoritmo da figura 3.6, página 61)
def alignSequences(s, t, a, b, c, d, start, end):
    global pref
    global suff
    global align_s
    global align_t

    # Caso base da recursão do algoritmo
    if a > b or c > d:
        while a <= b:
            align_s.append(s[a - 1])
            align_t.append('-')
            a += 1
        while c <= d:
            align_s.append('-')
            align_t.append(t[c - 1])
            c += 1
        end = start + max(len(s), len(t))

    else:
        i = math.floor((a + b) / 2)
        bestScore(s, t, a, i - 1, c, d)
        bestScoreReverse(s, t, i + 1, b, c, d)

        posmax = c - 1
        typemax = '-'
        vmax = pref[c - 1] + gap + suff[c - 1]
        middle = i

        for j in range(c, d + 1):
            if pref[j - 1] + p(s[i - 1], t[j - 1]) + suff[j + 1] > vmax:
                posmax = j
                typemax = 'symbol'
                vmax = pref[j - 1] + p(s[i - 1], t[j - 1]) + suff[j + 1]
            if pref[j] + gap + suff[j + 1] > vmax:
                posmax = j
                typemax = '-'
                vmax = pref[j] + gap + suff[j + 1]

        if typemax == '-':
            alignSequences(s, t, a, i - 1, c, posmax, start, middle)
            align_s[middle] = s[i - 1]
            align_t[middle] = '-'
            alignSequences(s, t, i + 1, b, posmax + 1, d, middle + 1, end)
        else:
            alignSequences(s, t, a, i - 1, c, posmax - 1, start, middle)
            align_s[middle] = s[i - 1]
            align_t[middle] = t[posmax - 1]
            alignSequences(s, t, i + 1, b, posmax + 1, d, middle + 1, end)


def main():
    global pref
    global suff
    global align_s
    global align_t

    s = input('s = ')
    t = input('t = ')

    s, t = localSequences(s, t)

    alignSequences(s, t, 1, len(s), 1, len(t), 0, len(s) + len(t))

    align_s = [char for char in align_s if char is not None]
    align_t = [char for char in align_t if char is not None]

    print('Um alinhamento local ótimo para as sequências s e t é: ')
    print(''.join(align_s))
    print(''.join(align_t))

main()
