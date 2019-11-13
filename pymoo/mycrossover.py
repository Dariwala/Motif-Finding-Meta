from pymoo.model.crossover import Crossover
from pymoo.operators.crossover.util import crossover_mask
import numpy as np
from random import random,randrange

def computeConsensus(solution,seqs,l_mer):
    consensus = ""
    j = -1

    for i in range(l_mer):
        j += 1
        a = c = g = t = 0
        for k in range(len(solution)):
            #print(solution)
            #print(int(solution[k])+j)
            if seqs[k][int(solution[k])+j] == 'A':
                a += 1
            elif seqs[k][int(solution[k])+j] == 'C':
                c += 1
            elif seqs[k][int(solution[k])+j] == 'G':
                g += 1
            else:
                t += 1
        if a == max(max(a,c),max(g,t)):
            consensus += 'A'
        elif c == max(max(a,c),max(g,t)):
            consensus += 'C'
        elif g == max(max(a,c),max(g,t)):
            consensus += 'G'
        elif t == max(max(a,c),max(g,t)):
            consensus += 'T'
    print(consensus)
    return consensus

def hammingDistance(str1,str2):
    count = 0
    for i in range(len(str1)):
        if str1[i] == str2[i]:
            count += 1
    return count/len(str1)

class MyCrossover(Crossover):

    def __init__(self, prob=0.5, **kwargs):
        super().__init__(2, 2, **kwargs)
        self.prob = prob

    def _do(self, problem, X, **kwargs):
        #print(X.shape)
        if self.prob is None:
            self.prob = 1/len(X[0])
        _X = np.copy(X)
        for i in range(len(_X[0])):
            p = random()
            if p <= self.prob:
                index = randrange(len(problem.seqs))
                aaaa = _X[0][i][index]
                _X[0][i][index] = _X[1][i][index]
                _X[1][i][index] = aaaa

        return _X
