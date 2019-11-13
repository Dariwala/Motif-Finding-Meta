import numpy as np

from pymoo.model.mutation import Mutation
from random import randrange
from random import random

def hammingDistance(str1,str2):
    count = 0
    for i in range(len(str1)):
        if str1[i] == str2[i]:
            count += 1
    return count/len(str1)

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
    #print(consensus)
    return consensus



class NayeemsirMutation(Mutation):

    def __init__(self, prob=None):
        super().__init__()
        self.prob = prob

    def _do(self, problem, X, **kwargs):
        if self.prob is None:
            self.prob = 1.0 / problem.n_var
        int_prob = 1.0 / problem.n_var
        for i in range(len(X)):
            for j in range(problem.n_var):
                take_seq = random()
                if take_seq <= int_prob:
                    if X[i][j] != len(problem.seqs[0]) - problem.l_mer and X[i][j] != 0:
                        p = random()
                        if p>=0.5:
                            X[i][j] += 1
                        else:
                            X[i][j] -=1
                    elif X[i][j] == 0:
                        X[i][j] +=1
                    else:
                        X[i][j] -= 1
        return X

class MyMutation(Mutation):

    def __init__(self, prob=None):
        super().__init__()
        self.prob = prob

    def _do(self, problem, X, **kwargs):
        #print(X.shape)
        if self.prob is None:
            self.prob = 1.0 / problem.n_var
        for i in range(len(X)):
            value = randrange(len(problem.seqs[0]) - problem.l_mer + 1)#problem.xu + 1
            #print(value)
            p = random()
            if p<= self.prob:
                X[i][randrange(problem.n_var)] = value
        return X


class RouletteWheelMutation(Mutation):

    def __init__(self, prob=None):
        super().__init__()
        self.prob = prob

    def _do(self, problem, X, **kwargs):
        if self.prob is None:
            self.prob = 1.0 / problem.n_var

        for i in range(len(X)):
            consensus = computeConsensus(X[i],problem.seqs,problem.l_mer)
            probabilities = []
            sum = 0.0
            for j in range(len(X[i])):
                probabilities.append(hammingDistance(consensus,problem.seqs[j][int(X[i][j]):int(X[i][j])+problem.l_mer]))
                sum += hammingDistance(consensus,problem.seqs[j][int(X[i][j]):int(X[i][j])+problem.l_mer])
            for p in probabilities:
                p/= sum
            for j in range(1,len(probabilities)):
                probabilities[j] += probabilities[j-1]
            p = random()
            for j in range(len(probabilities)):
                if p<= probabilities[j]:
                    index = j
                    break
            value = randrange(len(problem.seqs[0]) - problem.l_mer + 1)#problem.xu + 1
            #print(value)
            p = random()
            if p<= self.prob:
                #X[i][randrange(problem.n_var)] = value
                X[i][index] = value
        return X