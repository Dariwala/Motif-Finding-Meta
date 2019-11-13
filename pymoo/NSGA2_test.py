from problem import MotifFinding
from mymutation import MyMutation,NayeemsirMutation,RouletteWheelMutation
from mycrossover import MyCrossover
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.operators.crossover.simulated_binary_crossover import SimulatedBinaryCrossover
from pymoo.operators.mutation.polynomial_mutation import PolynomialMutation
from pymoo.factory import get_performance_indicator
import numpy as np
import matplotlib.pyplot as plt



input_file = "dm01g.txt"
seqs = open(input_file,'r')
seqs = seqs.readlines()
sequences = []
no_of_objs = 2

for line in seqs:
    sequences.append(line[:-1])
for l_mer in [18]:
    reference_points = np.ones((1,no_of_objs))*1.2
    #print(reference_points)
    motif_finding = MotifFinding(n_var = len(sequences),n_obj = no_of_objs,xl=0,xu = len(sequences[0])-l_mer,seqs = sequences,l_mer = l_mer)

    algorithm = NSGA2(pop_size=200,
                      elimate_duplicates=True,
                      mutation=MyMutation(prob=0.6)
                      #,crossover=SimulatedBinaryCrossover(prob = 0.9,eta = 15)
                      ,crossover=MyCrossover(prob=0.9)
                      )


    hv = get_performance_indicator("hv", ref_point=reference_points[0])

    hvs = []
    results = []
    for r in range(10):
        res = minimize(motif_finding,
                       algorithm,
                       ('n_gen', 1000),
                       seed=1,
                       verbose=False
                       , save_history=True)

        # plt.figure(1)
        # plt.clf()
        '''gens = []
        for g, a in enumerate(res.history):
            if g % 20 == 0 and g>50:
                a.opt = a.pop.copy()
                #a.opt = a.opt[a.opt.collect(lambda ind: ind.feasible)[:, 0]]
                #I = NonDominatedSorting().do(a.opt.get("F"),
                #                             only_non_dominated_front=True)
                #a.opt = a.opt[I]
                X, F, CV, G = a.opt.get("X", "F", "CV", "G")
                hvs.append(hv.calc(F))
                gens.append(g)
                #print("hv", hv.calc(A))
    
        plt.plot(gens,hvs)
        plt.show()'''
        print("Best solution found:F = %s" % (np.min(res.F[:,0])))
        results.append((np.min(res.F[:,0]),np.min(res.F[:,1])))
    #output = open(parent + str(l_mer) + ".txt", "w")
    #for r in results:
    #    output.write(str(-r[0]) + "\t" + str(-r[1]) + "\n")
