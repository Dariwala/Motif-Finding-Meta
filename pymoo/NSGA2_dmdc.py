from problem import MotifFinding
from mymutation import MyMutation,NayeemsirMutation
from mycrossover import MyCrossover
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.operators.crossover.simulated_binary_crossover import SimulatedBinaryCrossover
from pymoo.operators.mutation.polynomial_mutation import PolynomialMutation
from pymoo.factory import get_performance_indicator
import numpy as np
import matplotlib.pyplot as plt
import os,sys

dir = "NSGA-DefaultMutation-DefaultCrossover\\"

datasets = ["hm03","hm09g","yst04r","yst08r","dm01g","dm03g"]

if os.path.exists(dir) == False:
    os.makedirs(dir)
if os.path.exists(dir + "Hv-data\\") == False:
    os.makedirs(dir+ "Hv-data\\")

if os.path.exists(dir + "Best Outputs\\") == False:
    os.makedirs(dir+ "Best Outputs\\")
if os.path.exists(dir + "Pareto front\\") == False:
    os.makedirs(dir+ "Pareto front\\")

for dataset in datasets:
    if os.path.exists(dir+"Hv-data\\"+dataset)==False:
        os.makedirs(dir+"Hv-data\\"+dataset)
    if os.path.exists(dir+"Best Outputs\\"+dataset)==False:
        os.makedirs(dir+"Best Outputs\\"+dataset)
    if os.path.exists(dir+"Pareto front\\"+dataset)==False:
        os.makedirs(dir+"Pareto front\\"+dataset)

name = sys.argv[1]

input_file = name+".txt"
seqs = open(input_file,'r')
seqs = seqs.readlines()
sequences = []
no_of_objs = 2

for line in seqs:
    sequences.append(line[:-1])
parent = dir + "Best Outputs\\"+name+"\\"
for l_mer in [16,18,20,22]:
    reference_points = np.ones((1,no_of_objs))*1.2
    #print(reference_points)
    motif_finding = MotifFinding(n_var = len(sequences),n_obj = no_of_objs,xl=0,xu = len(sequences[0])-l_mer,seqs = sequences,l_mer = l_mer)

    algorithm = NSGA2(pop_size=200,
                      elimate_duplicates=True
                      #,mutation=MyMutation(prob=0.6)
                      #,crossover=MyCrossover(prob = 0.9)
                      )


    hv = get_performance_indicator("hv", ref_point=reference_points[0])

    hvs = []
    results = []
    for r in range(20):
        volume = open(dir+"Hv-data\\"+name+"\\run"+str(r)+".txt",'w')
        pareto = open(dir + "Pareto front\\" + name + "\\run" + str(r) + ".txt", 'w')
        res = minimize(motif_finding,
                       algorithm,
                       ('n_gen', 1000),
                       seed=1,
                       verbose=False
                       , save_history=True)

        # plt.figure(1)
        # plt.clf()
        gens = []
        for g, a in enumerate(res.history):
            if g % 40 == 0:
                a.opt = a.pop.copy()
                #a.opt = a.opt[a.opt.collect(lambda ind: ind.feasible)[:, 0]]
                #I = NonDominatedSorting().do(a.opt.get("F"),
                #                             only_non_dominated_front=True)
                #a.opt = a.opt[I]
                X, F, CV, G = a.opt.get("X", "F", "CV", "G")
                hvs.append(hv.calc(F))
                gens.append(g)

        for v in hvs:
            volume.write(str(v)+"\n")
        for re in res.F:
            pareto.write(str(-re[0])+"\t"+str(-re[1])+"\n")
        volume.close()
        pareto.close()
    
        #plt.plot(gens,hvs)
        #plt.show()
        print(name,l_mer,r,np.min(res.F[:,0]))
        results.append((np.min(res.F[:,0]),np.min(res.F[:,1])))
    output = open(parent + str(l_mer) + ".txt", "w")
    for r in results:
        output.write(str(-r[0]) + "\t" + str(-r[1]) + "\n")
    output.close()
