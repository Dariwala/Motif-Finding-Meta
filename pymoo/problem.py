from pymoo.model.problem import Problem
import autograd.numpy as anp
import numpy as np
import math

def similarity(str1,str2):
	count = 0
	for i in range(len(str1)):
		if str1[i] == str2[i]:
			count += 1
	return count/len(str1)

class MotifFinding(Problem):
	def __init__(self, n_var, n_obj,xl,xu,seqs,l_mer,**kwargs):
		self.seqs = seqs
		self.l_mer = l_mer
		super().__init__(n_var=n_var, n_obj=n_obj, xl=xl, xu=xu, type_var=int,**kwargs)

	def _evaluate(self, x, out, *args, **kwargs):
		f1 = np.zeros((len(x),1),dtype=np.double)
		f2 = np.zeros((len(x),1),dtype=np.double)
		f3 = np.zeros((len(x), 1), dtype=np.double)
		for i in range(len(x)):
			subseqs = []
			for j in range(len(x[0])):
				subseqs.append(self.seqs[j][int(x[i][j]):int(x[i][j])+self.l_mer])
			#print(len(subseqs),len(subseqs[0]))
			consensus = ""
			for j in range(self.l_mer):
				a = 0
				c = 0
				g = 0
				t = 0
				for k in range(len(subseqs)):
					#print(k,j)
					try:
						if subseqs[k][j] == "A":
							a+=1
						elif subseqs[k][j] == "C":
							c+=1

						elif subseqs[k][j] == "G":
							g+=1
						else:
							t+=1
					except IndexError:
						print(subseqs)
						print(k,j,len(subseqs))
				f1[i][0] += max(max(a,c),max(g,t))
				if a!=0:
					f3[i][0] += a/(a+c+g+t) * math.log(a/(a+c+g+t),2)
				if c!=0:
					f3[i][0] += c/(a+c+g+t) * math.log(c/(a+c+g+t),2)
				if g!=0:
					f3[i][0] += g/(a+c+g+t) * math.log(g/(a+c+g+t),2)
				if t!=0:
					f3[i][0] += t/(a+c+g+t) * math.log(t/(a+c+g+t),2)
				if a == max(max(a,c),max(g,t)):
					consensus += 'A'
				elif c == max(max(a,c),max(g,t)):
					consensus += 'C'
				elif g == max(max(a,c),max(g,t)):
					consensus += 'G'
				elif t == max(max(a,c),max(g,t)):
					consensus += 'T'
			f1[i][0] /= self.l_mer * len(self.seqs)
			f3[i][0] /= self.l_mer * 4
			#print(f1[i])
			for k in range(len(subseqs)):
				if similarity(subseqs[k],consensus)>=0.6:
					f2[i][0] += 1
			f2[i][0] /= len(self.seqs)
		out["F"] = np.column_stack([-f1,-f3])
		#print("lol")
			