from rdkit import Chem
import numpy as np
import time
import crossover as co
import scoring_functions as sc
import GB_GA as ga 
import sys
from multiprocessing import Pool
import pickle

scoring_function = sc.absorbance_target
target = 600. # nm
sigma = 50. # nm
threshold = 0.3
n_confs = 20
xtb_path = '/home/jhjensen/xtb4stda-1.0'
scoring_args = [n_confs, xtb_path, target, sigma, threshold]
max_score = 1.99

n_tries = 40 
n_cpus = n_tries
population_size = 20 
mating_pool_size = 20
generations = 50
mutation_rate = 0.05
co.average_size = 50. 
co.size_stdev = 5.
prune_population = True
seeds = np.random.randint(1_000_000, size=2*n_tries)

file_name = sys.argv[1]

print('* target', target)
print('* sigma', sigma)
print('* threshold', threshold)
print('* n_confs', n_confs)
print('* max_score', max_score)
print('* population_size', population_size)
print('* mating_pool_size', mating_pool_size)
print('* generations', generations)
print('* mutation_rate', mutation_rate)
print('* average_size/size_stdev', co.average_size, co.size_stdev)
print('* initial pool', file_name)
print('* prune population', prune_population)
print('* number of tries', n_tries)
print('* number of CPUs', n_cpus)
print('* seeds', ','.join(map(str, seeds)))
print('* ')
print('run,score,smiles,generations,representation,prune')

high_scores_list = []
count = 0

temp_args = [[population_size, file_name,scoring_function,generations,mating_pool_size,
                mutation_rate,scoring_args, max_score, prune_population] for i in range(n_tries)]
args = []
for x,y in zip(temp_args,seeds):
    x.append(y)
    args.append(x)
with Pool(n_cpus) as pool:
    output = pool.map(ga.GA, args)

for i in range(n_tries):   
    #(scores, population) = ga.GA([population_size, file_name,scoring_function,generations,mating_pool_size,mutation_rate,scoring_args],prune_population)
    (scores, population, high_scores, generation) = output[i]
    smiles = Chem.MolToSmiles(population[0], isomericSmiles=True)
    high_scores_list.append(high_scores)
    print(f'{i},{scores[0]:.2f},{smiles},{generation},Graph,{prune_population}')

pickle.dump(high_scores_list, open('GA_absorbance_600.p', 'wb' ))
