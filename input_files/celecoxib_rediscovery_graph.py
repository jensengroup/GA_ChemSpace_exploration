from rdkit import Chem
import numpy as np
import time
import crossover as co
import scoring_functions as sc
import GB_GA as ga 
import sys
from multiprocessing import Pool
import pickle

Celecoxib = 'O=S(=O)(c3ccc(n1nc(cc1c2ccc(cc2)C)C(F)(F)F)cc3)N'
target = Chem.MolFromSmiles(Celecoxib)
n_tries = 40
population_size = 100
mating_pool_size = 200
generations = 1000
mutation_rate = 0.5
co.average_size = target.GetNumAtoms() 
co.size_stdev = 5
scoring_function = sc.rediscovery
max_score = 1.0
scoring_args = [target]
n_cpus = n_tries
seeds = np.random.randint(1_000_000, size=2*n_tries)

file_name = sys.argv[1]

print('* target', Celecoxib)
print('* population_size', population_size)
print('* mating_pool_size', mating_pool_size)
print('* generations', generations)
print('* mutation_rate', mutation_rate)
print('* max_score', max_score)
print('* average_size/size_stdev', co.average_size, co.size_stdev)
print('* initial pool', file_name)
print('* number of tries', n_tries)
print('* number of CPUs', n_cpus)
print('* seeds', ','.join(map(str, seeds)))
print('*')
print('run,score,smiles,generations,representation,prune')

high_scores_list = []
count = 0
for prune_population in [True]:#, False]:
    index = slice(0,n_tries) if prune_population else slice(n_tries,2*n_tries)
    temp_args = [[population_size, file_name,scoring_function,generations,mating_pool_size,
                  mutation_rate,scoring_args, max_score, prune_population] for i in range(n_tries)]
    args = []
    for x,y in zip(temp_args,seeds[index]):
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

pickle.dump(high_scores_list, open('celecoxib_rediscovery_10K.p', 'wb' ))