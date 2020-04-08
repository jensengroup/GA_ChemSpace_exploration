'''
Written Jan H. Jensen 2019 
Requires import of string_crossover and string_mutate in GB_GA
'''

from rdkit import Chem
import numpy as np
import time

import string_scoring_functions as sc
import string_crossover as co
import string_GA as ga 
import sys
from multiprocessing import Pool
import pickle

from rdkit import rdBase

file_name = sys.argv[1]

target = 'O=C1NC(=O)SC1Cc4ccc(OCC3(Oc2c(c(c(O)c(c2CC3)C)C)C)C)cc4'
target = Chem.MolFromSmiles(target)
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
seeds = np.random.randint(1_000_000, size=3*2*n_tries)

print('* RDKit version', rdBase.rdkitVersion)
print('* target', target)
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
print('* ')
print('run,score,smiles,generations,representation,prune')

high_scores_list = []
count = 0
for s_type in ['SMILES', 'DeepSMILES', 'SELFIES']:
    co.string_type = s_type
    line_list = [f'{co.string_type}']
    for prune_population in [True, False]:
        args =  [population_size, file_name,scoring_function,generations,
                 mating_pool_size,mutation_rate,scoring_args, max_score, prune_population]
        args_list = []
        for i in range(n_tries):
            args_list.append(args + [seeds[count]])
            count += 1

        with Pool(n_cpus) as pool:
            output = pool.map(ga.GA, args_list)

        for i in range(n_tries):     
            #(scores, population) = ga.GA([population_size, file_name,scoring_function,generations,mating_pool_size,mutation_rate,scoring_args,max_score,prune_population])
            (scores, population, high_scores, generation) = output[i]
            high_scores_list.append(high_scores)
            print(f'{i},{scores[0]:.2f},{population[0]},{generation},{s_type},{prune_population}')

pickle.dump(high_scores_list, open('troglitazone_rediscovery.p ', 'wb' ))
