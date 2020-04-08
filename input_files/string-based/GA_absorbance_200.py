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

file_name = sys.argv[1]

scoring_function = sc.absorbance_target
target = 200. # nm
sigma = 50. # nm
threshold = 0.3
n_confs = 20
xtb_path = '/home/jhjensen/stda'
scoring_args = [n_confs, xtb_path, target, sigma, threshold]
max_score = 1.99

population_size = 20 
mating_pool_size = 20
generations = 50
mutation_rate = 0.05
co.average_size = 50. 
co.size_stdev = 5.
n_tries = 10
n_cpus = 8

print('population_size', population_size)
print('mating_pool_size', mating_pool_size)
print('generations', generations)
print('mutation_rate', mutation_rate)
print('max_score', max_score)
print('average_size/size_stdev', co.average_size, co.size_stdev)
print('initial pool', file_name)
print('target +/ sigma', target, sigma)
print('number of tries', n_tries)
print('number of CPUs', n_cpus)
print('')

line = []
for s_type in ['SMILES', 'DeepSMILES', 'SELFIES'][2:3]: #NB!
    co.string_type = s_type
    line_list = [f'{co.string_type}']
    for prune_population in [True, False][1:2]: #NB!
        print('string type', co.string_type)
        print('prune_population', prune_population)
        print('')

        results = []
        size = []
        t0 = time.time()
        all_scores = []
        generations_list = []
        args = n_tries*[[population_size, file_name,scoring_function,generations,
            mating_pool_size,mutation_rate,scoring_args, max_score, prune_population]]
        with Pool(n_cpus) as pool:
            output = pool.map(ga.GA, args)

        for i in range(n_tries):     
            #(scores, population) = ga.GA([population_size, file_name,scoring_function,generations,mating_pool_size,mutation_rate,scoring_args,max_score,prune_population])
            (scores, population, generation) = output[i]
            all_scores.append(scores)
            print(f'{i} {scores[0]:.2f} {co.string2smiles(population[0])} {generation}')
            results.append(scores[0])
            generations_list.append(generation)

        t1 = time.time()
        print('')
        print(f'max score {max(results):.2f}, mean {np.array(results).mean():.2f} +/- {np.array(results).std():.2f}')
        print(f'mean generations {np.array(generations_list).mean():.2f} +/- {np.array(generations_list).std():.2f}')
        print(f'time {(t1-t0)/60.0:.2f} minutes')
        print('')

        line_list.append(f'{max(generations_list)} & {np.array(generations_list).mean():.1f}$\pm${np.array(generations_list).std():.1f}')


    line.append(' & '.join(line_list)+' \\'+'\\')

print('')
for item in line:
    print(item)