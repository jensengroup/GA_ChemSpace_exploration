# How to use the input files

```
python GB_GA.py celecoxib_rediscovery_graph.py Celecoxib_guacamol.smi > celecoxib_rediscovery_graph.csv
```

GB_GA.py can be found [here](https://github.com/jensengroup/GB-GA), Celecoxib_guacamol.smi is in the initial_populations directory, and celecoxib_rediscovery_graph.csv is in the data directory. This will also create elecoxib_rediscovery_graph.p. To exactly reproduce the results set

```
seeds = [292273,265159,...,894088,306197]
```
where the seeds are taken from celecoxib_rediscovery_graph.csv
