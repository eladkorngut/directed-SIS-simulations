# directed-SIS-simulations
The following files contain numerical methods for simulating SIS dynamics on directed networks.

1. The file qsd_bimodal_probability_matrix.m, contains a method of finding the extinction time and the quasi-stationary distribution by solving
the master equations. The master equations are solved using the matlab only for bimodal networks with different rates

2. shooting.py solves the Hamilton equations using python ODE solver. The function is limited to bimodal Hamilton equations. The Data folder contains examples 
of paths found using this approach.

3. The files gillespierunhomo.py, nethomo.py, netinithomo.py rand_networks.py and slurm.serjob, are dedicated to runing Monte Carlo simulations. 
The file rand_networks.py creates directed networks with heterogenous degree according to the configuration model. The file netinithomo.py creates directed
graphs attributes and rates that can be heterogenous. gillespierunhomo.py the SIS dyamics on such networks while nethomo.py recieve parameters that determine 
which networks was choosen.
The slurm.serjob uses the slurm program to subumit jobs to a cluster.
