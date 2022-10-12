# This program creates graphs and submit the work to a cluster using the slurm file
# The graph created are either bimodal gauss or gamma generated using a modified networkx configuration model graph
# generator, so nodes don't have loops or double edges

import os
import networkx as nx
import numpy as np
import netinithomo
import rand_networks
import csv

if __name__ == '__main__':
    Epsilon_sus = [0.1]
    Epsilon_inf = [-0.2]
    epsilon = 0.1
    eps_din,eps_dout = 0.0,0.0
    eps_sus,eps_lam = 0.3,-0.3
    N = 2000
    k = 200
    x = 0.2
    Num_inf = int(x * N)
    Alpha = 1.0
    susceptibility = 'bimodal'
    infectability = 'bimodal'
    directed_model='uniform_c'
    prog = 'br' #can be either 'i' for the inatilization and reaching eq state or 'r' for running and recording fluc
    Lam = 0.2
    Time_limit = 1100
    Start_recording_time = 100.0
    Beta_avg = Alpha*Lam / k
    Num_different_networks= 1
    Num_inital_conditions= 1
    bank = 1000000
    parts = 1
    foldername ='syncro_N2000_k200_alpha10_startime100_timelimit1100_dt001_eps0_lam02'
    graphname  = 'GNull'
    count = 0
    susceptibility_avg = 1.0
    infectability_avg = 1.0
    sus_inf_correlation = 'c'

    if  prog=='ec' or prog=='ac' or prog=='bd' or prog=='co' or prog=='cr' or prog=='br':
        os.mkdir(foldername)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(foldername)
    if prog =='ec':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1+epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_correlated(N,epsilon_inf,epsilon_sus,1.0)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='ac':
        for epsilon_sus,epsilon_inf in zip(Epsilon_sus,Epsilon_inf):
            Beta=Beta_avg/(1-epsilon_sus*epsilon_inf)
            for n in range(Num_different_networks):
                beta_inf,beta_sus=netinithomo.bi_beta_anti_correlated(N,epsilon_inf,epsilon_sus,1.0)
                # G = nx.random_regular_graph(k, N)
                G = nx.complete_graph(N)
                G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
                infile = graphname + '_' + str(epsilon_sus).replace('.', '') + '_' + str(n)+'.pickle'
                nx.write_gpickle(G, infile)
                outfile ='o'+str(epsilon_sus).replace('.', '')
                for p in range(parts):
                    os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                              str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='bd':
        Beta=Beta_avg/(1+eps_din*eps_dout)
        d1_in, d1_out, d2_in, d2_out =int(k*(1-eps_din)),int(k*(1-eps_dout)),int(k*(1+eps_din)),int(k*(1+eps_dout))
        for n in range(Num_different_networks):
            G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out,d2_in,d2_out,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_d1in' + str(d1_in).replace('.', '') +'_o_d1out' + str(d1_out).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='co':
        for n in range(Num_different_networks):
            G = rand_networks.configuration_model_directed_graph(directed_model, eps_din,eps_dout,k,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_eps_in' + str(np.abs(eps_din)).replace('.', '') +'eps_dout' + str(np.abs(eps_dout)).replace('.', '')
            k_avg_graph = np.mean([G.in_degree(n) for n in G.nodes()])
            Beta_graph = Lam/k_avg_graph
            eps_in_graph = np.std([G.in_degree(n) for n in G.nodes()])/k_avg_graph
            eps_out_graph = np.std([G.out_degree(n) for n in G.nodes()])/k_avg_graph
            Beta = Beta_graph / (1 + np.sign(eps_din)*eps_in_graph * np.sign(eps_dout)* eps_out_graph)
            f = open('parameters_'+outfile + '.csv', "a+")
            with f:
                writer = csv.writer(f)
                writer.writerows([[k_avg_graph, np.sign(eps_din)*eps_in_graph,np.sign(eps_dout)*eps_out_graph]])
            f.close()
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='cr':
        for n in range(Num_different_networks):
            beta_inf,beta_sus=netinithomo.general_beta(N,eps_lam,eps_sus,directed_model,k)
            G = nx.random_regular_graph(k, N)
            # G = nx.complete_graph(N)
            G = netinithomo.intalize_lam_graph(G, N, beta_sus,beta_inf)
            infile = graphname + '_' + str(eps_sus).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            eps_sus_graph = np.std(beta_sus)/np.mean(beta_sus)
            eps_lam_graph = np.std(beta_inf)/np.mean(beta_inf)
            Beta = Beta_avg / (1 + np.sign(eps_sus)*eps_sus_graph * np.sign(eps_lam)* eps_lam_graph)
            outfile ='o'+str(eps_sus).replace('.', '')
            f = open('parameters_'+outfile + '.csv', "a+")
            with f:
                writer = csv.writer(f)
                writer.writerows([[np.mean(beta_sus), np.sign(eps_sus)*eps_sus_graph,np.sign(eps_lam)*eps_lam_graph]])
            f.close()
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' + str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta))
    elif prog =='br':
        # create a bimodal graph with different degrees that record the number of infected nodes up to timelimit
        Beta=Beta_avg/(1+eps_din*eps_dout)
        d1_in, d1_out, d2_in, d2_out =int(k*(1-eps_din)),int(k*(1-eps_dout)),int(k*(1+eps_din)),int(k*(1+eps_dout))
        for n in range(Num_different_networks):
            G = rand_networks.random_bimodal_directed_graph(d1_in, d1_out,d2_in,d2_out,N)
            G = netinithomo.set_graph_attriubute_DiGraph(G)
            infile = graphname + '_' + str(eps_din).replace('.', '') + '_' + str(n)+'.pickle'
            nx.write_gpickle(G, infile)
            outfile ='o_d1in' + str(d1_in).replace('.', '') +'_o_d1out' + str(d1_out).replace('.', '')
            for p in range(parts):
                os.system(dir_path + '/slurm.serjob python3 ' + dir_path + '/gillespierunhomo.py '+str(prog) + ' ' +
                          str(Alpha) + ' ' + str(bank) + ' ' + str(outfile)+ ' ' + str(infile) + ' ' +
                          str(Num_inital_conditions)+ ' ' + str(Num_inf) + ' ' +str(n)+ ' ' +str(Beta)+ ' ' +str(Start_recording_time)+ ' ' +str(Time_limit))

