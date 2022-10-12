# This program run a SIS dynamical process using the gillespie algorithm on different tyeps of networks. As a parameter
# the program receives a on which the SIS infection occurs


import numpy as np
import bisect
import netinithomo
import networkx as nx
import csv
import sys


def fluctuation_run_extinction(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    # Run a SIS process until extinction is reached
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_graph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.noes[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1


            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot - Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                    else:
                        Rates[person] = Rates[person] + Beta*(G.nodes[person]['contact_rate'] * G.nodes[Neighbor]['spread_rate'])
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G[person]:
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                        R_tot = R_tot + Beta*(G.nodes[Neighbor]['contact_rate'] * G.nodes[person]['spread_rate'])
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0


def fluctuation_run_extinction_DiGraph(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta):
    # Run a SIS process on a directed graph until extinction i reached and record the extinction time
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        Total_time = 0.0
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            else:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
        f = open(outfile+'.csv',"a+")
        with f:
            writer = csv.writer(f)
            writer.writerows([[Total_time,network_number]])
        f.close()
    return 0



def fluctuation_run_extinction_DiGraph_record(Alpha,bank,outfile,infile,runs,Num_inf,network_number,Beta,start_recording_time,Time_limit):
    # Run a SIS process on a directed graph until a timelimit is reached and record the number of infected at different time intrevals
    G = nx.read_gpickle(infile)
    seed_nodes = Num_inf
    for run_loop_counter in range(runs):
        T,I,I_type_above_avg,I_type_below_avg,runs_csv=[],[],[],[],[]
        runs_csv.append(run_loop_counter)
        Total_time = 0.0
        T.append(Total_time)
        count = 0
        Num_inf = seed_nodes
        r = np.random.uniform(0, 1, (bank, 2))
        for l in range(G.number_of_nodes()):
            G.nodes[l]['infected'] = False
        R_tot, Rates = netinithomo.inatlize_inf_DiGraph(G,Num_inf,G.number_of_nodes(),Alpha,Beta)
        net_num = []
        I.append(Num_inf)
        net_num.append(network_number)

        num_inf_above_avg, num_inf_below_avg = 0, 0
        for node_type_number in range(G.number_of_nodes()):
            if G.nodes[node_type_number]['infected'] == True:
                if G.nodes[node_type_number]['contact_rate'] > 1.0:
                    num_inf_above_avg = num_inf_above_avg + 1
                else:
                    num_inf_below_avg = num_inf_below_avg + 1
        I_type_above_avg.append(num_inf_above_avg)
        I_type_below_avg.append(num_inf_below_avg)


        ######################
        # Main Gillespie Loop
        ######################
        while Num_inf > 0 and Total_time<Time_limit:
            R_norm = np.cumsum(Rates)
            r_pos = R_tot * r[count, 1]
            person = bisect.bisect_left(R_norm, r_pos)
            tau= np.log(1 / r[count, 0]) / R_tot
            Total_time = Total_time + tau

            try:
                if G.nodes[person]['infected'] == True:
                  pass
            except:
                  print('Accessing G.node[person][infected] failed value of person is ',person)
                  if person == G.number_of_nodes():
                      person =G.number_of_nodes()-1

            if G.nodes[person]['infected'] == True and Num_inf>1:
                Num_inf = Num_inf - 1
                Rates[person] = 0.0
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] - Beta
                        R_tot = R_tot - Beta
                for Neighbor in G.predecessors(person):
                    if G.nodes[Neighbor]['infected'] == True:
                        Rates[person] = Rates[person] + Beta
                R_tot = R_tot + Rates[person] - Alpha
                G.nodes[person]['infected'] = False
            elif G.nodes[person]['infected'] == False:
                Num_inf = Num_inf + 1
                for Neighbor in G.successors(person):
                    if G.nodes[Neighbor]['infected'] == False:
                        Rates[Neighbor] = Rates[Neighbor] + Beta
                        R_tot = R_tot + Beta
                R_tot = R_tot - Rates[person] + Alpha
                Rates[person] = Alpha
                G.nodes[person]['infected'] = True
            count = count + 1
            if count >= bank:
                r = np.random.uniform(0, 1, (bank, 2))
                count = 0
            if Total_time - T[-1] >= 0.1 and Total_time >= start_recording_time:
                I.append(Num_inf)
                T.append(Total_time)
                net_num.append(network_number)
                runs_csv.append(run_loop_counter)
                num_inf_above_avg, num_inf_below_avg = 0, 0
                for node_type_number in range(G.number_of_nodes()):
                    if G.nodes[node_type_number]['infected'] == True:
                        if G.nodes[node_type_number]['contact_rate'] > 1.0:
                            num_inf_above_avg = num_inf_above_avg + 1
                        else:
                            num_inf_below_avg = num_inf_below_avg + 1
                I_type_above_avg.append(num_inf_above_avg)
                I_type_below_avg.append(num_inf_below_avg)
        f = open(outfile + '.csv', "a+")
        l = [T, I, I_type_above_avg, I_type_below_avg, runs_csv, net_num]
        l = zip(*l)
        with f:
            writer = csv.writer(f)
            writer.writerows(l)
        f.close()
    return 0


def actasmain():
    # The function allows to run the simulations on a PC instead of submitting the jobs to a cluster, one can insert
    # code using the above function and get numerical results
    Epsilon_sus = [0.0]
    Epsilon_inf = [0.0]
    Epsilon=[0.0]
    N = 500
    k = 200
    x = 0.2
    eps_din,eps_dout = 0.0,0.0
    eps_sus,eps_lam = 0.0,0.0
    Num_inf = int(x * N)
    Alpha = 1.0
    susceptibility = 'bimodal'
    infectability = 'bimodal'
    directed_model='gauss_c'
    prog = 'q' #can be either 'i' for the inatilization and reaching eq state or 'r' for running and recording fluc
    Lam = 1.5
    Time_limit = 200
    Start_recording_time = 0.0
    Beta_avg =Alpha* Lam / k
    Num_different_networks= 1
    Num_inital_conditions= 1
    bank = 1000000
    parts = 1
    graphname  = 'GNull'
    count = 0
    susceptibility_avg = 1.0
    infectability_avg = 1.0
    foldername ='base'
    graphname  = 'GNull'
    outfile ='o'
    sus_inf_correlation = 'a'


if __name__ == '__main__':
    # if submit equals True the program will submit the jobs to the cluster using the slurm file otherwise the program
    # run an infection process on the PC
    submit = True
    if submit==False:
        actasmain()
    else:
         if sys.argv[1]=='ec' or sys.argv[1]=='ac' or sys.argv[1] == 'cr':
             # Run and record extinction time for heterogeneous rate networks
             fluctuation_run_extinction(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))
         elif sys.argv[1] == 'bd' or sys.argv[1] == 'co':
             # Run and record extinction events for heterogeneous degree networks
             fluctuation_run_extinction_DiGraph(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),float(sys.argv[9]))

         elif sys.argv[1] == 'br':
            # Run and record number of infected for bimodal degree graph
             fluctuation_run_extinction_DiGraph_record(float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5],
                                                int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]),
                                                float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11]))


