import numpy as np
import cvxpy as cp
import random
import csv
from datetime import datetime
import time
import sys

np.random.seed(29)

def recurse_son(father, Spanning_tree, ratio_now, visited, gen_line_prop_mat):
    node_current = father
    visited[node_current] = True
    child_bus, child_connected_line = find_child(A_mat_directed, node_current)
    father_ratio=np.copy(ratio_now)
    #print("Current father", node_current)
    #print("Current son", child_bus)
    #print("Child connected line", child_connected_line)

    if child_bus==[]:
        return

    for i in range(len(child_bus)):
        current_son = child_bus[i]
        if visited[current_son]:
            #print("Already visited, working on son ", current_son)
            #print("Current father ", node_current)
            #print("Current ratio", father_ratio)
            ratio_now = line_prop_mat[child_connected_line[i], node_current] * father_ratio
            gen_line_prop_mat[child_connected_line[i], current_son] = ratio_now
            Spanning_tree[current_son] += bus_prop_vec[current_son] * ratio_now
            #print("New Ratio", ratio_now)
            recurse_son(current_son, Spanning_tree, ratio_now, visited, gen_line_prop_mat)

        if not visited[current_son]:
            #print("Working on son ", current_son)
            #print("Previous ratio", father_ratio)
            ratio_now = line_prop_mat[child_connected_line[i], node_current] * father_ratio
            gen_line_prop_mat[child_connected_line[i], current_son] = ratio_now
            Spanning_tree[current_son] += bus_prop_vec[current_son] * ratio_now
            #print("New Ratio", ratio_now)
            recurse_son(current_son, Spanning_tree, ratio_now, visited, gen_line_prop_mat)
    return



def find_child(directed_graph, current_vertex):
    son=[]
    line_index=[]
    for i in range(num_lines):
        if directed_graph[i][current_vertex]==-1:
            index=np.argwhere(directed_graph[i,:]==1)
            son.append(index)
            line_index.append(i)
    son = np.array(son).reshape(-1, 1)
    line_index=np.array(line_index).reshape(-1,1)
    return son, line_index


np.set_printoptions(threshold=sys.maxsize)
num_gen=6
num_buses=30
num_lines=41

with open('30bus_topology.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    rows = [row for row in reader]
    connection_all = np.array(rows[1:], dtype=float)

A_mat = np.zeros((num_lines, num_buses), dtype=float)
for i in range(num_lines):
    A_mat[i][int(connection_all[i][0])-1] = -1.0
    A_mat[i][int(connection_all[i][1])-1] = 1.0

carbon_emissions=np.array([3.0, 2.7, 2.0, 1.5, 1.0, 1.8])
Gen_node=np.array([0, 1, 12, 21, 22, 26])
C = np.array([1.0, 1.5, 2.4, 3.5, 5.0, 3.0])
D1= np.diag(np.full(6, 1))
D=np.stack((D1, -D1), axis=0)
D=np.reshape(D, (12, 6))
e = np.array([10.0, 15.0, 15.0, 15.0, 15.0, 15.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0])

b = np.zeros((30, 6), dtype=float)
b[0][0]=1.0
b[1][1]=1.0
b[12][2]=1.0
b[21][3]=1.0
b[22][4]=1.0
b[26][5]=1.0

x = cp.Variable(num_gen)
line_flow = cp.Variable(num_lines)
load = np.random.uniform(0.2, 1.8, (num_buses, 1))  # Uniform distribution now, change the scale here
load = np.round(load, 2)
line_flow_limit = np.ones((1, num_lines), dtype=float) * 4.0
neg_line_flow_limit = -np.ones((1, num_lines), dtype=float) * 4.0

cost = cp.sum(C @ x)
constraint = [b @ x + A_mat.T @ line_flow == load.reshape(-1, ),
              D @ x <= e,  # Generation constraint
              line_flow <= line_flow_limit.reshape(-1, ),
              line_flow >= neg_line_flow_limit.reshape(-1,)]  # line flow limits
prob = cp.Problem(cp.Minimize(cost), constraint)
prob.solve(solver=cp.CVXOPT)

print("Load, ", load.T)
print("Generation value", np.round(x.value,2))
print("Line flow value", np.round(line_flow.value,3))

A_mat_directed=np.copy(A_mat)

for i in range(num_lines):
    if line_flow.value[i]<0:
        A_mat_directed[i][int(connection_all[i][0]) - 1] = 1.0
        A_mat_directed[i][int(connection_all[i][1]) - 1] = -1.0


total_injection=np.dot(b, x.value)
print("Total injection", np.round(total_injection,2))
line_prop_mat=np.zeros((num_lines, num_buses), dtype=float)
bus_prop_vec=np.zeros((num_buses, 1), dtype=float)

for i in range(num_buses):
    total_power = np.copy(load[i])
    total_inflow = 0.0
    for j in range(num_lines):
        if A_mat_directed[j][i] == -1:
            total_power += np.abs(line_flow.value[j])
        elif A_mat_directed[j][i] == 1:
            total_inflow += np.abs(line_flow.value[j])

    #print("Node ", i)
    #print("Total power outflow and demand:", total_power)
    #print("Total generator injections at this node:", total_injection[i])
    #print("Total line flow injections at this node:", total_inflow)
    bus_prop_vec[i] = load[i] / total_power
    for j in range(num_lines):
        if A_mat_directed[j][i] == -1:
            line_prop_mat[j, i] = np.abs(line_flow.value[j]) / total_power

#print("FINISH!!!!!!!!!")
#print("Directed incidence matrix", A_mat_directed, flush=True)
#print("Line proportion matrix", line_prop_mat)
#print("Bus proportion vector", bus_prop_vec.T)

Gen_prop_mat=np.zeros((num_gen, num_lines, num_buses), dtype=float)
Gen_spanning_tree_all=np.zeros((num_buses, num_gen))

for i in range(num_gen):
    print("Generator", Gen_node[i])
    print("Generation Value", x.value[i])
    Gen_spanning_tree = np.zeros((num_buses, 1), dtype=float)
    current_Gen_prop_mat = np.zeros((num_lines, num_buses), dtype=float)
    visited = np.zeros(num_buses, dtype=bool)
    child = True
    current_ratio = 1.0
    current_node = Gen_node[i]
    Gen_spanning_tree[current_node] = bus_prop_vec[current_node] * current_ratio
    recurse_son(current_node, Gen_spanning_tree, current_ratio, visited, current_Gen_prop_mat)
    print("Generator's contribution to each load", np.round(Gen_spanning_tree, 3).T)
    print("Sum of generator's output", np.sum(Gen_spanning_tree))
    #print("Generator's contribution to each line:", current_Gen_prop_mat)
    Gen_prop_mat[i]=current_Gen_prop_mat
    Gen_spanning_tree_all[:,i]=Gen_spanning_tree.reshape(-1, )

load_vec=np.zeros((num_buses, 1), dtype=float)
carbon_vec=np.zeros((num_buses, 1), dtype=float)
for i in range(num_buses):
    for j in range(num_gen):
        load_vec[i] += Gen_spanning_tree_all[i,j] * x.value[j]
        carbon_vec[i] += Gen_spanning_tree_all[i,j] * x.value[j] * carbon_emissions[j]
print("Recovered Load vector", load_vec.T)
print("Original load vector", load.T)
print("Carbon emissions vector", carbon_vec.T)
average_emission_rate = carbon_vec/load_vec
print("Carbon emissions rate vector", average_emission_rate.T)






x = cp.Variable(num_gen)
line_flow = cp.Variable(num_lines)
load_new = np.copy(load)  # Uniform distribution now, change the scale here
load_new = np.round(load_new, 2)
load_new[0] += 0.1
line_flow_limit = np.ones((1, num_lines), dtype=float) * 4.0
neg_line_flow_limit = -np.ones((1, num_lines), dtype=float) * 4.0

cost = cp.sum(C @ x)
constraint = [b @ x + A_mat.T @ line_flow == load_new.reshape(-1, ),
              D @ x <= e,  # Generation constraint
              line_flow <= line_flow_limit.reshape(-1, ),
              line_flow >= neg_line_flow_limit.reshape(-1,)]  # line flow limits
prob2 = cp.Problem(cp.Minimize(cost), constraint)
prob2.solve(solver=cp.CVXOPT)

print("Load, ", load_new.T)
print("Generation value", np.round(x.value,2))
print("Line flow value", np.round(line_flow.value,3))

A_mat_directed=np.copy(A_mat)

for i in range(num_lines):
    if line_flow.value[i]<0:
        A_mat_directed[i][int(connection_all[i][0]) - 1] = 1.0
        A_mat_directed[i][int(connection_all[i][1]) - 1] = -1.0


total_injection=np.dot(b, x.value)
print("Total injection", np.round(total_injection, 2))
line_prop_mat=np.zeros((num_lines, num_buses), dtype=float)
bus_prop_vec=np.zeros((num_buses, 1), dtype=float)

for i in range(num_buses):
    total_power = np.copy(load_new[i])
    total_inflow = 0.0
    for j in range(num_lines):
        if A_mat_directed[j][i] == -1:
            total_power += np.abs(line_flow.value[j])
        elif A_mat_directed[j][i] == 1:
            total_inflow += np.abs(line_flow.value[j])

    #print("Node ", i)
    #print("Total power outflow and demand:", total_power)
    #print("Total generator injections at this node:", total_injection[i])
    #print("Total line flow injections at this node:", total_inflow)
    bus_prop_vec[i] = load_new[i] / total_power
    for j in range(num_lines):
        if A_mat_directed[j][i] == -1:
            line_prop_mat[j, i] = np.abs(line_flow.value[j]) / total_power

#print("FINISH!!!!!!!!!")
#print("Directed incidence matrix", A_mat_directed, flush=True)
#print("Line proportion matrix", line_prop_mat)
#print("Bus proportion vector", bus_prop_vec.T)

Gen_prop_mat=np.zeros((num_gen, num_lines, num_buses), dtype=float)
Gen_spanning_tree_all=np.zeros((num_buses, num_gen))

for i in range(num_gen):
    print("Generator", Gen_node[i])
    print("Generation Value", x.value[i])
    Gen_spanning_tree = np.zeros((num_buses, 1), dtype=float)
    current_Gen_prop_mat = np.zeros((num_lines, num_buses), dtype=float)
    visited = np.zeros(num_buses, dtype=bool)
    child = True
    current_ratio = 1.0
    current_node = Gen_node[i]
    Gen_spanning_tree[current_node] = bus_prop_vec[current_node] * current_ratio
    recurse_son(current_node, Gen_spanning_tree, current_ratio, visited, current_Gen_prop_mat)
    print("Generator's contribution to each load", np.round(Gen_spanning_tree, 3).T)
    print("Sum of generator's output", np.sum(Gen_spanning_tree))
    #print("Generator's contribution to each line:", current_Gen_prop_mat)
    Gen_prop_mat[i]=current_Gen_prop_mat
    Gen_spanning_tree_all[:,i]=Gen_spanning_tree.reshape(-1, )

load_vec=np.zeros((num_buses, 1), dtype=float)
carbon_vec=np.zeros((num_buses, 1), dtype=float)
for i in range(num_buses):
    for j in range(num_gen):
        load_vec[i] += Gen_spanning_tree_all[i,j] * x.value[j]
        carbon_vec[i] += Gen_spanning_tree_all[i,j] * x.value[j] * carbon_emissions[j]

print("Recovered Load vector", load_vec.T)
print("New load vector", load_new.T)
print("Original load vector", load.T)
print("Carbon emissions vector", carbon_vec.T)
average_emission_rate = carbon_vec/load_vec
print("Carbon emissions rate vector", average_emission_rate.T)














