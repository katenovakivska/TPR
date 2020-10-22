import numpy as np
from collections import defaultdict

with open ("conditions.txt", 'r') as data: #read data 
    dataset = [[int(x) for x in line.split()] for line in data]

make_R = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
R = make_R(dataset, 15)
print("Number of R: {}".format(len(R)))

class Graph:
    def __init__(self, vertices):
        self.graph = defaultdict(list)
        self.V = vertices

    def addEdge(self, u, v):
        self.graph[u].append(v)

    def isCyclicUtil(self, v, visited, recStack):
        # Mark current node as visited and
        # adds to recursion stack
        visited[v] = True
        recStack[v] = True

        # Recur for all neighbours
        # if any neighbour is visited and in
        # recStack then graph is cyclic
        for neighbour in self.graph[v]:
            if visited[neighbour] == False:
                if self.isCyclicUtil(neighbour, visited, recStack) == True:
                    return True
            elif recStack[neighbour] == True:
                return True

        # The node needs to be poped from
        # recursion stack before function ends
        recStack[v] = False
        return False

    # Returns true if graph is cyclic else false
    def isCyclic(self):
        visited = [False] * self.V
        recStack = [False] * self.V
        for node in range(self.V):
            if visited[node] == False:
                if self.isCyclicUtil(node, visited, recStack) == True:
                    return True
        return False


def createGraph(matrix):
    g = Graph(15)
    index = 0
    for row in matrix:
        index1 = 0
        for element in row:
            if(element == 1):
                g.addEdge(index, index1)
            index1 += 1
        index += 1
    return g

def check_x(S):
    max_k = np.array([])
    opt_k = np.array([])
    for i, x in zip(range(len(S)),S):
        if np.any(S-x > 0): #формуємо множину k-макс. елементів
            continue
        max_k = np.append(max_k, i+1)
        if np.sum(x) == 15: #формуємо множину k-опт. елементів
            opt_k = np.append(opt_k, i+1)
        
    return max_k, opt_k

def k_optimization(Rn):
    I = (Rn == Rn.T)*Rn #формуємо симетричну частину
    P = Rn-I #формуємо асиметричну частину
    N = (Rn == Rn.T)-I #формуємо відношення непорівнюваності
    
    for i in range(1,5):
        if i==1:
            S1=I+P+N
            max_1, opt_1 = check_x(S1)
        elif i==2:
            S2=P+N
            max_2, opt_2 = check_x(S2)
        elif i==3:
            S3=P+I
            max_3, opt_3 = check_x(S3)
        elif i==4:
            S4=P
            max_4, opt_4 = check_x(S4)
            
    parameters = {"1_max": max_1,
                 "1_opt": opt_1,
                 "2_max": max_2,
                 "2_opt": opt_2,
                 "3_max": max_3,
                 "3_opt": opt_3,
                 "4_max": max_4,
                 "4_opt": opt_4}
    
    return parameters

def print_information_k_opt(parameters, num_r): 
    print("____R{}____".format(num_r))
    print("1-max: {}".format(parameters["1_max"]))
    print("1-opt: {}".format(parameters["1_opt"]))
    print("2-max: {}".format(parameters["2_max"]))
    print("2-opt: {}".format(parameters["2_opt"]))
    print("3-max: {}".format(parameters["3_max"]))
    print("3-opt: {}".format(parameters["3_opt"]))
    print("4-max: {}".format(parameters["4_opt"]))
    print("4-opt: {}".format(parameters["4_opt"]))

def neumann_morgenshtern_optimization(Rn):
    S = np.array([])
    currentS = np.array([])
    for i, x in zip(range(len(Rn)),Rn.T):
        hasEnters = False
        for j, y in zip(range(len(x)),x):
            if(y != 0):
                hasEnters = True
        if(hasEnters == False):
            currentS = np.append(currentS, i + 1)
    S = np.append(S, [currentS])

    outS =np.array([])
    for i, x in zip(range(len(Rn)),Rn.T):
        isOut = False
        for j, y in zip(range(len(x)),x):
            if(y == 1 and any(el == j + 1 for el in currentS)):
                isOut = True
            elif(y == 1 and all(el != j + 1 for el in currentS)):
                isOut = False
                break
        if(isOut == True):
            outS = np.append(outS, i + 1)
            S = np.append(S, [outS])
    return S, outS


for num_r, r_i in zip(range(len(R)),R):
    Rn = np.array(r_i)
    graph = createGraph(Rn)
    if graph.isCyclic() == 1: 
        print ("Graph has a cycle")
        parameters = k_optimization(Rn)
        print_information_k_opt(parameters, num_r)
    else: 
        print ("Graph has no cycle")
        s, outs = neumann_morgenshtern_optimization(Rn)
        print(s)
        print(outs)