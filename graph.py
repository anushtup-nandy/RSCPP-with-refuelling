import numpy as np
import random
import sys
#List containing the nodes:
'''
    each element will be 
    another list containing: 
            [node number, price, x-coordinate, y-coordinate]
'''
# graph=[
#     [0,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [1,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [2,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [3,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [4,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [5,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [6,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [7,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [8,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [9,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [10,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [11,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [12,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [13,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [14,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [15,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [16,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [17,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [18,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [19,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     [20,random.uniform(1,3), random.uniform(200,1000), random.uniform(200,1000)],
#     ['start',0, random.uniform(200,1000), random.uniform(200,1000)],
#     ['finish', np.nan, random.uniform(200,1000), random.uniform(200,1000)]
# ]

# graph=[
#     [0, 1.46, 400, 500],
#     [1, 1.56, 450, 600],
#     [2, 2.67, 340, 550],
#     [3, 2.2, 800, 250],
#     [4, 1.2, 740, 650],
#     ['start', 0,  950, 800],
#     ['finish', np.nan, 210, 350]
# ]

graph=[
    [0,1,0,0],
    [1,2,3,2],
    [2,3,3,3],
    [3,4,5,4],
    ['start', 0, ]
]

#The pre-requesite details:
# Kilo_per_gas = 10.0
# START_GAS = 20
# MAX_GAS = 40
# START_RANGE = START_GAS * Kilo_per_gas
# MAX_RANGE = MAX_GAS * Kilo_per_gas
# MAX_STOP = 10 # limit on how many times driver can stop

Kilo_per_gas = 1.0
START_GAS = 2
MAX_GAS = 4
START_RANGE = START_GAS * Kilo_per_gas
MAX_RANGE = MAX_GAS * Kilo_per_gas
MAX_STOP = 2 # limit on how many times driver can stop


#distance:
def distance(x_list,y_list):
    dist=np.zeros((len(x_list), len(y_list)))
    for i in range(len(x_list)):
        for j in range(len(x_list)):
            dist[i][j]=(np.sqrt((x_list[i]-x_list[j])**2 + (y_list[i]-y_list[j])**2))
    return dist


#loading the data from the table of nodes:
def load(graph):
    num_points=len(graph)-1
    x_list=[graph[i][2] for i in range(len(graph))]
    y_list=[graph[i][3] for i in range(len(graph))]
    z=[graph[i][1] for i in range(len(graph)-1)]
    price=[i/Kilo_per_gas for i in z]
    idx=[graph[i][0] for i in range(len(graph))]
    return (x_list, y_list, price, num_points, idx)

#discretizing since it has been asked in 2.1
'''
    Generates the set that contains the possible values of u
'''
def set_of_u(dist, price):
    # UV={ MAX_RANGE - d[u,v] | for all v where price[v] < price[u]
    #                            and d[u,v] <= MAX_RANGE}
    range_values = {}
    C = MAX_RANGE * np.ones((1, num_points))
    C[0,-1] = START_RANGE
    temp = C - dist
    for u in range(num_points - 1):
        mask = np.logical_and(np.array(price) < price[u],temp[u,:] >= 0) 
        range_values[u] = np.append(temp[u,mask],0)
        range_values[u].sort()
    return range_values

#initializing the cost
'''
    This will initialize the cost table
    Gives the output of a dictionary from (u, range values)-->(minimum cost, path taken)

        C[u,rv]={ (d[u,t]-rv)*price[u] if rv<= d[u,t]<=maximum range the car can travel}
'''
def init_cost_dict(dist, range_val, price, num_points):
    cost={}
    for u in range(num_points-1):
        for rv in range_val[u]:
            if (rv<=distance_from_target[u] and distance_from_target[u]<=MAX_RANGE):
                cost[(u, rv)]=((distance_from_target[u]-rv)*price[u], [num_points])
            else:
                cost[(u, rv)]=(np.inf, [])

    return cost


#cost mapping for the starting node! (earlier function finds the cost mapping for 0-20 nodes)
def cost_final(old_cost, dist, price, range_val, num_points):
    cost={}
    for u in range(num_points-1):
        indep1v=[] #corresponds to C_old[v,0] + d[u,v] * (price[u])
        indep2v=[] #corresponds to C_old[v,MAX_RANGE - d[u,v]] + MAX_RANGE * price[u] 

        for v in range(num_points-1):
            d=dist[u,v]
            if (d> MAX_RANGE or u==v):
                continue

            if (price[v]<=price[u]):
                c, prev=old_cost[(v,0)]

                if (prev):
                    prev=list(prev) #copying prev
                    prev.append(v) #adding the next node
                    indep1v.append((c + d*(price[u]), prev))
            else:
                c, prev=old_cost[(v, MAX_RANGE-d)]
                if (prev):
                    prev=list(prev)
                    prev.append(v)
                    indep2v.append((c + MAX_RANGE*price[u], prev))

        #finding the minimum of {indep1v, indep2v} 
        indep1v.sort(reverse=True)
        if indep2v:
            min_indep2v=min(indep2v)
        else:
            min_indep2v=(np.inf, [])
        
        #boolean that returns True if minimum is from indep2v
        if indep1v:
            from_indep2v= min_indep2v <= indep1v[-1]
        else:
            from_indep2v=True

        # fill in the range dependent values and deactivate terms in
        # ind1 list once rv > d[u,v]
        for rv in range_val[u]:
            if (from_indep2v):
                cost[(u,rv)]=(min_indep2v[0] - rv*price[u], min_indep2v)
            else:
                while(indep1v and dist[indep1v[-1][1][-1], u] < rv):
                    del indep1v[-1]

                if (indep1v and min_indep2v > indep1v[-1]):
                    cost[(u,rv)]=(indep1v[-1][0] - rv*price[u], indep1v[-1][1])
                else:
                    from_indep2v=True
                    cost[(u,rv)]=(min_indep2v[0] - rv*price[u], min_indep2v[1])
    return cost

# this will return the cost it takes to go from start point s to final point t
def start_finish(cost, dist, num_points):
    min_cost=np.inf
    prev=[]
    for v in range(num_points-1):
        d=dist[v,-1] #from start
        if (d > START_RANGE):
            continue
        c_cost, prev_v=cost[(v, START_RANGE-d)]
        if (c_cost < min_cost):
            min_cost=c_cost
            prev= list(prev_v) #copying
            prev.append(v)
            prev.pop(0) #popping the finishing index
    return (min_cost, prev)


#ouput starts here:
x_list, y_list, price, num_points, idx=load(graph)
dist=distance(x_list, y_list)
distance_from_target=dist[-1,:-1]
dist=dist[:-1,:-1]
range_values=set_of_u(dist, price)
old_cost={}
cost=init_cost_dict(dist, range_values,price,num_points)


if (distance_from_target[num_points - 1] <= START_RANGE):
    print ('Problem not well defined, can reach target from start '
           ' with initial range')
    quit()

for n in range(1, MAX_STOP + 1):
    total_cost, prev=start_finish(cost, dist, num_points)
    if (total_cost !=np.inf):
        print("path found with cost :", total_cost, " ! ")
    old_cost.clear()
    old_cost=cost.copy()
    
    cost=cost_final(old_cost, dist, price, range_values, num_points)
    sys.stdout.write('Computing cost for {} stop'.format(n))
    sys.stdout.write('\n')

