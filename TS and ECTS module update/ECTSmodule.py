import math
import random
import DSSPmodule as DSSP 
import time
#Enhanced Continuous Tabu Search combined with DSSPmodule
#Inputs: -------------------------------
# m = GUROBI model -- specifically a linear relaxation of the Fixed Cost Network Flow problem
# maxIter = maximum number of iterations for tabu search
# p: parameter to modify DSSP algorithm
# tabulength = length of tabu list
# neighbors = number of neighbors search
# arcs = set of arcs (i,j) in the network model; data type: tuple
# varcost: dictionary of variable costs for each commodity k on arc (i,j); denoted mathematically: c_ijk
# fixedcost: dictionary of fixed costs for each arc (i,j); denoted mathematically: f_ij
# totalSupply: dictionary of total supply of commodity k in the network, denoted mathematically: M_k
# K: the number of commodities in the network (NOTE: up until now all work has K = 1)
# flow: dictionary of GUROBI model variables which represent the flow of commodity k on arc (i,j); denoted mathematically: x_ijk
# hk: radius of the external ball
# h0: radius of the internal ball
# e: radius of the tabu list balls, equal to h0

#Outputs: ------------------------------
# iterations: total iterations completed
# pbest: best p value found
# gbestObj: best objective value found

def Choose_neighbors(p,k, h0, hk, tabulist):
    neighborslist = []
    
    h = []
    for i in xrange(0,k+1):
        h.append(0)    

    #caculate the h[0],h[1],h[2],...,h[k]
    h[0] = h0
    for i in xrange(1,k+1):
        h[k-i+1] = hk/(2**(i-1))
    
    #add neighbors into list between h[x], h[x+1]
    for i in xrange(1,k+1):
        flag = False
        while(flag == False):
            if random.random()<0.5:
                temp = random.uniform(p+h[i-1],p+h[i])   #look for neighbor between p+h[i-1],p+h[i]
            else: temp = random.uniform(p-h[i-1],p-h[i])   #look for neighbor between p-h[i-1],p-h[i]
            
            #compare with pmin and pmax
            if temp > 2:
                temp =  temp - 1
            elif temp <= 0:
                temp = abs(temp)         
            
            for t in xrange(0,len(tabulist)):
                flag = abs(tabulist[t]-temp) > h0 #check the tabustates of neighborlist
            if flag == True:
                break
            
        neighborslist.append(temp)
        
    return neighborslist

def update_tabulist(a,b):
    a.remove(a[0])
    a.append(b)
    return a

def main (m, maxIter, max_time, p, tabulength, neighbors, arcs, varcost, fixedcost, totSupply, K, flow, decision, h0, hk,seed):
    
    random.seed(seed)
    
    #solve the problem with original p
    t0 = time.clock()
    result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, p, decision)
    DSSPtime = time.clock()-t0
    iterations = result[0]
    gbestObj = result[2]
    pbest = p
    curp = p
    proball = []
    proball.append(result[2])
    
    neighborslist = Choose_neighbors(curp,neighbors, h0, hk, [curp])
    for i in range (0,len(neighborslist)):
        result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, neighborslist[i], decision)
        proball.append(result[2])
    
    proBest = sum(proball)  #sum of all points in the promising ball
    procenter = p  #center of the promising area
    proBestObj_track = proBest+1000
    
    #creat promising list
    prolist = []
    prolist.append(curp)
    
    long_term_memory = {curp:0} #creat long-term memory 
    long_term_max = 2           #define the max frequency of long_term memory    
    
    early_termin_cnt = 0   #count of iterations gbestObj is not improved
    gbestObj_track = gbestObj+1000   #initialize the gbestObj track
    
    #define the total iteraions of search
    DSSPIter = min(maxIter,int(max_time/DSSPtime))
    Iter = int(DSSPIter/neighbors)    
    
    #diversified search for promising area
    for i in xrange (0,int(0.3*Iter)):
        #early termination condition check
        if proBestObj_track == proBest:
            early_termin_cnt += 1
            if early_termin_cnt >= 3:
                break    #stop loop
        else: early_termin_cnt = 0
                        
        proBestObj_track = proBest                
        
        proball = []
        
        flag = False
        while(flag == False):
            curp = random.uniform(0,2)      
            for t in xrange(0,len(prolist)):
                flag = abs(prolist[t]-curp) > hk #check if the random point is in any existed promising ball
            if flag == True:
                break
        
        prolist.append(curp)
        m.reset()
        result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, curp, decision)
        iterations += result[0]
        proball.append(result[2])
        neighborslist = Choose_neighbors(curp,neighbors, h0, hk, prolist)
        for i in range (0,len(neighborslist)):
            m.reset()
            result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, neighborslist[i], decision)
            proball.append(result[2])            
        proballObj = sum(proball)
        if proballObj < proBest:
            proBest = proballObj
            procenter = curp
        
    curp = random.uniform(procenter-hk,procenter+hk)  #find a random number in the most promising area
    #initialize tabu list with curp
    tabulist=[]
    for i in xrange(0,tabulength):
        tabulist.append(curp)    
        
    #inversitificated search based on tabu list
    for i in xrange(0,Iter):
        
        ##check elapsed time
        #iteration_time = time.clock()                                   
        #elapsed_time = time.clock() - t0
        #if (elapsed_time > max_time): 
            #break            
        
        #early termination condition check
        if gbestObj_track == gbestObj:
            early_termin_cnt += 1
            if early_termin_cnt >= 10:
                break    #stop loop
        else: early_termin_cnt = 0
                
        gbestObj_track = gbestObj        

        neighborlist = Choose_neighbors(curp,neighbors, h0, hk, tabulist)  #choose neighborhoods around current p
        
        candidate = []
        candidate_obj = []                
        #caculate the objective value of element in neighborlist
        for n in xrange(0,len(neighborlist)):
            m.reset()
            result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, neighborlist[n], decision)
            iterations += result[0]
            candidate.append(result[1])
            candidate_obj.append(result[2])
        
        if candidate_obj == []:
            pass
            #result = DSSP.DSSP (m, arcs, varcost, fixedcost, totSupply, K, flow, random.random()*2, decision)
            #iterations += result[0]
            #candidate.append(result[1])
            #candidate_obj.append(result[2])  
            
        #find the best from candidate_obj and get its p value        
        bestcandidate = min(candidate_obj)
        index = candidate_obj.index(bestcandidate)
        curp = candidate[index]
        
        #compare with best so far  
        if bestcandidate <= gbestObj:
            gbestObj=bestcandidate
            pbest = curp
  
        #update short tabu list with FIFO rule whatever curp is better then best so far
        tabulist = update_tabulist(tabulist,curp)
        
        #long term memory based on frequency of short-term memory
        #difference<0.1, was regarded same p and recoreded in the long term memory
        for pre_curp in long_term_memory.keys():
            if abs(curp-pre_curp) <= h0:
                long_term_memory[pre_curp] += 1
                if long_term_memory[pre_curp] >= long_term_max and bestcandidate >= gbestObj:
                    curp = random.random()*(2-0)
            else:
                long_term_memory.update({curp:1})
                
    TSresult=[iterations,pbest,gbestObj]    
    return TSresult