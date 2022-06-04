#!/usr/bin/env python
# coding: utf-8

import gurobipy as gp
import numpy as np
from gurobipy import GRB
import time
import matplotlib.pyplot as plt
import math
import networkx as nx
import sys


#### This algorithm is not presented in the paper

def cpm(N,S,T):
    temp=[]
    for i in range(N):
        for k in S[i]:
            temp.append((i+1,k+1,{"cost":T[i,0]}))
    for i in range(1,N+1):
        temp.append((0,i,{"cost":T[i-1,0]}))
    DG = nx.DiGraph(temp)
    temp = nx.dag_longest_path(DG)
    s=0
    for i in temp:
        if i>0:
            s=s+T[i][0]
    return int(s)




####################### PROBLEM SOLVED USED GUROBI##################################
##### This is a slightly modification of the problem formulation of the paper, which as some problems
#J is the relaxation on the bound of W, H changes the objective function, debug masks outputFlag
def modelgurobi(J,H,debug=0):
    np.random.seed(2022)
######################DISABLE OUTPUT######################
    env = gp.Env(empty=True)
    if debug==0:
        env.setParam("OutputFlag",0)
    env.start()
###########DATA TAKEN FROM THE PAPER##############ààà
    Mg=9000000
    model=gp.Model(env=env)
    FOC=1000
    MOB=10000
    OV=0.1
    AP=15000
    LP=2
    MU=2.64
    R=2
    N=3
    RP=0.25
    V=4
    rW=0.00123
    MU=2.64
    BC=1656
    BP=84456
    S=np.zeros(3,dtype=object)
    M=np.zeros(3,dtype=object)

    S[0]=[1,2]
    S[1]=[]
    S[2]=[]

    C=np.array([[15000,12000,9000],[10000,6000,500000],[7000,500000,500000]])
    T=np.array([[3,4,6],[5,6,500000],[7,500000,500000]])
    M[0]=[0,1,2]
    M[1]=[0]
    M[2]=[0]
    MMF=np.zeros([N,3],dtype=np.int64)
    for i in range(N):
        for j in range(3):
            MMF[i,j]=C[i,0]/C[i,j]
    Wmin=cpm(N,S,T)
    Wmax=Wmin+J
    if H==0:
        W=model.addVar(name="W")
        CL=50000
    
    if H==1:
        CL=model.addVar(name="CL",lb=0)
        W=model.addVar(name="W")

    if H==2:
        W=Wmax
        CL=50000

    if H==3:
        W=Wmax
        CL=50000 
    
    if H==4 or H==5:
        W=Wmax
        CL=50000

    

        
    #Equation from 1 to 20 are not in the model
    x = model.addVars(N,N,Wmax+1, vtype=GRB.BINARY, name="x") #equation 45
    DC = model.addVars(N,Wmax+1, name="dc")
    EV = model.addVars(N,math.ceil(Wmax/R)*R+LP+1, name="ev")
    TE = model.addVars(Wmax+LP+1, name="te")
    I = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="I")
    P = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="P")
    K=math.ceil(Wmax/R)
    y = model.addVars(K+1, vtype=GRB.BINARY, name="y")#equation 47
    ####Mg is the debt limit
    Mg=90000000
    Z=model.addVar(name="Z")
    IB = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="IB")
    B = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="B")
    IL = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="IL") 
    alpha= model.addVars(Wmax+LP+1, vtype=GRB.BINARY, name="alpha") #equation 48
    model.addConstr(W<=Wmax)# equation 22
############ TIME CONSTRAINTS###########
    model.addConstrs(gp.quicksum(x[i,j,k] for j in M[i] for k in range(1,Wmax+1))==1 for i in range(N))#equation 23
    for i in range(N):
        for q in S[i]:
            model.addConstr(gp.quicksum((k+T[i,j])*x[i,j,k] for j in M[i] for k in range(1,Wmax+1))+1<=gp.quicksum(k*x[q,j,k] for j in M[q] for k in range(1,Wmax+1))) #equation 24
        if len(S[i])==0:
            model.addConstr(gp.quicksum((k+T[i,j])*x[i,j,k] for j in M[i] for k in range(1,Wmax+1))<=W) #equation 25
################## CASH OUTFLOW#############
    for i in range(N):
         for k in range(1,Wmax+1):
             model.addConstr(DC[i,k]==gp.quicksum(C[i,j]/T[i,j]*x[i,j,t] for t in range(max(1,k-T[i,j]),k+1) for j in M[i])  ) #equation 26
    model.addConstr(TE[0]==MOB+BC)
    for k in range(1,Wmax+1):
         model.addConstr(TE[k]==FOC+(1+OV)*gp.quicksum(DC[i,k] for i in range(N))) #equation 27
    for k in range(Wmax+1,math.ceil(Wmax/R)*R+1):
         model.addConstr(TE[k]==FOC)
#######CASH INFLOW#####################
    for k in range(1,Wmax+1):
        for i in range(N):
            model.addConstr(EV[i,k]==gp.quicksum(MMF[i,j]*DC[i,k] for j in M[i] for t in range(max(1,k-T[i,j]),k+1)))
            #equation 28
    for k in range(Wmax+1,math.ceil(Wmax/R)*R+1):
        for i in range(N):
            model.addConstr(EV[i,k]==0) #equation 28
    templist=[k for k in range(R,math.ceil(Wmax/R)*R+1,R)]
    
    for k in templist:
        model.addConstr(I[k]==MU*gp.quicksum(EV[i,t] for t in range(k-R+1,k+1) for i in range(N))) #equation 29
    for k in set(range(math.ceil(Wmax/R)*R))-set(templist):
       model.addConstr(I[k]==0)    
    templist=[k for k in range(R+LP,math.ceil(Wmax/R)*R+LP,R)]
    for k in templist:
        model.addConstr(P[k]==(1-RP)*I[k-LP]-gp.quicksum(AP*y[n]/n for n in range(1,K+1))) #equation 30
    for k in set(range(math.ceil(Wmax/R)*R+LP))-set(templist+[0]):    
       model.addConstr(P[k]==0)
    model.addConstr(P[0]==AP)
    model.addConstr(P[math.ceil(Wmax/R)*R+LP]==RP*BP)
########## CASH BALANCE##################
    model.addConstr(R*Z>=W) #equation 31
    model.addConstr(R*Z<=W+R-1) #equation 32
    model.addConstr(gp.quicksum(n*y[n] for n in range(K+1))==Z) #equation 33-34
    model.addConstr(B[0]==AP-TE[0]) #equation 35
    for k in range(1,Wmax+LP+1):
        model.addConstr((B[k]==B[k-1]-TE[k]+P[k-1]-IB[k])) #equation 36
    templist=[k for k in range(V,Wmax,V)]
    templist.append(Wmax)
    for k in templist:
          model.addConstr(IB[k]==rW*gp.quicksum(IL[t] for t in range(max(1,k-V+1),k+1))) #equation 37:
    for k in range(1,Wmax+1+LP):
        model.addConstr(B[k-1]-TE[k-1]<=(1-alpha[k])*Mg )#equation 38
    for k in range(1,Wmax+1+LP):
        model.addConstr(B[k-1]-TE[k-1]>=-alpha[k]*Mg )#equation 39
    for k in range(1,Wmax+1+LP):
         model.addConstr(IL[k]>=0)
         model.addConstr(IL[k]>=-(B[k-1]-TE[k])+(alpha[k]-1)*Mg ) #equation 41
         model.addConstr(IL[k]<=-(B[k-1]-TE[k])+(1-alpha[k])*Mg ) #equation 42
         model.addConstr(IL[k]<=alpha[k]*Mg)#equation 43
         model.addConstr(B[k]<=CL)
    model.addConstr(Z>=0) 
    model.addConstr(Z<=Mg)
    
 
    if H==0:
        model.setObjective(W,GRB.MINIMIZE)   #equation 21

    if H==1:
        model.setObjective(CL,GRB.MINIMIZE)   
    
    if H==2:
        model.setObjective(gp.quicksum(TE[k] for k in range(1,Wmax+1)), GRB.MINIMIZE)
        
    if H==3:
        model.setObjective(B[Wmax+LP], GRB.MAXIMIZE)
    
    
    model.optimize()
    if H==4:
        return len(model.getConstrs())

    if H==5:
        return len(model.getVars())
    
    
    
########################PROBLEM SOLVED USING THE HEURISTIC ALGORITHM###############

##Not present in the paper
#Calculates the activities that can be scheduled
def checkpossibility(statact, S):
    possibleActivities=[]
    possibleActivitiesbase=[i for i in range(len(S)) if statact[i]>0]
    possibleActivitiestemp=[i for i in range(len(S)) if statact[i]==sys.maxsize]
    for i in possibleActivitiestemp:
        flag=1
        for j in possibleActivitiesbase:
            if i in S[j]:
                flag=0
        if flag==1:
            possibleActivities.append(i)

    return possibleActivities
##Not present in the paper
#Calculates the slack
def slack(T,i,Wmax):

    return Wmax-np.max(T[i,:])

##Not present in the paper
#Select the activity to schedule          
def selectActivity(T,possibleActivities,C,Wmax):
    #print(C.astype(int).tolist)
    return possibleActivities[np.lexsort((np.array([slack(T,i,Wmax) for i in possibleActivities]).astype(int).tolist,C.astype(int).tolist))]
       
##Auxiliary function present only in my implementation
def choosemode(Activity,Criteria,C,T):
    if Criteria=='time':
        return np.argmin(T[Activity,:])
    else:
        return np.argmin(C[Activity,:])

##Not present in the paper
def selectCriteria(Activity,F,statact,S,T,Wmax):
    if slack(T,Activity,Wmax)<=0 or len(F)==0:
        return 'time'
    if slack(T,Activity,Wmax)>0:
        for i in F:
            if slack(statact,S,T,i,Wmax)<0:
                return 'time'
        return 'cash'

###This is the heuristic algorithm of the paper with some differences, the most important one are 1)Wmax is fixed, so there no need of terminationCount, and feasibility is checked at every turn
### Also I had to write the internal functions, which are not present in the paper
def heuristic(J):
    Mg=9000000
    FOC=1000
    MOB=10000
    AP=15000
    LP=2
    MU=2.64
    R=2
    N=3
    RP=0.25
    V=4
    rW=0.00123
    MU=2.64
    BC=1656

    S=np.zeros(3,dtype=object)
    M=np.zeros(3,dtype=object)

    S[0]=[1,2]
    S[1]=[]
    S[2]=[]

    C=np.array([[15000,12000,9000],[10000,6000,500000],[7000,500000,500000]])
    T=np.array([[3,4,6],[5,6,500000],[7,500000,500000]])
    M[0]=[0,1,2]
    M[1]=[0]
    M[2]=[0]
    MMF=np.zeros([N,3],dtype=np.int64)
    for i in range(N):
        for j in range(3):
            MMF[i,j]=C[i,0]/C[i,j]
    Wmin=cpm(N,S,T)
    Wmax=Wmin+J
    Mg=9000000
    U=set([i for i in range(len(S))])
    completed=set([])
    F=set([])
    
    
    #Time remaining to the end of the activities, is infinity if they have yet to start
    statact=np.ones(len(S))*(sys.maxsize)
    #Setting financial variables to 0
    IB=np.zeros(math.ceil(Wmax/R)*R+LP+1)
    B=np.zeros(math.ceil(Wmax/R)*R+LP+1)
    P=np.zeros(math.ceil(Wmax/R)*R+LP+1)
    EV=np.zeros([N,math.ceil(Wmax/R)*R+LP+1])
    TE=np.zeros(math.ceil(Wmax/R)*R+LP+1)
    IL=np.zeros(math.ceil(Wmax/R)*R+LP+1) 
    I=np.zeros(math.ceil(Wmax/R)*R+LP+1) 
    TE[0]=MOB+BC
    P[0]=AP
    B[0]=AP-TE[0]
    CL=18000
    timep=0
    while len(U.union(F))!=0 and timep<=Wmax:
        #possibily select and activity
        possibleActivities=checkpossibility(statact, S)
        temp=np.min(C,axis=1)
        if (possibleActivities!=[]):
            nextActivity=selectActivity(T,possibleActivities,temp[possibleActivities],Wmax)
            criteria=selectCriteria(nextActivity,F,statact,S,T,Wmax)
            mode=choosemode(nextActivity, criteria, C, T)
            Tetemp=FOC+C[nextActivity,mode]/T[nextActivity,mode]
            flag=1
        else:
            flag=0
        #updata finances
        if (timep-LP)//R==0:
            P[timep]=(1-RP)*I[timep-LP]-AP/math.ceil(Wmax/R)  
        else:
            P[timep]=0

        if timep//V==0:
            Iltemp=max(Tetemp-B[timep-1],0)
            IBtemp=rW*sum(IL[t] for t in range(timep//V*(V-1), timep))+rW*Iltemp
        else:
            IBtemp=0
        
        Btemp=B[timep-1]+P[timep-1]-Tetemp-IBtemp
        #check financial feasibility
        if Btemp>CL:
            flag=0
        if Btemp<CL and Tetemp-B[timep-1]>Mg:
            flag=0
        # if not feasible there is no direct cost
        if flag==0:
            if possibleActivities!=[]:
                U=U-set([nextActivity])
                F=F-set([nextActivity])
                F=F.union(set([nextActivity]))
            TE[timep]=FOC
            EV[nextActivity,timep]=0
            IL[timep]=0

        #if feasible there is direct cost
        if flag==1:
            TE[timep]=Tetemp
            statact[nextActivity]=T[nextActivity,mode]+1
            EV[nextActivity,timep]=C[nextActivity,0]/T[nextActivity,mode]
            IL[timep]=Iltemp
        ####time passes
        B[timep]=B[timep-1]+P[timep-1]-TE[timep]
        IB[timep]=rW*sum(IL[t] for t in range(timep//V*(V-1), timep))
        IB[timep]=IBtemp
        I[timep]=MU*sum(EV[i,t] for i in range(N) for t in range(max(timep-R-1,0),timep))
        P[timep]=(1-RP)*I[timep-LP]-AP/math.ceil(Wmax/R)
        for i in range(N):
            if statact[i]<sys.maxsize and statact[i]>0:
                #remaining time of running activities decreases
                statact[i]=statact[i]-1
                if statact[i]==0:
                    U=U-set([i])
                    F=F-set([i])
                    completed=completed.union(set([i]))
        timep=timep+1
    if timep>Wmax:
        return -1
    else:
        return timep
                    


######################BENCHMARKING PART ##############
'''
print("Computing scalability of number of constrains and variables on the gurobi model")


constrs=np.zeros(10)
var=np.zeros(10)
for i in range(10):
    print(i)
    constrs[i]=modelgurobi(200*i,4)
    var[i]=modelgurobi(200*i,5)

fig1, axs1 = plt.subplots(2, 1)
axs1[0].plot([i for i in range(0,2000,200)],constrs)
axs1[1].plot([i for i in range(0,2000,200)],var)
axs1[0].set(xlabel="Wmax", ylabel="Number of constraints")
axs1[1].set(xlabel="Wmax", ylabel="Number of variables")
plt.tight_layout()
plt.plot()
fig1.savefig("fig1.png",dpi=900)

'''

print("Computing scalability of time on the gurobi model")

times=np.zeros([10,4])
for i in range(0,10):
    print(i)
    for j in range(4):
        timeold=time.time()
        temp=modelgurobi(i*200,j)
        times[i][j]=time.time()-timeold


fig2, axs2 = plt.subplots(2, 2)
for j in range(4):
    axs2[j//2, j % 2].plot([i for i in range(0, 2000, 200)], times[0:10, j])
    axs2[j//2, j % 2].set(xlabel="Wmax", ylabel="Time")
    axs2[j//2, j % 2].title.set_text('Model '+str(j))
plt.tight_layout()
plt.plot()
fig2.savefig("fig2.png", dpi=900)

fig2, axs2 = plt.subplots(2, 2)
for j in range(4):
    axs2[j//2, j % 2].plot([i for i in range(0, 2000, 200)], np.log(times[0:10, j]))
    axs2[j//2, j % 2].set(xlabel="Wmax", ylabel="Time (log)")
    axs2[j//2, j % 2].title.set_text('Model '+str(j))
plt.tight_layout()
plt.plot()
fig2.savefig("fig2log.png", dpi=900)


print("Computing scalability of time on the heuristic model")
'''
times=np.zeros(10)
for i in range(0,10):
    timeold=time.time()
    temp=heuristic(i*200)
    times[i]=time.time()-timeold
fig3, axs3 = plt.subplots(1, 1)
axs3.plot([i for i in range(0,2000,200)],times)
axs3.set(xlabel="Wmax", ylabel="Time")
plt.plot()
fig3.savefig("fig3.png",dpi=900)

'''