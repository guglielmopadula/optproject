#!/usr/bin/env python
# coding: utf-8

import gurobipy as gp
import numpy as np
from gurobipy import GRB
import time
import matplotlib.pyplot as plt
import math
import networkx as nx
from sklearn.linear_model import LinearRegression



def cpm(N,S,T):
    temp=[]
    for i in range(N-1):
        for k in S[i]:
            temp.append((i,k,{"cost":T[i,0]}))
    DG = nx.DiGraph(temp)
    temp=nx.dag_longest_path(DG)
    s=0
    for i in temp:
        s=s+T[i][0]
    return int(s)



#Model has fixed size
def model(J,H,N):
    np.random.seed(2022)
######################DATA PREPARATION######################
    env = gp.Env(empty=True)
    env.setParam("OutputFlag",0)
    env.start()

    Mg=90000
    model=gp.Model(env=env)
    FOC=100+np.random.rand()*(1000-100)
    MOB=2000+np.random.rand()*(20000-2000)
    OB=np.random.rand()
    OP=np.random.rand()
    OV=np.random.rand()
    AP=5000+(50000-5000)*np.random.rand()
    
   
#Parameters taken by the data (r1)
    S=np.zeros(N, dtype=object)
    for i in range(N-1):
        temp=np.random.randint(0,N-i-1,1)
        S[i]=sorted(np.random.randint(i+1,N,temp).tolist())
    S[N-1]=[]
    M=np.zeros(N, dtype=object)
    for i in range(N):
        M[i]=[0,1,2]
    T=np.zeros([N,3],dtype=np.int64)
    for i in range(N):
        T[i,0]=np.random.randint(1,5,1)
        for j in range(1,3):
            T[i,j]=np.random.randint(T[i,0],5,1)
        C=np.zeros([N,3],dtype=np.int64)
    for i in range(N):
        for j in range(3):
            C[i][j]=np.random.randint(1,5,1)
    MMF=np.zeros([N,3],dtype=np.int64)
    for i in range(N):
        for j in range(3):
            MMF[i,j]=C[i,0]/C[i,j]
    Wmin=cpm(N,S,T)
    Wmax=Wmin+J

    if H==0:
        W=model.addVar(name="W")
        CL=18000
    
    if H==1:
        CL=model.addVar(name="CL",lb=0)
        W=model.addVar(name="W")

    if H==2:
        W=Wmax
        CL=18000

    if H==3:
        W=Wmax
        CL=18000 
    
    if H==4 or H==5:
        W=Wmax
        CL=18000

    
    
    LP=np.random.randint(1,Wmax)
    LR=np.random.randint(1,Wmax)
    R=np.random.randint(1,Wmax)
    RP=np.random.rand()
    rW=np.random.rand()
    V=np.random.randint(1,Wmax)


    DCsum=np.sum(C[:,0])#Equation 1 
    VOCsum=OV*DCsum# Equation 2
    FOCsum=DCsum*Wmax# Equation 3
    PM=OP*(DCsum+VOCsum+FOCsum+MOB)#Equation 4
    BC=OB*(DCsum+VOCsum+FOCsum+MOB+PM)#Equation 5
    BP=DCsum+VOCsum+FOCsum+MOB+PM+BC#Equation 6
    MU=BP/DCsum#Equation 7

          
    
    #Equation from 7 to 20 are not in the model and equation from 21 to 22 are not used because of the fixed sizde
    x = model.addVars(N,N,Wmax+1, vtype=GRB.BINARY, name="x") #equation 45
    DC = model.addVars(N,Wmax+1, name="dc")
    EV = model.addVars(N,math.ceil(Wmax/R)*R+1, name="ev")
    TE = model.addVars(Wmax+1, name="te")
    I = model.addVars(math.ceil(Wmax/R)*R+1, name="I")
    P = model.addVars(math.ceil(Wmax/R)*R+LP+1, name="P")
    K=math.ceil(Wmax/R)
    y = model.addVars(K+1, vtype=GRB.BINARY, name="y")#equation 47
    temp=model.addVar(vtype=GRB.BINARY, name="temp")
    Mg=model.addVar(name="Mg")
    Z=model.addVar(name="Z")
    IB = model.addVars(Wmax+1, name="IB")
    B = model.addVars(Wmax+1, name="B")
    IL = model.addVars(Wmax+1, name="IL") 
    alpha= model.addVars(Wmax+1, vtype=GRB.BINARY, name="alpha") #equation 48

    model.addConstrs(gp.quicksum(x[i,j,k] for j in M[i] for k in range(1,Wmax+1))==1 for i in range(N))#equation 23
    model.addConstrs(gp.quicksum((k+T[i,j])*x[i,j,k] for j in M[i] for k in range(1,Wmax+1))<gp.quicksum(k*x[q,j,k] for j in M[q] for k in range(1,Wmax+1)) for q in S[i] for i in range(N)) #equation 24
    model.addConstrs(gp.quicksum((k+T[i,j])*x[i,j,k] for j in M[i] for k in range(1,Wmax+1))<=W for i in range(N) if len(S[i])==0) #equation 25
    model.addConstrs(DC[i,k]==gp.quicksum(C[i,j]/T[i,j]*x[i,j,t] for t in range(max(1,k-T[i,j]),k+1) for j in M[i]) for i in range(N) for k in range(1,Wmax+1)) #equation 26
    model.addConstrs(TE[k]==FOC+(1+OV)*gp.quicksum(DC[i,k] for i in range(N)) for k in range(1,Wmax+1)) #equation 27
    model.addConstr(TE[0]==MOB+BC)
    model.addConstrs(EV[i,k]==gp.quicksum(MMF[i,j]*DC[i,k] for j in M[i] for t in range(max(1,k-T[i,j]),k+1)) for i in range(N) for k in range(1,Wmax+1)) #equation 28
    model.addConstrs(EV[i,k]==0 for i in range(N) for k in range(Wmax+1,math.ceil(Wmax/R)*R+1)) #equation 28
    templist=[k for k in range(R,math.ceil(Wmax/R)*R,R)]
    templist.append(math.ceil(Wmax/R)*R)
    model.addConstrs(I[k]==MU*gp.quicksum(EV[i,t] for t in range(k-R+1,k+1) for i in range(N)) for k in templist) #equation 29
    model.addConstrs(I[k]==0 for k in set(range(math.ceil(Wmax/R)*R))-set(templist))    #templist=[k for k in range(R+LP,math.ceil(Wmax/R)*R,R)]
    templist=[k for k in range(R+LP,math.ceil(Wmax/R)*R,R)]
    templist.append(math.ceil(Wmax/R)*R+LP)
    model.addConstrs(P[k]==RP*I[k-LP]-gp.quicksum(AP*y[n]/n for n in range(1,K+1)) for k in templist) #equation 30
    model.addConstrs(P[k]==0 for k in set(range(math.ceil(Wmax/R)*R+LP+1))-set(templist)) 
    model.addConstr(R*Z>=W) #equation 31
    model.addConstr(R*Z<=W+R-1) #equation 32
    model.addConstr(gp.quicksum(n*y[n] for n in range(K+1))==Z) #equation 33-34
    model.addConstr(B[0]==AP-TE[0]) #equation 35
    model.addConstrs(B[k]==B[k-1]-TE[k]+P[k-1]-IB[k] for k in range(1,Wmax+1)) #equation 36
    templist=[k for k in range(V,Wmax,V)]
    templist.append(Wmax)
    model.addConstrs(IB[k]==rW*gp.quicksum(IL[t] for t in range(max(1,k-V+1),k+1)) for k in templist) #equation 37
    model.addConstrs(B[k-1]-TE[k-1]<=(1-alpha[k])*Mg for k in range(1,Wmax+1))#equation 38
    model.addConstrs(B[k-1]-TE[k-1]>=-alpha[k]*Mg for k in range(1,Wmax+1))#equation 39
    model.addConstrs(IL[k]>=0  for k in range(1,Wmax+1))
    model.addConstrs(IL[k]>=-(B[k-1]-TE[k])+(alpha[k]-1)*Mg for k in range(1,Wmax+1)) #equation 41
    model.addConstrs(IL[k]<=-(B[k-1]-TE[k])+(1-alpha[k])*Mg for k in range(1,Wmax+1)) #equation 42
    model.addConstrs(IL[k]<=alpha[k]*Mg  for k in range(1,Wmax+1))#equation 43
    model.addConstrs(B[k]<=CL for k in range(1,Wmax+1))
    model.addConstr(Z>=0) 
    model.addConstr(Z<=Mg)
 
    if H==0:
        model.setObjective(W,GRB.MINIMIZE)   

    if H==1:
        model.setObjective(CL,GRB.MINIMIZE)   
    
    if H==2:
        model.setObjective(gp.quicksum(TE[k] for k in range(1,Wmax+1)), GRB.MINIMIZE)
        
    if H==3:
        model.setObjective(B[Wmax], GRB.MAXIMIZE)
    
        
    model.optimize()

    if H==4:
        return len(model.getConstrs())

    if H==5:
        return len(model.getVars())





times=np.zeros([10,4])
for i in range(0,10):
    print(i)
    for j in range(4):
        timeold=time.time()
        temp=model(i*100,j,20)
        times[i][j]=time.time()-timeold

print("done1")

constrs=np.zeros(10)
var=np.zeros(10)
for i in range(10):
    constrs[i]=model(100*i,4,20)
    var[i]=model(100*i,5,20)

fig1, axs1 = plt.subplots(2, 1)
axs1[0].plot([i for i in range(0,1000,100)],constrs)
axs1[1].plot([i for i in range(0,1000,100)],var)
axs1[0].set(xlabel="Wmax", ylabel="Number of constraints")
axs1[1].set(xlabel="Wmax", ylabel="Number of variables")
plt.tight_layout()
plt.plot()
fig1.savefig("fig1.png",dpi=900)


fig2, axs2 = plt.subplots(2, 2)
for j in range(4):
    axs2[j//2,j%2].plot([i for i in range(0,1000,100)],times[0:10,j])
    axs2[j//2,j%2].set(xlabel="Wmax", ylabel="Time")
    axs2[j//2,j%2].title.set_text('Model '+str(j))
plt.tight_layout()
plt.plot()
fig2.savefig("fig2.png",dpi=900)

print("done2")


times2=np.zeros([10,4])
for i in range(0,10):
    print(i)
    for j in range(4):
        timeold=time.time()
        temp=model(20,j,i*10+4)
        times2[i][j]=time.time()-timeold
plt.tight_layout()
plt.plot()

print("done3")


constrs2=np.zeros(10)
var2=np.zeros(10)
for i in range(10):
    constrs2[i]=model(20,4,i*10+4)
    var2[i]=model(20,5,i*10+4)
    
print("done4")


fig3, axs3 = plt.subplots(2, 1)
axs3[0].plot([i for i in range(0,100,10)],constrs2)
axs3[1].plot([i for i in range(0,100,10)],var2)
axs3[0].set(xlabel="N", ylabel="Number of constraints")
axs3[1].set(xlabel="N", ylabel="Number of variables")
plt.tight_layout()
plt.plot()
fig3.savefig("fig3.png",dpi=900)



fig4, axs4 = plt.subplots(2, 2)
for j in range(4):
    axs4[j//2,j%2].plot([i for i in range(0,100,10)],times2[0:10,j])
    axs4[j//2,j%2].set(xlabel="N", ylabel="Time")
    axs4[j//2,j%2].title.set_text('Model '+str(j))
plt.tight_layout()
plt.plot()
fig4.savefig("fig4.png",dpi=900)



