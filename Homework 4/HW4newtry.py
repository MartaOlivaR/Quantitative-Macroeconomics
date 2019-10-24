# QUANTITATIVE MACROECONOMICS: PROBLEM SET 4 - VALUE FUNCTION ITERATION
# Marta Oliva Riera
import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe

#%% ===========================================================================
# 1. VFI with inelastic labour supply
# =============================================================================
# parameters + labour supply:
theta = 0.679
beta = 0.988
delta = 0.013
kappa = 5.24
nu = 2.0

h = 1   # inelastic labour supply

# (a) Brute force Value Function Iteration 
#------------------------------------------------------------------------------
# Discretize the state space: 
qe.tic()
kSS = ((1/beta-1+delta)/(1-theta))**(-1/theta)
nk = 200
k = np.linspace(1, 1.5*kSS, nk)   # nk evenly spaced points

# Returns matrix:
M = np.empty([nk,nk])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                M[i,j]  = np.log(ki**(1- theta) + (1-delta)*ki - kj)
            else:
                M[i,j] = -1000000

# Value function iteration:
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            X[i,j] = M[i,j] + beta*Vi[j]
    for i in range(nk):
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

# Policy functions:
gk = np.empty(nk)
gc = np.empty(nk)
for i in range(nk):
    gk[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gk[i]


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')       
plt.xlabel('k')
plt.ylabel('V(k)')
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  
    

#%% (b) Monotonicity of the optimal decision rule 
#------------------------------------------------------------------------------
# Value function iteration:
qe.tic()
epsilon = 0.01 
convergence = False
s = 0
S = 500
g = np.ones([nk, 1])    
X = np.zeros([nk, nk])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])   # empty updated value function

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            gLB = int(g[i])    # lower bound of the decision rule
            if j+gLB < nk:     # monotonicity of the decision rule
                X[i,j+gLB] = M[i,j+gLB] + beta*Vi[j+gLB]
            else:
                continue
    for i in range(nk):
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')     
plt.xlabel('k')
plt.ylabel('V(k)')  
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  
    
#%% (c) Concavity of the value function 
#------------------------------------------------------------------------------
qe.tic()
# Value function iteration:
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.zeros([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            X[i,j] = M[i,j] + beta*Vi[j]
            if X[i,j] < X[i,j-1]:   # concavity property
                break
    for i in range(nk):
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')   
plt.xlabel('k')
plt.ylabel('V(k)')    
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  


#%% (d) Local search on the decision rule 
#------------------------------------------------------------------------------
qe.tic()
# Value function iteration:
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.zeros([nk, nk])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            if (j >= g[i]) and (j <= g[i] + 5) :    # local search
                X[i,j] = M[i,j] + beta*Vi[j]
    for i in range(nk):
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')
plt.xlabel('k')
plt.ylabel('V(k)')       
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  


#%% (e) Concavity + Monotonicity 
#------------------------------------------------------------------------------
# Value function iteration:
qe.tic()
epsilon = 0.01 
convergence = False
s = 0
S = 500
g = np.ones([nk, 1])    
X = np.zeros([nk, nk])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])   # empty updated value function

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            gLB = int(g[i])    # lower bound of the decision rule
            if j+gLB < nk:     # monotonicity of the decision rule
                X[i,j+gLB] = M[i,j+gLB] + beta*Vi[j+gLB]
                if X[i,j] < X[i,j-1]:   # concavity property
                    break
            else:
                continue
    for i in range(nk):
        Vj[i] = np.max(X[i,:])
        g[i] = np.argmax(X[i,:])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True

T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')  
plt.xlabel('k')
plt.ylabel('V(k)')     
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  


#%% (f) Howard's policy function iterations 
#------------------------------------------------------------------------------
qe.tic()

# discretize the state space:
kSS = ((1/beta-1+delta)/(1-theta))**(-1/theta)
nk = 100
k = np.linspace(1, 1.5*kSS, nk)   # nk evenly spaced points

# Initial guesses and empty matrices:
g0 = np.transpose(np.arange(0,100))    # Initial decision rule guess
g = np.empty([nk, 1])  # empty matrix to update the guess

gki = np.empty([nk, 1])
gkj = np.empty([nk, 1])
gc = np.empty([nk, 1])

# Policy function iteration:
epsilon = 0.01 
convergence = False
s = 0
S = 500
X = np.empty([nk, nk])
V0 = np.zeros([nk,1])
V = np.empty([nk,1])

"""# 1. use the inital guess to compute the values associated to it (gk and gc) + use 
# those in the value function:

for i in range(nk):
    gki[i] = k[int(g0[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gki[i]
    V[i] = np.max(np.log(gc[i]) + beta*V0[i])
    g[i] = np.argmax(np.log(gc[i]) + beta*V0[i])   # new decision rule
  
# 2. find the new policy function (with the new decision rule):
for i in range(nk):
    gkj[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gkj[i]    

# 3. Check for convergence of the policy functions:
np.max(np.abs(gkj - gki)) < epsilon

# REPEAT:
gki = np.copy(gkj)
gkj = np.empty([nk,1])
g0 = np.copy(g)
g = np.empty([nk, 1])

# 1. use the inital guess to compute the values associated to it (gk and gc) + use 
# those in the value function:

for i in range(nk):
    gki[i] = k[int(g0[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gki[i]
    V[i] = np.max(np.log(gc[i]) + beta*V0[i])
    g[i] = np.argmax(np.log(gc[i]) + beta*V0[i])   # new decision rule
  
# 2. find the new policy function (with the new decision rule):
for i in range(nk):
    gkj[i] = k[int(g[i])]
    gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gkj[i]    

# 3. Check for convergence of the policy functions:
np.max(np.abs(gkj - gki)) < epsilon """


# Policy iteration: 
while (s < S) and (convergence == False):
    for i in range(nk):
        gki[i] = k[int(g0[i])]
        gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gki[i]
        V[i] = np.max(np.log(gc[i]) + beta*V0[i])
        g[i] = np.argmax(np.log(gc[i]) + beta*V0[i])   # new decision rule
    for i in range(nk):
        gkj[i] = k[int(g[i])]
        gc[i] = k[i]**(1- theta) + (1-delta)*k[i] - gkj[i]    
    if np.max(gkj - gki)>= epsilon:
        gki = np.copy(gkj)
        gkj = np.empty([nk,1])
        s += 1
    else:
        convergence = True


#%% ===========================================================================
# 2. VFI with labour choice
# =============================================================================
# parameters + labour supply:
theta = 0.679
beta = 0.988
delta = 0.013
kappa = 5.24
nu = 2.0

# Brute force Value Function Iteration 
#------------------------------------------------------------------------------
# Discretize the state space, now also for h: 
qe.tic()
kSS = ((1/beta-1+delta)/(1-theta))**(-1/theta)
nk = 100
nh = 50
k = np.linspace(0.5, 1.5*kSS, nk)   # nk evenly spaced points
h = np.linspace(0.01, 1, nh)   # nh evenly spaced points

# Returns matrix:
M = np.empty([nk,nk, nh])
for i,ki in enumerate(k):
        for j, kj in enumerate(k):
            for l, hl in enumerate(h):
                if kj <= ki**(1- theta) + (1-delta)*ki:   # using only feasible values (c>0)
                    M[i,j,l]  = np.log(ki**(1- theta) + (1-delta)*ki - kj) + kappa*(hl**(1+(1/nu))/(1+(1/nu)))
                else:
                    M[i,j,l] = -1000000

# Value function iteration:
epsilon = 0.01 
convergence = False
s = 0
S = 800
X = np.empty([nk, nk, nh])
g = np.ones([nk, 1])
Vi = np.zeros([nk,1])   # initial value function guess
Vj = np.empty([nk,1])
gk = np.empty(nk)
gh = np.empty(nh)
gc = np.empty(nk)

while (s<S) and (convergence == False):
    for i in range(nk):
        for j in range(nk):
            X[i,j] = M[i,j] + beta*Vi[j]
    for i in range(nk):
        Vj[i] = np.max(X[i,:,:])
        gk[i] = k[np.unravel_index(np.argmax(X[i,:,:], axis=None), X[i,:,:].shape)[0]]
        gc[i] =  k[i]**(1- theta) + (1-delta)*k[i] - gk[i]    
    for i in range(nh):
        gh[i] = h[np.unravel_index(np.nanargmax(X[i,:,:], axis=None), X[i,:,:].shape)[1]]
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nk,1])
        s += 1
    else:
        convergence = True
        


T = qe.toc()

# Plotting:
plt.plot(k, Vj, label='V(k)')     
plt.xlabel('k')
plt.ylabel('V(k)')    
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  
