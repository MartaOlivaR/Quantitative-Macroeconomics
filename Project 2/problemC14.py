# QUANTITATIVE MACROECONOMICS: PROJECT 2
# Marta Oliva Riera
import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe

#%% ===========================================================================
# PROBLEM C: VALUE FUNCTION ITERATION
# =============================================================================
# parameters + labour supply:
r = 0.02
rho = 0.03
g = 0.01
theta = 1
beta = 1/(1+rho)

# (a) Brute force Value Function Iteration 
#------------------------------------------------------------------------------
# Discretize the state space: 
qe.tic()

nx = 50
grd = np.empty(nx)
def makegrid(x1, x2, nc, curv):
    scale = x2-x1
    grd[0] = x1
    grd[nc-1] = x2
    for ic in range(1,nc-1):
        grd[ic] = x1 + scale*((ic-1.0)/(nc-1.0))**curv
    return grd


x = makegrid(0.01, 30, nx, 3.0)   # nk evenly spaced points

# Create the income shocks:
ne = 7
varepsi = 0.01
muepsi = -varepsi/2
epsi, probepsi = qe.quad.qnwnorm(ne, muepsi, varepsi)
epsi = np.exp(epsi)

# Utility function:
def u(c):
    if theta == 1:
        u = np.log(c)
    else:
        u = (1/(1-theta))*c**(1-theta)
    return u

# Returns matrix:
       
M = np.empty([nx,nx, ne])
for i,xi in enumerate(x):
        for j, xj in enumerate(x):
            for I, eI in enumerate(epsi):
                for J, eJ in enumerate(epsi):
                    c = xi - (xj - eJ)*((1+g)/(1+r))    # consumption in terms of cash at hand, from the BC 
                    if  c >= 0:
                        if xi >= np.any(c) :   # using only values where the borr. constraint holds (x>c)
                            M[i,j,I]  = u(c)
                        else:
                            M[i,j,I] = u(xi) # if not binding, make x = c, so u(c) = u(x)
                    else:
                        M[i,j] = -1000000

# Value function iteration:
epsilon = 1e-4 
convergence = False
s = 0
S = 800
X = np.empty([nx, nx, ne])
G = np.ones([nx, ne])
Vi = np.zeros([nx,ne])   # initial value function guess
Vj = np.empty([nx,ne])
V = np.empty([nx,ne])

while (s<S) and (convergence == False):
    for i in range(nx):
        for j in range(nx):
            for I in range(ne):
                V = np.mean(probepsi[I]*Vi[j,I])    # the value function will be the mean over all states of the shock
                X[i,j,I] = M[i,j,I] + (1+g)**(1-theta)*beta*(V)
    for i in range(nx):
        for I in range(ne):
            Vj[i,I] = np.max(X[i,:,I])
            G[i,I] = np.argmax(X[i,:,I])
    if np.max(Vj - Vi) >= epsilon:
        Vi = np.copy(Vj)
        Vj = np.empty([nx,ne])
        s += 1
    else:
        convergence = True
        
T = qe.toc()

# Plotting:
plt.plot(x, Vj, label='V(x)')       
plt.xlabel('x')
plt.ylabel('V(x)')
plt.show()

print("Convergence of the value functions took "+str(s)+" iterations, in " + str(T) +" seconds.")  

#%% Policy functions:
gx = np.empty(nx)
gc = np.empty(nx)

for i in range(nx):
    gx[i] = x[int(G[i,I])]
for i in range(1, nx):
    for j in range(1, ne):
        gc[i] = gx[i] - (gx[i-1] - epsi[j])*((1+g)/(1+r))
    
#plt.plot(x, gx)
plt.plot(gc)
    



    
