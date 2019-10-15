# QUANTITATIVE MACROECONOMICS: PROBLEM SET 4
# Marta Oliva Riera
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import quantecon as qe

#%% ===========================================================================
# QUESTION 1: VALUE FUNCTION ITERATION
# =============================================================================
# 2.1 Recursive formulation without productivity shocks
#------------------------------------------------------------------------------
# Parameters: 
theta = 0.679
beta = 0.988
delta = 0.013
kappa = 5.24
nu = 2

# Inelastic labour supply:
h = 1

# Functions (in steady state)
def f(k, h):
    return k**(1-theta)*h**theta

def u(c, h):
    return  np.log(c) - kappa*((h**(1+(1/nu)))/(1+(1/nu)))

# (a) Brute Force iterations --------------------------------------------------
qe.tic()
## Discretize the state space
N_k = 100
k = np.linspace(0.01, 50, N_k)      # evenly spaced grid points
ki, kj = np.meshgrid(k,k)           # grid with all possible combinations of k values

## Define the value function matrix and make an initial guess:
V = np.empty([N_k, 300])
V[:,0] = np.zeros((N_k))      #initial guess


## Define the initial return matrix (with feasible values of c and k)
def c(k1, k2):
    return k1**(1-theta)*h**(theta) + (1-delta)*k1 - k2
C = c(ki, kj)   #evaluated at the combinations of kt and kt+1

def U(k1,k2):
    for i in range(N_k):
        for j in range(N_k):
            if C[i,j] >= 0 :
                return np.log(k1**(1-theta)*h**(theta) + (1-delta)*k1 - k2) - kappa*((h**(1+(1/nu)))/(1+(1/nu)))

M = U(ki,kj)
M[np.isnan(M)] = -1000

## Define the value function iteration:
X = np.empty([N_k, 300])
G = np.empty([N_k, 300])
n = 0

for s in range(299):    # loop over number of iterations
    e = 0.01
    # computing matrix X:
    for i in range(N_k):
        for j in range(N_k):
            X[i,j] = M[i,j] + (beta*(V[:,s][j]))
    # updating the value function, getting the policy functions:
    for i in range(N_k):
        V[:,s+1] = np.amax(X[i,:])
        G[:, s] = np.argmax(X[i,:])
        #checking convergence of the value functions:
        for i in range(N_k):
            if abs(V[:,s+1][i] - V[:,s][i]) > e:
                continue
            else:
                n += 1
                break

qe.toc()

plt.plot(V, label='V(k)')
plt.show()

