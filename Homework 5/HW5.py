# QUANTITATIVE MACROECONOMICS: PROBLEM SET 5
# Marta Oliva Riera
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% ===========================================================================
# 1. FACTOR INPUT MISALLOCATION
# =============================================================================
# 1.1 Simulate the population:
n = 10000000
mu = np.transpose([1, 1])
ro = 0
sigma = np.array([[1,ro], [ro,1]])
X = np.random.multivariate_normal(mu, sigma, n)
lnk = X[:,0]
lnz = X[:,1]
k = np.exp(lnk)
z = np.exp(lnz)

plt.hist(k, alpha = 0.5)
plt.hist(z, alpha = 0.5)
plt.title('Joint density in levels')
plt.show()

plt.hist(lnk, bins=10, alpha = 0.5)     # this one doesn't look great but not sure how to fix it
plt.hist(lnz, bins=10, alpha = 0.5)
plt.title('Joint density in logs')
plt.show()

#%% 1.2 Output:
gamma = 0.6
y = z + np.abs(k)**gamma

# 1.3 Solve the maximization problem of output:
K = sum(k)  # aggregate capital from the population data

d = np.empty(n)
kE = np.empty(n)

for i in range(n):
    d[i] = np.abs(z[0]/z[i])**(gamma-1)
kE[0] = K/sum(d)

for j in range(1,n):
    kE[j] = (z[0]/z[j])**(gamma-1)*kE[0] 

# 1.4 Compare ke to the actual values:
plt.hist(k, alpha=0.5)
plt.hist(kE, alpha=0.5)
plt.show()

# 1.5 Output gains from relocation:
Y = sum(y)
yE = z + np.abs(kE)**gamma
YE = sum(yE)

ygain0 = (YE/Y - 1)*100

#%% 1.6 Repeat with different correlation coefficients:
# correlation = 0.5 -----------------------------------------------------------
# 1.1 Simulate the population
ro = 0.5
sigma = np.array([[1,ro], [ro,1]])
X = np.random.multivariate_normal(mu, sigma, n)
k = X[:,0]
z = X[:,1]
lnk = np.exp(k)
lnz = np.exp(z)

# 1.2 Output:
gamma = 0.6
y = z + np.abs(k)**gamma

# 1.3 Solve the maximization problem of output:
K = sum(k)  # aggregate capital from the population data

d = np.empty(n)
kE = np.empty(n)
    
for i in range(n):
    d[i] = np.abs(z[0]/z[i])**(gamma-1)
kE[0] = K/sum(d)

for j in range(1,n):
    kE[j] = (z[0]/z[j])**(gamma-1)*kE[0]    

# 1.4 Compare ke to the actual values:
plt.hist(k, alpha=0.5)
plt.hist(kE, alpha=0.5)
plt.show()

# 1.5 Output gains from relocation:
Y = sum(y)
yE = z + np.abs(kE)**gamma
YE = sum(yE)

ygain1 = (YE/Y - 1)*100

#%% correlation = -0.5 --------------------------------------------------------
# 1.1 Simulate the population
ro = -0.5
sigma = np.array([[1,ro], [ro,1]])
X = np.random.multivariate_normal(mu, sigma, n)
k = X[:,0]
z = X[:,1]
lnk = np.exp(k)
lnz = np.exp(z)

# 1.2 Output:
gamma = 0.6
y = z + np.abs(k)**gamma

# 1.3 Solve the maximization problem of output:
K = sum(k)  # aggregate capital from the population data

d = np.empty(n)
kE = np.empty(n)
    
for i in range(n):
    d[i] = np.abs(z[0]/z[i])**(gamma-1)
kE[0] = K/sum(d)

for j in range(1,n):
    kE[j] = (z[0]/z[j])**(gamma-1)*kE[0] 

# 1.4 Compare ke to the actual values:
plt.hist(k, alpha=0.5)
plt.hist(kE, alpha=0.5)
plt.show()

# 1.5 Output gains from relocation:
Y = sum(y)
yE = z + np.abs(kE)**gamma
YE = sum(yE)

ygain2 = (YE/Y - 1)*100


#%% ===========================================================================
# 2. HIGHER SPAN OF CONTROL
# =============================================================================
# Repeat everything with a higher span of control value, now gamma = 0.8
# 1.1 Simulate the population:
n = 10000000
np.random.seed(42)
mu = np.transpose([1, 1])
ro = 0
sigma = np.array([[1,ro], [ro,1]])
X = np.random.multivariate_normal(mu, sigma, n)
k = X[:,0]
z = X[:,1]
lnk = np.exp(k)
lnz = np.exp(z)

# 1.2 Output:
gamma = 0.8
y = z + np.abs(k)**gamma

# 1.3 Solve the maximization problem of output:
K = sum(k)  # aggregate capital from the population data

d = np.empty(n)
kE = np.empty(n)
    
for i in range(n):
    d[i] = np.abs(z[0]/z[i])**(gamma-1)
kE[0] = K/sum(d)

for j in range(1,n):
    kE[j] = (z[0]/z[j])**(gamma-1)*kE[0] 

# 1.4 Compare ke to the actual values:
plt.hist(k, alpha=0.5)
plt.hist(kE, alpha=0.5)
plt.show()

# 1.5 Output gains from relocation:
Y = sum(y)
yE = z + np.abs(kE)**gamma
YE = sum(yE)

ygain3 = (YE/Y - 1)*100

#%% ===========================================================================
# 3. FROM COMPLETE DISTRIBUTIONS TO RANDOM SAMPLES
# =============================================================================
# 3.1 Take a sample of 10000 observations from the population:
# -----------------------------------------------------------------------------
from random import sample

lnkS = sample(lnk, 10000)
lnzS = sample(lnz, 10000)

sigmakS = np.var(lnkS)
sigmazS = np.var(lnzS)
roS = np.cov(lnkS, lnzS)

# 3.2 Recompute 1.3-1.5 with the sample:
# -----------------------------------------------------------------------------
# 1.3 Solve the maximization problem of output:
gamma = 0.6
K = sum(k)  # aggregate capital from the population data

d = np.empty(n)
kE = np.empty(n)
    
for i in range(n):
    d[i] = np.abs(z[0]/z[i])**(gamma-1)
kE[0] = K/sum(d)

for j in range(1,n):
    kE[j] = (z[0]/z[j])**(gamma-1)*kE[0] 
    
# 1.4 Compare ke to the actual values:
plt.hist(k, alpha=0.5)
plt.hist(kE, alpha=0.5)
plt.show()

# 1.5 Output gains from relocation:
y = z + np.abs(k)**gamma
Y = sum(y)
yE = z + np.abs(kE)**gamma
YE = sum(yE)

ygainS = (YE/Y - 1)*100

# 3.3 Repeat the two previous steps 1000 times, to get different measures of misallocation:
# -----------------------------------------------------------------------------------------
m = 1000
ygainS = np.empty()

for i in range(m):
    # 3.1 Get a new random sample:
    lnkS = sample(lnk, 10000)
    lnzS = sample(lnz, 10000)
    sigmakS = np.var(lnkS)
    sigmazS = np.var(lnzS)
    roS = np.cov(lnkS, lnzS)
    # 3.2 Efficient output, comparison to actual and output gains:
    K = sum(k)
    d = np.empty(n)
    kE = np.empty(n)
        
    for i in range(n):
        d[i] = np.abs(z[0]/z[i])**(gamma-1)
    kE[0] = K/sum(d)
    
    for j in range(1,n):
        kE[j] = (z[0]/z[j])**(gamma-1)*kE[0] 
    
    plt.hist(k, alpha=0.5)
    plt.hist(kE, alpha=0.5)
    plt.show()
    
    y = z + np.abs(k)**gamma
    Y = sum(y)
    yE = z + np.abs(kE)**gamma
    YE = sum(yE)
    
    ygainS[i] = (YE/Y - 1)*100

plt.hist(ygainS)
ygdf = pd.DataFrame(ygainS)
print(ygdf.describe())

# 3.4 Probability that output gain from the samples is within 10% of the actual output gain

