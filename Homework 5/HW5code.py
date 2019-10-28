# QUANTITATIVE MACROECONOMICS: PROBLEM SET 5
# Marta Oliva Riera
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%% ===========================================================================
# 1. FACTOR INPUT MISALLOCATION
# =============================================================================
ygain1 = np.empty(3) # computing the three correlation cases at the same time
rho = [0,0.5,-0.5]
for ri,r in enumerate(rho):
    print("Rho = " + str(r))
    # 1.1 Simulate the population:
    n = 10000000
    mu = [1, 1]
    sigma = np.array([[1,r], [r,1]])
    X = np.random.multivariate_normal(mu, sigma, n)
    lnk = X[:,0]
    lnz = X[:,1]
    k = np.exp(lnk)
    z = np.exp(lnz)
    
    plt.hist(k, bins=15, alpha = 0.5)   # this one doesn't look great but not sure how to fix it
    plt.hist(z, bins=15, alpha = 0.5)
    plt.title('Joint density in levels')
    plt.show()
    
    plt.hist(lnk, alpha = 0.5)
    plt.hist(lnz, alpha = 0.5)
    plt.title('Joint density in logs')
    plt.show()
    
    # 1.2 Output:
    gamma = 0.6
    y = z*((k)**gamma)
    
    # 1.3 Solve the maximization problem of output:
    K = sum(k)  # aggregate capital from the population data
    
    d = np.empty(n)
    kE = np.empty(n)
    
    for i in range(n):
        d[i] = (z[0]/z[i])**(1/(gamma-1))
    
    kE[0] = K/sum(d)
    kE = d*kE[0]
    
    # 1.4 Compare ke to the actual values:
    plt.plot(kE-k)
    plt.show()
    
    # 1.5 Output gains from relocation:
    Y = sum(y)
    yE = z*((kE)**gamma)
    YE = sum(yE)
    
    ygain1[ri] = ((YE/Y) - 1)*100

#%% ===========================================================================
# 2. HIGHER SPAN OF CONTROL
# =============================================================================
# Repeat everything with a higher span of control value, now gamma = 0.8
ygain2 = np.empty(3) # computing the three correlation cases at the same time
ro = [0,0.5,-0.5]
for ri,r in enumerate(ro):
    # 1.1 Simulate the population:
    n = 10000000
    mu = [1, 1]
    sigma = np.array([[1,r], [r,1]])
    X = np.random.multivariate_normal(mu, sigma, n)
    lnk = X[:,0]
    lnz = X[:,1]
    k = np.exp(lnk)
    z = np.exp(lnz)
    
    plt.hist(k, bins=15, alpha = 0.5)   # this one doesn't look great but not sure how to fix it
    plt.hist(z, bins=15, alpha = 0.5)
    plt.title('Joint density in levels')
    plt.show()
    
    plt.hist(lnk, alpha = 0.5)
    plt.hist(lnz, alpha = 0.5)
    plt.title('Joint density in logs')
    plt.show()
    
    # 1.2 Output:
    gamma = 0.8
    y = z*((k)**gamma)
    
    # 1.3 Solve the maximization problem of output:
    K = sum(k)  # aggregate capital from the population data
    
    d = np.empty(n)
    kE = np.empty(n)
    
    for i in range(n):
        d[i] = (z[0]/z[i])**(1/(gamma-1))
    
    kE[0] = K/sum(d)
    kE = d*kE[0]
    
    # 1.4 Compare ke to the actual values:
    plt.plot(kE-k)
    plt.show()
    
    # 1.5 Output gains from relocation:
    Y = sum(y)
    yE = z*((kE)**gamma)
    YE = sum(yE)
    
    ygain2[ri] = ((YE/Y) - 1)*100

#%% Table summarizing all output gains in Q1 + Q2:
# -----------------------------------------------------------------------------
YGain = {'Output gains': ['gamma = 0.6', 'gamma = 0.8'], 'ro = 0': [ygain1[0], ygain2[0]], 'ro = 0.5': [ygain1[1], ygain2[1]], 'ro = -0.5': [ygain1[2], ygain2[2]]}
df1 = pd.DataFrame(YGain)
print(df1)

#%% ===========================================================================
# 3. FROM COMPLETE DISTRIBUTIONS TO RANDOM SAMPLES
# =============================================================================
# 3.1 Take a sample of 10000 observations from the population:
# -----------------------------------------------------------------------------
# Resetting the population values to those with no correlation:
n = 10000000
ro = 0
mu = [1, 1]
sigma = np.array([[1,ro], [ro,1]])
X = np.random.multivariate_normal(mu, sigma, n)
lnk = X[:,0]
lnz = X[:,1]
k = np.exp(lnk)
z = np.exp(lnz)

# Take the random sample:
from numpy.random import choice
m = 10000
lnkS = choice(lnk, m)
kS = np.exp(lnkS)
lnzS = choice(lnz, m)
zS = np.exp(lnzS)

# Check the variance and covariance of the sample:
sigmakS = np.var(lnkS)
sigmazS = np.var(lnzS)
roS = np.cov(lnkS, lnzS)
print("The mean of k and z in the sample are: " + str(sigmakS) + ", " + str(sigmazS))
print("The covariance matrix of the sample is: " + str(roS))

# 3.2 Recompute 1.3-1.5 with the sample:
# -----------------------------------------------------------------------------
# 1.3 Solve the maximization problem of output:
gamma = 0.6
KS = sum(kS)  
d = np.empty(m)
kES = np.empty(m)

for i in range(m):
    d[i] = (zS[0]/zS[i])**(1/(gamma-1))

kES[0] = KS/sum(d)
kES = d*kES[0]

# 1.4 Compare ke to the actual values:
plt.plot(kES-kS)
plt.show()

# 1.5 Output gains from relocation:
yS = zS*((kS)**gamma)
YS = sum(yS)
yES = zS*((kES)**gamma)
YES = sum(yES)

ygainS = ((YES/YS) - 1)*100

# 3.3 Repeat the two previous steps 1000 times, to get different measures of misallocation:
# -----------------------------------------------------------------------------------------
N = 1000
ygainS = np.empty(N)

for j in range(N):
    # 3.1 Get a new random sample:
    lnkS = choice(lnk, m)
    kS = np.exp(lnkS)
    lnzS = choice(lnz, m)
    zS = np.exp(lnzS)
    # 3.2 Efficient output and output gains:
    KS = sum(kS)  
    d = np.empty(m)
    kES = np.empty(m)
    
    for i in range(m):
        d[i] = (zS[0]/zS[i])**(1/(gamma-1))
    
    kES[0] = KS/sum(d)
    kES = d*kES[0]
                                # skipping the plotting here because it'd take too long
    yS = zS*((kS)**gamma)
    YS = sum(yS)
    yES = zS*((kES)**gamma)
    YES = sum(yES)
    
    ygainS[j] = (YES/YS - 1)*100

plt.hist(ygainS)
ygdf = pd.DataFrame(ygainS)
print(ygdf.describe())

# 3.4 Probability that output gain from the samples is within 10% of the actual output gain:
# --------------------------------------------------------------------------------------------
# create an interval of +-10% around output gains from the population:
LB = 0.9*ygain1[0]
UB = 1.1*ygain1[0]

# count the amount of times the output gain in the sample is within the interval:
count = 0
for i in ygainS:
    if (i>=LB) and (i<=UB):
        count += 1

# compute the probability:
prob = count/N
print("Probability of being in a 10% interval of the actual output gains: " + str(prob) +" %")

#%% 3.5 Repeat with different sample-to-population ratios:
#------------------------------------------------------------------------------
samplesize = [100, 1000, 10000, 100000]
N = 1000
ygainS = np.empty([1000,4])
prob = np.empty(4)
yg = np.empty([N,4])

for si, s in enumerate(samplesize):
    # 3.3 Repeat 1000 times to get different measures of output gains from relocation:
    for j in range(N):
        # 3.1 Get a new random sample:
        lnkS = choice(lnk, s)
        kS = np.exp(lnkS)
        lnzS = choice(lnz, s)
        zS = np.exp(lnzS)
        # 3.2 Efficient output and output gains:
        KS = sum(kS)  
        d = np.empty(s)
        kES = np.empty(s)
        
        for i in range(s):
            d[i] = (zS[0]/zS[i])**(1/(gamma-1))
        
        kES[0] = KS/sum(d)
        kES = d*kES[0]

        yS = zS*((kS)**gamma)
        YS = sum(yS)
        yES = zS*((kES)**gamma)
        YES = sum(yES)
        
        ygainS[j,si] = (YES/YS - 1)*100
    # Plot the output gains for the different random samples + print descriptive statistics:
    plt.hist(ygainS[:,si])
    plt.show()
    print("Descriptive statistics: " + str(pd.DataFrame(yg[:,si]).describe()))
    
    #3.4 Probability that output gain from the samples is within 10% of the actual output gain:
    # create an interval of +-10% around output gains from the population:
    LB = 0.9*ygain1[0]
    UB = 1.1*ygain1[0]
    # count the amount of times the output gain in the sample is within the interval:
    count = 0
    for i in ygainS[:,si]:
        if (i>=LB) and (i<=UB):
            count += 1
    # compute the probability:
    prob[si] = (count/N)*100
    print("Probability of being in a 10% interval of the actual output gains: " + str(prob[si]) +" %")

#%% Summary of the results:
plt.hist(ygainS[:,0], label='sample 100', alpha=0.8)
plt.hist(ygainS[:,1], label='sample 1000', alpha=0.8)
plt.hist(ygainS[:,2], label='sample 10000', alpha=0.8)
plt.hist(ygainS[:,3], label='sample 10000', alpha=0.8)
plt.legend()
plt.show()

YGainSamples = {'Results': ['Mean', 'Variance', 'Probability 10%'], 's = 100': [np.mean(ygainS[:,0]), np.var(ygainS[:,0]), prob[0]], 's = 1000': [np.mean(ygainS[:,1]), np.var(ygainS[:,1]), prob[1]], 's = 10000': [np.mean(ygainS[:,2]), np.var(ygainS[:,2]), prob[2]], 's = 100000': [np.mean(ygainS[:,3]), np.var(ygainS[:,3]), prob[3]]}
df2 = pd.DataFrame(YGainSamples)
print(df2)
