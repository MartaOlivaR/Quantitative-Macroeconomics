# QUANTITATIVE MACROECONOMICS: PROBLEM SET 3
# Marta Oliva Riera
import pandas as pd
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


#%% ===========================================================================
# QUESTION 1: COMPUTING TRANSITIONS IN A REPRESENTATIVE AGENT ECONOMY
# =============================================================================
# (a). Steady state, where k/y=4 and i/y=0.25:
#------------------------------------------------------------------------------
# I analytically wrote the maximization problem as a lagrangean, took FOCs and obtained some conditions in SS.
# Assuming the following:
theta = 0.67
h = 0.31
y1 = 1       # this is a normalization

i1 = 0.25   # resulting from the normalization and the assumed ratios in the question
k1 = 4

# Conditions to obtain the rest of steady state values:
delta = i1/k1
c1 = y1- delta*k1
z1 = (y1/(k1**(1-theta)*h**theta))**(1/theta)
beta = 1/((1-theta)*((h*z1)/k1)**theta + 1 - delta)

SS1 = {'Variables': ['theta', 'beta', 'delta', 'c', 'i', 'h', 'y', 'z', 'k'], 'Steady State 1': [theta, beta, delta, c1, i1, h, y1, z1, k1]}
df1 = pd.DataFrame(SS1)
print(df1)


# (b) Double z and compute the new steady state:
# -----------------------------------------------------------------------------
# Keeping the same parameters (delta, theta, beta) as well as the same labour and doubling z,
# these are the new steady state values:
z2 = 2*z1
k2 = ((((h*z2)**(theta))/(((1/beta) - (1 - delta))/(1-theta))))**(1/theta)
i2 = delta*k2
y2 = k2**(1-theta) * (z2*h)**theta
c2 = y2 - delta*k2

SS12 = {'Variables': ['theta', 'beta', 'delta', 'c', 'i', 'h', 'y', 'z', 'k'], 'Steady State 1': [theta, beta, delta, c1, i1, h, y1, z1, k1], 'Steady State 2': [theta, beta, delta, c2, i2, h, y2, z2, k2]}
df12 = pd.DataFrame(SS12)
print(df12)

# (c) Transiton between the steady states: 
# -----------------------------------------------------------------------------
c = 1.02443     # initial guess for c0, the immediate jump when the shock in z hits - for now updated manually by me,
                # but should add an iteration algorithm which did this itself!
t = 0
ct = []         # setting up the different variables as arrays and determining their initial values
ct.append(c)
it = []
it.append(i1)
ht = []
ht.append(h)
yt = []
yt.append(y1)
kt = []
kt.append(k1)

# Loop creating the transition path for each of the variables of interest: 
while k2 - kt[t] > 0.1:
    kt.append(kt[t]**(1-theta)*(z2*h)**theta + (1-delta)*kt[t]- ct[t])
    ct.append(((1-theta)*(h*z2)**theta*kt[t+1]**(-theta) + (1 - delta))*ct[t]*beta)
    yt.append(kt[t+1]**(1-theta) * (z2*h)**theta)
    it.append(yt[t+1] - ct[t+1])
    t += 1

# Plotting these time paths: 
plt.figure()
plt.subplot(221)
plt.plot(ct)
plt.title('Consumption transition path')
plt.ylabel('Consumption')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(kt)
plt.title('Capital transition path')
plt.ylabel('Capital')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(yt)
plt.title('Output transition path')
plt.ylabel('Output')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(it)
plt.title('Investment transition path')
plt.ylabel('Investments')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show()

print('Last capital value in the time path: '+ str(kt[-1]))
print('Last consumption value in the time path: '+ str(ct[-1]))


# (d) Unexpected shocks: 
# -----------------------------------------------------------------------------
# Resetting the previously defined initial values:
c = 1.02443     #0.938  #again, this is manually chosen for the transition to converge, which isn't ideal
t = 0
ct = []
ct.append(c)
it = []
it.append(i1)
ht = []
ht.append(h)
yt = []
yt.append(y1)
kt = []
kt.append(k1)

# Loop of the transition now includes an unexpected shock at t=10: 
while (k2 - kt[t] > 0.1) and (t < 10):
    kt.append(kt[t]**(1-theta)*(z2*h)**theta + (1-delta)*kt[t]- ct[t])
    ct.append(((1-theta)*(h*z2)**theta*kt[t+1]**(-theta) + (1 - delta))*ct[t]*beta)
    yt.append(kt[t+1]**(1-theta) * (z2*h)**theta)
    it.append(yt[t+1] - ct[t+1])
    t += 1

while (kt[t] - k1 > 0.1) and (t >= 10):
    if t == 10: 
        ct.append(0.93)     # guessing again what the reaction at the shock is
        kt.append(kt[t]**(1-theta)*(z1*h)**theta + (1-delta)*kt[t]- ct[t])
        yt.append(kt[t+1]**(1-theta) * (z1*h)**theta)
        it.append(yt[t+1] - ct[t+1])
        t += 1
    else: 
        kt.append(kt[t]**(1-theta)*(z1*h)**theta + (1-delta)*kt[t]- ct[t])
        ct.append(((1-theta)*(h*z1)**theta*kt[t+1]**(-theta) + (1 - delta))*ct[t]*beta)
        yt.append(kt[t+1]**(1-theta) * (z1*h)**theta)
        it.append(yt[t+1] - ct[t+1])
        t += 1

# Plotting these new time paths: 
plt.figure()
plt.subplot(221)
plt.plot(ct)
plt.title('Consumption transition path')
plt.ylabel('Consumption')
plt.xlabel('Time')
plt.subplot(222)
plt.plot(kt)
plt.title('Capital transition path')
plt.ylabel('Capital')
plt.xlabel('Time')
plt.subplot(223)
plt.plot(yt)
plt.title('Output transition path')
plt.ylabel('Output')
plt.xlabel('Time')
plt.subplot(224)
plt.plot(it)
plt.title('Investment transition path')
plt.ylabel('Investments')
plt.xlabel('Time')
plt.subplots_adjust(top=2, bottom=0.08, left=0, right=2, hspace=0.3, wspace=0.2)
plt.show()

print('Last capital value in the time path: '+ str(kt[-1]))
print('Last consumption value in the time path: '+ str(ct[-1]))



#%%  ============================================================================================
# QUESTION 2: MULTICOUNTRY MODEL WITH FREE MOBILITY OF CAPITAL AND PROGRESSIVE LABOUR INCOME TAX
# ===============================================================================================
## (2.1) Closed Economy: 
# -----------------------------------------------------------------------------
# Defining the parameter values:
kappa = 5.0
nu = 1.0
sigma = 0.8
etaAH = 5.5     # switched the H and L provided
etaAL = 0.5
etaBH = 3.5
etaBL = 2.5
Z = 1.0
theta = 0.6
kbar = 2.0
lambdaA = 0.95 
lambdaB = 0.84
phi = 0.2

kAH, kAL = 1.0, 1.0   # making an assumption on how the capital endowment of the countries is distributed
kBH, kBL = 1.0, 1.0   
Kd = kbar


'''Notation: 
    r = x0
    w = x1
    hL = x2
    hH = x3
    cL = x4
    cH = x5'''

# Solved for the equilibrium analytically and obtained 6 equations which we need to solve, which
# I put together into the following function to solve as a system of equations:
# Country A:
def f(x): 
    f1=(-x[0]+(1-theta)*Z*(x[2]*etaAL+x[3]*etaAH)**(theta)*Kd**(-theta))        # Firm's FOC wrt K
    f2=(-x[1]+(theta)*Z*Kd**(1-theta)*(x[2]*etaAL+x[3]*etaAH)**(theta-1))       # Firm's FOC wrt H
    f3=(-x[2]+((1-phi)*lambdaA*(1/kappa)*x[4]**(-1*sigma)*(x[1]*etaAL)**(1-phi))**(nu/(1+nu*phi)))      # HH's FOCs, productivity = L
    f4=(-x[3]+((1-phi)*lambdaA*(1/kappa)*x[5]**(-1*sigma)*(x[1]*etaAH)**(1-phi))**(nu/(1+nu*phi)))      # HH's FOCs, productivity = H
    f5=(-x[4]+lambdaA*(x[1]*x[2]*etaAL)**(1-phi)+x[0]*kAL**etaAL)               # Budget constraint, productivity = L
    f6=(-x[5]+lambdaA*(x[1]*x[3]*etaAH)**(1-phi)+x[0]*kAH**etaAH)               # Budget constraint, productivity = H
    return (f1, f2, f3, f4, f5, f6)
solA = fsolve(f,[1,1,1,1,1,1])

# Country B: 
def g(x): 
    g1=(-x[0]+(1-theta)*Z*(x[2]*etaBL+x[3]*etaBH)**(theta)*Kd**(-theta))
    g2=(-x[1]+(theta)*Z*Kd**(1-theta)*(x[2]*etaBL+x[3]*etaBH)**(theta-1))                               
    g3=(-x[2]+((1-phi)*lambdaB*(1/kappa)*x[4]**(-1*sigma)*(x[1]*etaBL)**(1-phi))**(nu/(1+nu*phi)))
    g4=(-x[3]+((1-phi)*lambdaB*(1/kappa)*x[5]**(-1*sigma)*(x[1]*etaBH)**(1-phi))**(nu/(1+nu*phi)))      
    g5=(-x[4]+lambdaB*(x[1]*x[2]*etaBL)**(1-phi)+x[0]*kBL**etaBL)                                       
    g6=(-x[5]+lambdaB*(x[1]*x[3]*etaBH)**(1-phi)+x[0]*kBH**etaBH)                                       
    return (g1, g2, g3, g4, g5, g6)
solB = fsolve(g,[1,1,1,1,1,1])

ab = {'Variables':['r', 'w', 'hL', 'hH', 'cL', 'cH'], 'Country A':[solA[0], solA[1], solA[2], solA[3], solA[4], solA[5]], 'Country B':[solB[0], solB[1], solB[2], solB[3], solB[4], solB[5]]}
dfab = pd.DataFrame(ab)
print(dfab)

## (2.2) Union Economy: 
# -----------------------------------------------------------------------------
# We now have 16 unkowns:
'''Notation:
    rA = x0
    rB = x1
    wA = x2
    wB = x3
    hAL = x4
    hAH = x5
    hBL = x6
    hBH = x7 
    cAL = x8
    cAH = x9
    cBL = x10
    cBH = x11
    kAL = x12
    kAH = x13
    kBL = x14
    kBH = x15'''

# Assumptions on the distribution of capital endowments:
kAH, kAL = 1.0, 1.0
kBH, kBL = 1.0, 1.0
Kd = kbar

# Having solved for the equilibrium conditions analytically, I have 16 conditions, which I write below in 
# a function equivalent to the one I wrote in the previous question, and solve it with the same procedure.
# The only difference is that now both countries are solved at the same time:
def h(x): 
    h1=-x[0]+(1-theta)*Z*(x[5]*etaAH+x[4]*etaAL)**(theta)*Kd**(-1*theta)    # Firm's FOCs
    h2=-x[1]+(1-theta)*Z*(x[6]*etaBH+x[7]*etaBL)**(theta)*Kd**(-1*theta)
    h3=-x[2]+(theta)*Z*Kd**(1-theta)*(x[5]*etaAH+x[4]*etaAL)**(theta-1)
    h4=-x[3]+(theta)*Z*Kd**(1-theta)*(x[6]*etaBH+x[7]*etaBL)**(theta-1)
    h5=-x[4]+((1-phi)*lambdaA*(1/kappa)*x[8]**(-1*sigma)*(x[2]*etaAL)**(1-phi))**(nu/(1+nu*phi))    # HH's FOC wrt h
    h6=-x[5]+((1-phi)*lambdaA*(1/kappa)*x[9]**(-1*sigma)*(x[2]*etaAH)**(1-phi))**(nu/(1+nu*phi))
    h7=-x[6]+((1-phi)*lambdaB*(1/kappa)*x[10]**(-1*sigma)*(x[3]*etaBL)**(1-phi))**(nu/(1+nu*phi))
    h8=-x[7]+((1-phi)*lambdaB*(1/kappa)*x[11]**(-1*sigma)*(x[3]*etaBH)**(1-phi))**(nu/(1+nu*phi))
    h9=-x[8]+lambdaA*(x[2]*x[4]*etaAL)**(1-phi)+x[0]*kAL**etaAL+x[1]*(Kd-kAL)       # Budget constraints
    h10=-x[9]+lambdaA*(x[2]*x[5]*etaAH)**(1-phi)+x[0]*kAH**etaAH+x[1]*(Kd-kAH)
    h11=-x[10]+lambdaB*(x[3]*x[6]*etaBL)**(1-phi)+x[1]*kBL**etaBL+x[0]*(Kd-kBL)
    h12=-x[11]+lambdaB*(x[3]*x[7]*etaBH)**(1-phi)+x[1]*kBH**etaBH+x[0]*(Kd-kBH)
    h13=-x[12]+(x[1]/(x[0]*etaAL))**(1/(etaAL-1))     # HH's FOC wrt k
    h14=-x[13]+(x[1]/(x[0]*etaAH))**(1/(etaAH-1))
    h15=-x[14]+(x[0]/(x[1]*etaBL))**(1/(etaBL-1))
    h16=-x[15]+(x[0]/(x[1]*etaBH))**(1/(etaBH-1))
    return(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15, h16)
solAB = fsolve(h, [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
AB = {'Union':['rA', 'rB', 'wA', 'wB', 'hAL', 'hAH', 'hBL', 'hBH', 'cAL', 'cAH', 'cBL', 'cBH', 'kAL', 'kAH', 'kBL', 'kBH'], 'Solutions':[solAB[0], solAB[1], solAB[2], solAB[3], solAB[4], solAB[5], solAB[6], solAB[7], solAB[8], solAB[9], solAB[10], solAB[11], solAB[12], solAB[13], solAB[14], solAB[15]]}
dfAB = pd.DataFrame(AB)
print(dfAB)