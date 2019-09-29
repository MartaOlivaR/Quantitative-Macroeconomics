# QUANTITATIVE MACROECONOMICS: PROBLEM SET 2
# Marta Oliva Riera

#%% ===========================================================================
# QUESTION 1: UNIVARIATE FUNCTION APPROXIMATION
# =============================================================================
## 1. Approximate f(X)=x^0.321
# -----------------------------------------------------------------------------
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt

# Defining the function of interest:
x = sy.Symbol('x')
def f(x):
    return x**0.321

# Create the necessary functions:
def fact(n):
    """returns the factorial of n"""
    if n <= 0:
        return 1
    else:
        return n*fact(n-1)

def taylor(funct, x0, n):
    """taylor series expansion of function func at point x0, of order n"""
    i = 0
    t = 0
    while i <= n:
        t = t + (funct(x).diff(x, i).subs(x, x0))/(fact(i))*(x-x0)**i
        i += 1
    return t

# Plot the approximations of different orders:
def plot(F, X, LIM, ORD, LAB):
    """plots taylor approximations of a function F evaluated at X along with its
    true value, in several orders (ORD) along the specified range (LIM).
    The true function is labeled with LAB"""
    x_lims = LIM
    x1 = np.linspace(x_lims[0],x_lims[1],800)       # get rid of magic numbers!!!
    y1 = []
    for j in ORD:     # For loop orders 1, 2, 5 and 20
        func = taylor(F,1,j)
        print('Taylor expansion at n='+str(j),func)
        for k in x1:
            y1.append(func.subs(x,k))
        plt.plot(x1,y1,label='order '+str(j))
        y1 = []
    # Plot the actual function, and display all graphs:
    plt.plot(x1,F(x1),label=LAB)
    plt.xlim(x_lims)
    plt.ylim([-2,7])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title('Taylor series approximation')
    plt.show()  

plot(f, 1, [0,4], [1,2,5,20], '$x^{0.321}$')

#%% ---------------------------------------------------------------------------
# 2. Approximate the ramp function:
# -----------------------------------------------------------------------------
# Define the ramp function:
x = sy.Symbol('x', real=True)
def g(x):
    return (x + abs(x))/2

# Using the same fact, taylor and plot functions, updated for the new case:
plot(g, 2, [-2, 6], [1,2,5,20], '$(x+|x|)/2')

#%% ---------------------------------------------------------------------------
# 3. Approximate three functions with different methods:
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import numpy.polynomial.chebyshev as cheb
import sympy as sy

xlim = [-1,1]
N = 11                                      # number of nodes
xx = np.linspace(xlim[0], xlim[1], 200)     # domain

# (1) Evenly spaced interpolation nodes:
x = np.linspace(xlim[0], xlim[1], N)      # linspace returns evenly spaced points

# (2) Chebyshev interpolation nodes:
xC = np.cos((2 * np.arange(1, N + 1) - 1) * np.pi / (2 * N))

## EXP FUNCTION:
# =============================================================================
def h(x):
    h = np.exp(1/x)
    return h

 # Evenly spaced nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(x, h(x), 3)   # polyfit calculates the coefficients of interpolation
fitted3 = poly.polyval(xx, theta3)  # polyval evaluates the function with the previously calculated coefficients, it results in the interpolated function
err3 = abs(h(xx) - fitted3)         # error = true function - interpolation
theta5 = poly.polyfit(x, h(x), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(h(xx) - fitted5)
theta10 = poly.polyfit(x, h(x), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(h(xx) - fitted10)

plt.subplot(121)
plt.plot(x, h(x), 'o', label='nodes')       # interpolation nodes
plt.plot(xx, h(xx), label='True')           # true values of the function
plt.plot(xx, fitted3, label='Order 3')      # interpolation of order 3
plt.plot(xx, fitted5, label='Order 5')      # interpolation of order 5
plt.plot(xx, fitted10, label='Order 10')    # interpolation of order 10
plt.ylim([-20000,100000])
plt.title('Evenly spaced nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')         # error of the order 3 interpolation
plt.plot(xx, err5, label='Order 5')         # error of the order 5 interpolation
plt.plot(xx, err10, label='Order 10')       # error of the order 10 interpolation
plt.ylim([-20000,100000])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

 # Chebysheb nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(xC, h(xC), 3)
fitted3 = poly.polyval(xx, theta3)
err3 = abs(h(xx) - fitted3)
theta5 = poly.polyfit(xC, h(xC), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(h(xx) - fitted5)
theta10 = poly.polyfit(xC, h(xC), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(h(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, h(xC), 'o', label='nodes')
plt.plot(xx, h(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-20000,100000])
plt.title('Chebyshev nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([-20000,100000])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

# Chebyshev nodes, Chebyshev polynomials:
 # ----------------------------------------------------------------------------
theta3 = cheb.chebfit(xC, h(xC), 3)
fitted3 = cheb.chebval(xx, theta3)
err3 = abs(h(xx) - fitted3)
theta5 = cheb.chebfit(xC, h(xC), 5)
fitted5 = cheb.chebval(xx, theta5)
err5 = abs(h(xx) - fitted5)
theta10 = cheb.chebfit(xC, h(xC), 10)
fitted10 = cheb.chebval(xx, theta10)
err10 = abs(h(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, h(xC), 'o', label='nodes')
plt.plot(xx, h(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-20000,100000])
plt.title('Chebyshev nodes, Chebyshev basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([-20000,100000])
plt.title('Errors')
plt.legend(loc='best')
plt.show()



## RUNGE FUNCTION:
# =============================================================================
def i(x):
	return 1.0 / (1.0 + 25.0 * x**2)

 # Evenly spaced nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(x, i(x), 3)
fitted3 = poly.polyval(xx, theta3)
err3 = abs(i(xx) - fitted3)
theta5 = poly.polyfit(x, i(x), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(i(xx) - fitted5)
theta10 = poly.polyfit(x, i(x), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(i(xx) - fitted10)

plt.subplot(121)
plt.plot(x, i(x), 'o', label='nodes')
plt.plot(xx, i(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,2])
plt.title('Evenly spaced nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,2])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

 # Chebysheb nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(xC, i(xC), 3)
fitted3 = poly.polyval(xx, theta3)
err3 = abs(i(xx) - fitted3)
theta5 = poly.polyfit(xC, i(xC), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(i(xx) - fitted5)
theta10 = poly.polyfit(xC, i(xC), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(i(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, i(xC), 'o', label='nodes')
plt.plot(xx, i(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,2])
plt.title('Chebyshev nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,2])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

# Chebyshev nodes, Chebyshev polynomials:
 # ----------------------------------------------------------------------------
theta3 = cheb.chebfit(xC, i(xC), 3)
fitted3 = cheb.chebval(xx, theta3)
err3 = abs(i(xx) - fitted3)
theta5 = cheb.chebfit(xC, i(xC), 5)
fitted5 = cheb.chebval(xx, theta5)
err5 = abs(i(xx) - fitted5)
theta10 = cheb.chebfit(xC, i(xC), 10)
fitted10 = cheb.chebval(xx, theta10)
err10 = abs(i(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, i(xC), 'o', label='nodes')
plt.plot(xx, i(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,2])
plt.title('Chebyshev nodes, Chebyshev basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,2])
plt.title('Errors')
plt.legend(loc='best')
plt.show()



## RAMP FUNCTION:
# =============================================================================
def g(x):
    return (x + abs(x))/2

 # Evenly spaced nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(x, g(x), 3)
fitted3 = poly.polyval(xx, theta3)
err3 = abs(g(xx) - fitted3)
theta5 = poly.polyfit(x, g(x), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(g(xx) - fitted5)
theta10 = poly.polyfit(x, g(x), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(g(xx) - fitted10)

plt.subplot(121)
plt.plot(x, g(x), 'o', label='nodes')
plt.plot(xx, g(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,1])
plt.title('Evenly spaced nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,0.5])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

 # Chebysheb nodes, monomial basis:   
 # ----------------------------------------------------------------------------
theta3 = poly.polyfit(xC, g(xC), 3)
fitted3 = poly.polyval(xx, theta3)
err3 = abs(g(xx) - fitted3)
theta5 = poly.polyfit(xC, g(xC), 5)
fitted5 = poly.polyval(xx, theta5)
err5 = abs(g(xx) - fitted5)
theta10 = poly.polyfit(xC, g(xC), 10)
fitted10 = poly.polyval(xx, theta10)
err10 = abs(g(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, g(xC), 'o', label='nodes')
plt.plot(xx, g(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,1])
plt.title('Chebyshev nodes, monomial basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,0.5])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

# Chebyshev nodes, Chebyshev polynomials:
 # ----------------------------------------------------------------------------
theta3 = cheb.chebfit(xC, g(xC), 3)
fitted3 = cheb.chebval(xx, theta3)
err3 = abs(g(xx) - fitted3)
theta5 = cheb.chebfit(xC, g(xC), 5)
fitted5 = cheb.chebval(xx, theta5)
err5 = abs(g(xx) - fitted5)
theta10 = cheb.chebfit(xC, g(xC), 10)
fitted10 = cheb.chebval(xx, theta10)
err10 = abs(g(xx) - fitted10)

plt.subplot(121)
plt.plot(xC, g(xC), 'o', label='nodes')
plt.plot(xx, g(xx), label='True')
plt.plot(xx, fitted3, label='Order 3')
plt.plot(xx, fitted5, label='Order 5')
plt.plot(xx, fitted10, label='Order 10')
plt.ylim([-0.5,1])
plt.title('Chebyshev nodes, Chebyshev basis')
plt.legend(loc='best')

plt.subplot(122)
plt.plot(xx, err3, label='Order 3')
plt.plot(xx, err5, label='Order 5')
plt.plot(xx, err10, label='Order 10')
plt.ylim([0,0.5])
plt.title('Errors')
plt.legend(loc='best')
plt.show()

    
    # exponential is not right, but can't find what's wrong
    # 2nd & 3rd graphs look very similar
    


#%% ===========================================================================
# QUESTION 2: MULTIVARIATE FUNCTION APPROXIMATION 
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# CES function:
def f(alpha, sigma, k, h):
    return ((1-alpha)*k**((sigma-1)/sigma)+alpha*h**((sigma-1)/sigma))**(sigma/(sigma-1))

# Setting the parameters
alpha = 0.5                 # relative input share parameter
sigma = 0.25                # elasticity of substitution between capital and labour
N = 20                      # interpolation nodes
kd = np.linspace(0, 10, N)   # capital
hd = np.linspace(0, 10, N)   # labour
a,b = min(kd), max(kd)
c,d = min(hd), max(hd)

# (1) Generate Chebyshev nodes:
z = np.cos((2 * np.arange(1, N + 1) - 1)* np.pi / (2 * N))

# (2) Adjust the nodes' intervals to (0,10)x(0,10):
k = []
for i in z:
    kk = (i+1)*((b-a)/2)+a
    k.append(kk)
k = np.asarray(k)

h = []
for i in z:
    hh = (i+1)*((d-c)/2)+c
    h.append(hh)
h = np.asarray(h)

K,H = np.meshgrid(k,h)

# (3) Evaluate the function at the nodes:
w = f(alpha, sigma, K, H)

# (4) Compute the interpolation coefficients using Chebyshev's basis of degreee 3 to 15:
degree = 3

def T(nodes, degree):
    """function that results in Chebyshev's polynomials of a specified degreee
    using the provided interpolation nodes"""
    psi = []
    psi.append(np.ones(N))              # 1st polynomial
    psi.append(nodes)                   # 2nd polynomial
    for i in range(1, degree):          # further polynomials
        p = 2*nodes*psi[i-1]-psi[i-2]
        psi.append(p)
        psi = np.matrix(psi[degree])
    return psi

def coeff(nodes, fnodes, degree):
    """function returning the Chebyshev coefficients associated with the 
    Chebyshev polynomial basis T()"""
    theta=np.empty((degree+1)*(degree+1))
    theta.shape = (degree+1, degree+1)
    for i in range(degree+1):
        for j in range(degree+1):
            theta[i,j] = np.sum(np.array(fnodes)*np.array(np.dot(np.transpose(T(nodes, i)),T(nodes, j))))/np.array((T(nodes, i)*np.transpose(T(nodes, i)))*(T(nodes, j)*np.transpose(T(nodes,j))))
                        #(np.sum(np.array(fnodes)*np.array(((T(i,enodes).T @ T(j,enodes)))))/np.array((T(i,enodes)*T(i,enodes).T)*(T(j,enodes)*T(j,enodes).T)))
    return theta

def fapprox(x, y, theta, degree):
    """function resulting in the approximation of function f(x, y) for values of
    x in the interval [a,b] and y in the interval [c,d]"""
    f = []
    Zi = (2*(x-a)/(b-a))-1
    Zj = (2*(y-c)/(d-c))-1
    for i in range(degree):
        for j in range(degree):
            f.append(np.array(theta[i,j]*np.dot(np.transpose(T(Zi, i))*T(Zj, j))))
    f = sum(f)
    return f

# Run the approximation of degree 3 and plot it, along with its errors:
theta3 = coeff(z, w, 3)
approx3 = fapprox(k, h, theta3, 3)
true = f(alpha, sigma, k, h)
error3 = abs(true - approx3)

fig = plt.figure()
fig.suptitle('σ=0.25, Order 3')
# true function
ax = fig.add_subplot(131, projection='3d')
ax.set_title("True function")
K,H = np.meshgrid(k,h)
ax.set_xlabel('Capital')
ax.set_ylabel('Labour')
ax.set_zlabel('Ouput')
ax.plot_surface(K, H, true)
# approximation
ax = fig.add_subplot(132, projection='3d')
ax.set_title("Approximated function")
K,H = np.meshgrid(k,h)
ax.set_xlabel('Capital')
ax.set_ylabel('Labour')
ax.set_zlabel('Ouput')
ax.plot_surface(K, H, approx3)
# errors
ax = fig.add_subplot(133, projection='3d')
ax.set_title("Errors")
K,H = np.meshgrid(k,h)
ax.set_xlabel('Capital')
ax.set_ylabel('Labour')
ax.set_zlabel('Ouput')
ax.plot_surface(K, H, error3)

plt.show()

