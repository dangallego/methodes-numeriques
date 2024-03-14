'''
Consider a square potential well of width w with walls of height V:
Using Schroödinger's equation, it can be shown that the allowed energies E of a single quantum particle of mass m trapped in the well are solutions of

\tan{\sqrt{w^2 m E / 2 \hbar^2}} = \sqrt{(V-E)/E}    for even numbered states
\tan{\sqrt{w^2 m E / 2 \hbar^2}} = \sqrt{E/(V-E)}    for odd numbered states

where the states are numbered starting from 0, with the ground state being state 0, the first excited state being state 1, and so forth.


    For an electron (mass 9.1094x10-31 kg) in a well with V = 20eV and w = 1nm, write a Python program to plot the three quantities

    y_1 = \tan{\sqrt{w^2 m E / 2 \hbar^2}}, \,\,\,\, y_2 = \sqrt{V-E\over{E}}, \,\,\,\, y_3 = -\sqrt{E\over{V-E}}

    on the same graph, as a function of E from E = 0 to E = 20eV. From your plot make approximate estimates of the energies of the first six energy levels of the particle.
    Write a second program to calculate the values of the first six energy levels in electron volts to an accuracy of 0.001 eV using binary search.

'''
import numpy as np 
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u

m_e = c.m_e.value
hbar = c.hbar.value
V = 20*1.60218e-19    # 20 eV ~= 3.20435e-18  --- CONVERSION FACTOR is 1ev ~= 1.60218e-19 Joules
w = 1e-9           # 1nm == 1e-9 meters

def conv_joules(eV):
    return eV*1.60218e-19

def tan_funct(E):
    E_joules = conv_joules(E) 
    y1 = np.tan(np.sqrt((w**2*m_e*E_joules)/(2*hbar**2)))
    return y1

E = np.linspace(0,20,1000)
#E_joules = conv_joules(E)

def Y2(E):
    E_joules = conv_joules(E) 
    return np.sqrt((V-E_joules)/E_joules)

def Y3(E):
    E_joules = conv_joules(E) 
    return -1*np.sqrt(E_joules/(V-E_joules))

y1 = tan_funct(E)
y2 = Y2(E)
y3 = Y3(E)


def plot1():
    #plt.scatter(E, y1, s = 5, alpha = 0.5, c ='k') # x values can stay in eV (since they should correspond to same y anyways)
    y1[:-1][np.diff(y1) < 0] = np.NaN # gets rid of lines to infinity from tangent values
    plt.plot(E, y1, color = 'k', lw = 1.4, alpha = 0.9, label = 'y1')
    plt.plot(E, y2, alpha = 0.7, label  = 'y2')
    plt.plot(E, y3, alpha = 0.7, label = 'y3')
    plt.axhline(y = 0, color = 'grey', linestyle = '--') 
    plt.xlabel('E')
    plt.xlim(0,20)
    plt.ylim(-4,4)
    plt.legend()
    plt.show()

plot1()


# FIRST BISECTION POINT
def plot2(funct1, funct2):
    fxprime = funct1 - funct2
    fxprime[:-1][np.diff(fxprime) < 0] = np.NaN
    plt.plot(E, fxprime)
    plt.axhline(y = 0, color = 'grey', linestyle = '--') 
    plt.ylim(-1,1)
    plt.show()

plot2(y1, y2) # from plot we see f(x2) == + and f(x1) == - 

# we see that f(x2) at x = 15 should be positive 
# and f(x1) at x = 14.5 should be negative (in terms of y-value)

'''normally would index fxprime array by E1 and E1 values we are looking at
but since we are looking for decimal values, can't index therefore use interpolation '''

# BUT FIRST lets make a function for to get fxprime singular values given x1 and x2
def transformed_funct(E):
    ''' Returns y value at given E for our transformed function, y1 - y2. '''
    #E = conv_joules(Ex) # why no conversion here?? 
    return tan_funct(E) - Y2(E) 

'''
# accuracy threshold
eps = 0.001 #eV
# give initial pair of points (eyeball guess)
x1 = 14
x2 = 16

fx1 = transformed_funct(x1) # I am aware this is sloppy
fx2 = transformed_funct(x2)
print("f(x1)=" , f'{fx1:.3}', " and f(x2)=", f'{fx2:.3}') # voila! Good guess (-0.49 and +0.57) # voila! Good guess (-0.49 and +0.57)
'''

# calculate the midpoint
def midpoint(x1, x2):
    xprime = 0.5*(x1+x2)
    return xprime

def sign(num):
    ''' 
    Function used to check if sign of values are the same. 
    If two variables have the same sign then 
            sign(num1) = 1 == 1 = sign(num2) 
    '''
    if num > 0: 
        return 1
    if num < 0:
        return -1
    else: 
        return 0 


def bisection_method(x1, x2, target_accuracy = 0.001):
    ''' Performs bisection method on a number two points, x1 and x2, which should be opposite in sign. 

    ===========================================================================================
    Uses following steps: 
    1. Given initial pair of points x1 and x2, checks that f(x1) and f(x2) have opposite signs. 
        Accuracy is also defined as argument into function. 
    2. Calculates the midpoint  x' = 1/2 (x1+x2) and evaluates f(x')
    3. If f(x') has same sign as f(x1), sets x1 == x'. 
        Otherwise sets x2 == x' .
    4. If |x1 - x2| > target_accuracy, repeat starting from step 2. 
        Otherwise calculates midpoint final time and returns final estimate of position of root  
     '''
    eps = target_accuracy
    fx1 = transformed_funct(x1) # evaluating x1 value on transformed function (y1 - y2)
    fx2 = transformed_funct(x2) # evaluating x2 value on transformed function 

    # immediate sign check -- fx1 should not be same sign as fx2
    if sign(fx1) == sign(fx2):
        warning = " f(x1) and f(x2) have the same sign. Try again "
        return warning

    while np.abs(x1 - x2) > eps: 
        # step 2: calculate midpoint
        xprime = midpoint(x1,x2)
        # evalue f(x')
        fxprime = transformed_funct(xprime)
        # step 3: check if f(x') sign has changed 
        if sign(fxprime) == sign(fx1): 
            x1 = xprime         # sign of f(x') == f(x1), set x1 = x'
        else: 
            x2 = xprime         # sign of f(x') != f(x1), set x2 = x'
    else:                  
        xprime = midpoint(x1,x2)        # calculates final midpoint estimate
    return xprime

eps = 0.001
x1 = 14
x2 = 16

print(bisection_method(14,16))