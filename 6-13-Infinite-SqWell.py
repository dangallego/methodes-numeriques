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
    ''' Converts from eV to joules. '''
    return eV*1.60218e-19

def tan_funct(E):
    ''' Calculates y1 using tangent function. '''
    E_joules = conv_joules(E) 
    y1 = np.tan(np.sqrt((w**2*m_e*E_joules)/(2*hbar**2)))
    return y1

E = np.linspace(0,20,1000) # creates array of Energy values from 0 to 20 
#E_joules = conv_joules(E)

def Y2(E):
    ''' Calculates y2 quantity. '''
    E_joules = conv_joules(E) 
    return np.sqrt((V-E_joules)/E_joules)

def Y3(E):
    ''' Calculates y3 quantity. '''
    E_joules = conv_joules(E) 
    return -1*np.sqrt(E_joules/(V-E_joules))

y1 = tan_funct(E)
y2 = Y2(E)
y3 = Y3(E)


def plot1():
    ''' Plots all three quantites on one plot, fixing limit and boundary issues that arise from plotting tangent functions. '''
    #plt.scatter(E, y1, s = 5, alpha = 0.5, c ='k') # x values can stay in eV (since they should correspond to same y anyways)
    y1[:-1][np.diff(y1) < 0] = np.NaN # gets rid of lines to infinity from tangent values
    plt.plot(E, y1, color = 'k', lw = 1.4, alpha = 0.9, label = 'y1')
    plt.plot(E, y2, alpha = 0.7, label  = 'y2')
    plt.plot(E, y3, alpha = 0.7, label = 'y3')
    plt.axhline(y = 0, color = 'grey', linestyle = '--') 
    plt.xlabel('E (eV)')
    plt.xlim(0,20)
    plt.ylim(-4,4)
    plt.legend()
    plt.show()


def plot2(funct1, funct2):
    ''' Plots a transformed function. '''
    fxprime = funct1 - funct2       # not the same fxprime as in our bisection calculation (just a new function)
    fxprime[:-1][np.diff(fxprime) < 0] = np.NaN
    plt.plot(E, fxprime, c = 'steelblue', alpha = 0.9)
    plt.axhline(y = 0, color = 'grey', linestyle = '--', alpha = 0.7) 
    plt.xlabel('E (eV)')
    plt.xlim(0,20)
    plt.ylim(-1,1)
    plt.show()

# first plot to visualize x1 and x2 points
#plot2(y1, y2) # we see that f(x2) at x = 15 should be positive and f(x1) at x = 14.5 should be negative (in terms of y-value)
    
'''normally would index fxprime array by E1 and E1 values we are looking at
but since we are looking for decimal values, can't index therefore use interpolation '''

# first define function to get fxprime singular values given x1 and x2
def transformed_funct(E):
    ''' Returns y value at given E for our transformed function, y1 - y2. '''
    return tan_funct(E) - Y2(E) 

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


def bisection_method(x1, x2, funct, target_accuracy = 0.001):
    ''' Performs bisection method on a number two points, x1 and x2, which should be opposite in sign. 
    -------------------------------------------------------------------------------------------
    x1 -- float; value used to evaulate the given function at this point 
    x2 -- float; value used to evaluate the given function at this point, should be opposite in sign of x1
    funct -- function that should be predefined and called to evaluate points x1 and x2
    target_accuray -- float; the desired accuracy for performing binary search / bisection method. 
    ===========================================================================================
    Performs following steps: 
    1. Given initial pair of points x1 and x2, checks that f(x1) and f(x2) have opposite signs. 
        Returns warning otherwise.
    2. Calculates the midpoint  x' = 1/2 (x1+x2) and evaluates f(x') .
    3. If f(x') has same sign as f(x1), sets x1 == x'. 
        Otherwise sets x2 == x' .
    4. If |x1 - x2| > target_accuracy, repeat starting from step 2. 
        Otherwise calculates midpoint final time and returns final estimate of position of root  
     '''
    eps = target_accuracy
    fx1 = funct(x1) # evaluating x1 value on transformed function (y1 - y2)
    fx2 = funct(x2) # evaluating x2 value on transformed function 

    # immediate sign check -- fx1 should not be same sign as fx2
    if sign(fx1) == sign(fx2):
        warning = " f(x1) and f(x2) have the same sign. Try again "
        return warning

    while np.abs(x1 - x2) > eps: 
        # step 2: calculate midpoint
        xprime = midpoint(x1,x2)
        # evalue f(x')
        fxprime = funct(xprime)
        # step 3: check if f(x') sign has changed 
        if sign(fxprime) == sign(fx1): 
            x1 = xprime         # sign of f(x') == f(x1), set x1 = x'
        else: 
            x2 = xprime         # sign of f(x') != f(x1), set x2 = x'
    else:                  
        xprime = midpoint(x1,x2)        # calculates final midpoint estimate
    return xprime


def bisection_plot(point1, point2, point3, y_n): 
    ''' 
    Plots roots obtained from bisection function. 
    y_n -- the y_n function being subtracted from y1 (so either y2 or y3 are acceptable arguments). 
    point1,2,3 -- the bisection point (root) calculated.
    '''
    plt.scatter(point1, 0, marker = 'x', c = 'r', label = 'bisection point', s = 50)
    plt.scatter(point2, 0, marker = 'x', c = 'r', s = 50)
    plt.scatter(point3, 0, marker = 'x', c = 'r', s = 50)
    plt.legend()
    plot2(y1, y_n) # y2 or y3
    plt.show()
############################
# For first transformed function find all the bisection midpoints 
eps = 0.001

point1 = bisection_method(2.59,2.95, transformed_funct)
point2 = bisection_method(7.45,8.10, transformed_funct)
point3 = bisection_method(14,16, transformed_funct)

#plt.ylabel('y1 - y2')
#bisection_plot(point1, point2, point3, y2)

# for second transformed function (first need to define function)
def transformed_funct2(E):
    ''' Returns y value at given E for second transformed function, y1 - y3. '''
    return tan_funct(E) - Y3(E) 

# second plot --> to visualize for finding possible x1 and x2 values
#plt.ylabel('y1 - y3') ; plot2(y1, y3)

point4 = bisection_method(0.97, 1.8, transformed_funct2)
point5 = bisection_method(4.6, 5.7, transformed_funct2)
point6 = bisection_method(10.9, 12, transformed_funct2)

#plt.ylabel('y1 - y3')
#bisection_plot(point4, point5, point6, y3)


def plot1():
    ''' Plots all three quantites on one plot, fixing limit and boundary issues that arise from plotting tangent functions. '''
    #plt.scatter(E, y1, s = 5, alpha = 0.5, c ='k') # x values can stay in eV (since they should correspond to same y anyways)
    y1[:-1][np.diff(y1) < 0] = np.NaN # gets rid of lines to infinity from tangent values
    plt.plot(E, y1, color = 'k', lw = 1.4, alpha = 0.9, label = 'y1')
    plt.plot(E, y2, alpha = 0.7, label  = 'y2')
    plt.plot(E, y3, alpha = 0.7, label = 'y3')
    plt.axhline(y = 0, color = 'grey', linestyle = '--') 
    plt.xlabel('E (eV)')
    plt.xlim(0,20)
    plt.ylim(-4,4)
    plt.legend()
    plt.show()


def plot2(funct1, funct2):
    ''' Plots a transformed function. '''
    fxprime = funct1 - funct2       # not the same fxprime as in our bisection calculation (just a new function)
    fxprime[:-1][np.diff(fxprime) < 0] = np.NaN
    plt.plot(E, fxprime, c = 'steelblue', alpha = 0.9)
    plt.axhline(y = 0, color = 'grey', linestyle = '--', alpha = 0.7) 
    plt.xlabel('E (eV)')
    plt.xlim(0,20)
    plt.ylim(-1,1)
    plt.show()



def plotall():
    # Create figure and axes objects
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 5)) # nrows=1 for one row, ncols=3 for three columns

    plt.xlim(0,20)
    plt.ylim(-4,4)
    # Plot data on each subplot
    y1[:-1][np.diff(y1) < 0] = np.NaN
    axs[0].plot(E, y1, color = 'k', lw = 1.4, alpha = 0.9, label = 'y1')
    axs[0].plot(E, y2, alpha = 0.7, label  = 'y2')
    axs[0].plot(E, y3, alpha = 0.7, label = 'y3')
    axs[0].axhline(y = 0, color = 'grey', linestyle = '--') 
    axs[0].set_title('y1, y2, y3')
    axs[0].set_xlim([0, 20])
    axs[0].set_ylim([-4, 4])

    # plot 2
    fxprime = y1 - y2    
    axs[1].set_title('y1 - y2')   
    axs[1].plot(E, fxprime, c = 'steelblue', alpha = 0.9)
    axs[1].axhline(y = 0, color = 'grey', linestyle = '--', alpha = 0.7) 
    axs[1].set_xlabel('E (eV)')
    axs[1].scatter(point1, 0, marker = 'x', c = 'r', label = 'bisection point', s = 50)
    axs[1].scatter(point2, 0, marker = 'x', c = 'r', s = 50)
    axs[1].scatter(point3, 0, marker = 'x', c = 'r', s = 50)
    axs[1].set_xlim([0, 20])
    axs[1].set_ylim([-1, 1])

    #plot 3
    fxprime2 = y1 - y3     
    axs[2].set_title('y1 - y3')  
    axs[2].plot(E, fxprime2)
    axs[2].plot(E, fxprime2, c = 'steelblue', alpha = 0.9)
    axs[2].axhline(y = 0, color = 'grey', linestyle = '--', alpha = 0.7) 
    axs[2].set_title('y1 - y3')
    axs[2].scatter(point4, 0, marker = 'x', c = 'r', label = 'bisection point', s = 50)
    axs[2].scatter(point5, 0, marker = 'x', c = 'r', s = 50)
    axs[2].scatter(point6, 0, marker = 'x', c = 'r', s = 50)
    axs[2].set_xlim([0, 20])
    axs[2].set_ylim([-1, 1])

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    plt.show()

print("Visualization of the first first six energy levels for the functions y1, y2, and y3 (for the case of a square potential well), ",
      "\n", "as well as the calculated first six energy levels with an accuracy of 0.001 (eV) using bisection method. ", "\n", 
      "The first four values are: ", f'{point1:.2}', f'{point2:.2}' , f'{point3:.2}', f'{point4:.2}'  )
plotall()

