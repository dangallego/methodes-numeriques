import numpy as np 
import time
import timeit


def default_funct(x): # needs to be passed x 
    ''' Default function for class examples so far.'''
    function = x**4-2*x+1
    return function


def funct4_4(x): 
    '''Function for exercise 4.4'''
    function = np.sqrt(1-x**2)
    return function

 # Riemann method 
def riemann_integral(function, a,b, N):
    '''Integrates using Riemann definition of integral'''
    h = (b-a)/N
    Sum = 0 
    for k in range(1, N):
        Sum += function(a+k*h)
    F = h*(function(a) + function(b) + Sum)
    return F


def trapezoid_int(function, a, b, Nmax):
    '''Integrates a given function from limits a to b,
        using trapzoidal rule.'''
    delx = (b-a)/Nmax # gives delta_x; how small incriments are between limits a and b (often called h)

    sum_term = 0 # need to initialize sum_term variable 
    for k in range(1, Nmax): # from 1 to Nmax
        sum_term = sum_term + function(a+k*delx) # this serves as "summation term"

    F = delx*(0.5*function(a) + 0.5*function(b) + sum_term)
    

def fast_trapezoid_int(function, a,b, Nmax):
    '''Integrates a given function from limits a to b,
        using trapzoidal rule without using loop.'''
    N = np.arange(Nmax)
    delx = (b-a)/Nmax
    sum_term = np.sum(function(a+N*delx))
    F = delx*(0.5*function(a)+0.5*function(b)+sum_term)
    return F


def simpsons_int(function, a, b, Nmax=100):
    '''Integrates same function as before using Simpson's rule'''
    delx = (b-a) / Nmax

    sum1 = 0 
    for k in range(1, Nmax//2): # from 1 to Nmax // 2 (summing from k to N/2) --> N//2 to do integer division because should be integers for summation! 
        sum1 += function(a+(2*k-1)*delx) #sum1 += is same as sum1 = sum1 + ... 
    sum2 = 0 
    for k in range(1, (Nmax//2)-1):
        sum2 += function(a+2*k*delx)
    
    F = (delx)/3 * (function(a) + function(b) +4*sum1 + 2*sum2)
    return F


