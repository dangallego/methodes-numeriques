'''Exercise 4.3: Calculating derivatives

Suppose we have a function f(x) and we want to calculate its derivative at a point x. 
We can do that with pencil and paper if we know the mathematical form of the function, 
or we can do it on the computer by making use of the definition of the derivative:

    df/dx = lim (del->0) f(x+del) - f(x) / del 

On the computer we can't actually take the limit as δ goes to zero, but we can get a reasonable approximation just by making δ small.

    a)  Write a program that defines a function f(x) returning the value x(x-1), then calculates the derivative of the function at the point x = 1 
    using the formula above with δ = 10-2. Calculate the true value of the same derivative analytically and compare with the answer your program gives. 
    The two will not agree perfectly. Why not?
    b)  Repeat the calculation for δ = 10-4, 10-6, 10-8, 10-10, 10-12, and 10-14. You should see that the accuracy of the calculation initially gets better 
    as δ gets smaller, but then gets worse again. Plot your result as a function of δ. Why is this?
'''
# part a
import numpy as np
import matplotlib.pyplot as plt
import time 

start = time.time()

def funct(x):
    return x*(x-1)

def derivative(x=1, delta = 1e-2):
    return (funct(x+delta) - funct(x)) / delta

def derivatives(delta, x=1):
    '''Here delta should be array '''
    return (funct(x+delta) - funct(x)) / delta # returns multiple values

#print(derivative())

#part b 
delta_arr = np.array((1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14))
derivs = derivatives(delta_arr)

end = time.time()

print(derivs, "\n")
print(f"Program took :{end-start:.2e} seconds to run", "\n")

#plot 
#plt.scatter(derivs, delta_arr) ; plt.axvline(x = 1, c='k', label='Actual Derivative') ; plt.gca().invert_yaxis() ; plt.ylabel('$\delta$') ; plt.xlabel('Calculated Derivative')
#plt.scatter(delta_arr, derivs) ; plt.axhline(y=1,c='k', label='Actual Derivative') ; plt.gca().invert_xaxis() ; plt.xlabel('$\delta$') ; plt.ylabel('Calculated Derivative')

#LOG 
plt.scatter(np.log(delta_arr), derivs) ; plt.axhline(y=1,c='k', label='Actual Derivative') ; plt.gca().invert_xaxis() ; plt.xlabel('log($\delta$)') ; plt.ylabel('Calculated Derivative')
plt.legend()
plt.show()