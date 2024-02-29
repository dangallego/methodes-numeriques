'''Exercise 5.12: The Stefan-Boltzmann constant
The Planck theory of thermal radiation tells us that in the (angular) frequency interval ω to ω + dω,
 a black body of unit area radiates electromagnetically an amount of thermal energy per second equal to I(ω) dω, where 
                    (FORMULA).
Here \hbar is Planck's constant over 2π, c is the speed of light, and kB is Boltzmann's constant.

    a)  Show that the total energy per unit area radiated by a black body is: 
                    (FORMULA)
    b)  Write a program to evaluate the integral in this expression. Explain what method you used, and how accurate you think your answer is.
    c)  Even before Planck gave his theory of thermal radiation around the turn of the 20th century, it was known that the total energy W given off by a black body per unit area per second followed Stefan’s law: W = σT4, where σ is the Stefan–Boltzmann constant. Use your value for the integral above to compute a value for the Stefan–Boltzmann constant (in SI units) to three significant figures. Check your result against the known value, which you can find in books or on-line. You should get good agreement.

'''
import numpy as np 
import astropy.constants as c
import astropy.units as u
import CalcFunctions as cf

#help(c) # cool !

def funct(z):
    '''Function with new limits. 
    tan(z)^3 / e^{tan(z)} - 1
    '''
    num = np.tan(z)**3
    denom = np.exp(np.tan(z)) - 1
    return (num / denom) * (1/np.cos(z))**2 # no secant defined 

# limits for tan cannot be exactly pi/2 or 0, so have them be ever so slightly below/above
A = 0+1E-6 ; B = (np.pi/2) - 1E-6

int = cf.simpsons_int(funct, a=A, b = B, Nmax=10000)

#print(int) # should be 6.49394 (wolfram

# Now find the total energy per unit area radiated by a black body (which is the integrand we just did multiplied by some constants)
W1 = c.k_B**4 / (4*np.pi**2 * c.c**2 * c.hbar**3 )
#help(c)

W = W1 * int
print(W)

print(c.sigma_sb.value)