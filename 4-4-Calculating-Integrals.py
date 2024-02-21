'''Suppose we want to calculate the value of the integral

I = \int_{-1}^1 \sqrt{1-x^2} \> dx.

The integrand looks like a semicircle of radius 1:
and hence the value of the integral—the area under the curve—must be  ½π = 1.57079632679 . . .

Alternatively, we can evaluate the integral on the computer by dividing the domain of integration into a large number N of slices of 
width h = 2/N each and then using the Riemann definition of the integral:

I = \lim_{N\to\infty} \sum_{k=1}^N hy_k

where
y_k = \sqrt{1-x_k^2} and x_k = -1 + hk.

We cannot in practice take the limit N → ∞, but we can make a reasonable approximation by just making N large.

    a)  Write a program to evaluate the integral above with N = 100 and compare the result with the exact value. The two will not agree very well, because N = 100 is not a sufficiently large number of slices.
    b)  Increase the value of N to get a more accurate value for the integral. If we require that the program runs in about one second or less, how accurate a value can you get?

'''
import numpy as np 
import timeit
import CalcFunctions as cf

def funct(x): # needs to be passed x 
    function = np.sqrt(1-x**2)
    return function

 # Riemann method 
N = 10
print(cf.riemann_integral(funct, a=-1, b=1,N= N))
print(timeit.timeit("cf.riemann_integral(cf.funct4_4,a=-1,b=1,N=10)", "import CalcFunctions as cf")) # for timeit.timeit, arguments have to be defined in string separately from above
# therefore change N in the string above for accurate time 