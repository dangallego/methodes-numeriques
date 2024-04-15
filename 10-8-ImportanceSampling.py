'''
Exercise 10.8: Calculate a value for the integral. 

using the importance sampling formula, Eq. (10.42), with w(x) = x^-1/2, as follows.
a) Show that the probability distribution p(x) from which the sample points should be drawn is given by
    p(x) = 1/(2*sqrt(x))
and derive a transformation formula for generating random numbers between zero and one from this distribution.
b) Using your formula, sample N = 1,000,000 random points and hence evaluate the integral. You should get a value around 0.84.
'''
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 

'''
Using transformation function, we get the integral I = 1/N sum(2/e^x +1 ). 
Derived from z = integral of p(x) = 1/2 x^{1/2}
'''

rng = np.random.default_rng()

N = 10_000
z = rng.random(N)
x = 0.5*z**0.5

def funct(x):
    return 2 / (np.exp(x) + 1 )

I = np.sum(funct(x)) / N

message = "Using our probability distribution and transformation function, we evaluate \
            this integral to be (expected ~0.84): "

print(message,I)