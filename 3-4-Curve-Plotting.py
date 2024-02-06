## EXERCISE 3.2 CURVE PLOTTING 
'''Although the plot function is designed primarily for plotting standard xy graphs, it can be adapted for other kinds of plotting as well.

    a)  Make a plot of the so-called deltoid curve, which is defined parametrically by the equations, x = 2 cos θ + cos 2θ, y = 2 sin θ - sin 2θ, where 0 ≤ θ < 2π. 
    Take a set of values of θ between zero and 2π and calculate x and y for each from the equations above, then plot y as a function of x.
    b)  Taking this approach a step further, one can make a polar plot r = f(θ) for some function f by calculating r for a range of values of θ and then converting 
    r and θ to Cartesian coordinates using the standard equations x = r cos θ, y = r sin θ. Use this method to make a plot of the Galilean spiral, r=θ2 for 0 ≤ θ ≤ 10π.
    c)  Using the same method, make a polar plot of “Fey's function”
'''
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot(x,y, alpha = 1):
    '''Plots x and y values.'''
    plt.plot(x,y, color = 'k')
    plt.xlabel('x') ; plt.ylabel('y')
    plt.show()

def polarplot(r,theta):
    '''Polar representation of values.'''
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta, r)
    ax.grid(True)
    plt.show()

def deltoid_curve(num_theta_vals=100):
    '''Takes a set of values of theta between zero and 2pi and calculates / plots x and y.''' 
    theta_vals = np.linspace(0, 2*np.pi, num_theta_vals) 
    sin = np.sin ; cos = np.cos
    x = 2*cos(theta_vals) + cos(2*theta_vals) ; y = 2*sin(theta_vals) - sin(2*theta_vals)
    return x, y

x,y = deltoid_curve(1000)
#plot(x,y)

def Galilean_spiral(num_theta_vals=100): 
    '''Takes x and y values to compute convert from Cartesian to Polar coordinates. 
        Plots Galilean spiral '''
    theta_vals = np.linspace(0, 10*np.pi, num_theta_vals)  # 0 <= theta <= 10pi
    sin = np.sin ; cos = np.cos 
    # polar coordinates
    r = theta_vals**2 
    x = r*cos(theta_vals) ; y = r*sin(theta_vals)
    return x, y, r, theta_vals

x,y,r,thetas = Galilean_spiral(10000)
#plot(x,y)

# POLAR REPRESENTATION
#polarplot(r,thetas)

def Feys_funct(num_theta_vals=100): 
    '''Uses same method as above to make a polar plot of Fey's function.
    r = e^(cos(θ)) - 2*cos4(θ) + sin^5(θ/12)'''
    sin = np.sin ; cos = np.cos 
    theta_vals = np.linspace(0, 24*np.pi, num_theta_vals)
    # polar coordinates
    r = np.exp(cos(theta_vals)) - 2*cos(4*theta_vals) + (sin(theta_vals/12))**5
    x = r*cos(theta_vals) ; y = r*sin(theta_vals)
    return x, y, r, theta_vals

x,y,r,thetas = Feys_funct(1000)
plot(x,y, 0.8)