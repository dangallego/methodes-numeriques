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
    return x, y

x,y = Galilean_spiral(10000)
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
    return x, y

x,y = Feys_funct(1000)
#plot(x,y, 0.8)

# Create subplot example 
N = 1000
x1,y1 = deltoid_curve(N)
x2,y2 = Galilean_spiral(N)
x3,y3 = Feys_funct(N)


def plot_3(N=1000):
    '''Plots all three plots next to each other. '''
    x1,y1 = deltoid_curve(N)
    x2,y2 = Galilean_spiral(N)
    x3,y3 = Feys_funct(N)
    # makes the plot 
    figure, axs = plt.subplots(1,3, figsize= (14,4))
    axs[0].plot(x1,y1, color = 'darkslategrey', lw = 5)
    axs[0].set_title('Deltoid curve')
    axs[0].set_ylabel('y')
    axs[1].plot(x2,y2, color = 'darkslategrey', lw = 5)
    axs[1].set_title('Galilean spiral')
    axs[1].set_xlabel('x')
    axs[2].plot(x3,y3, color = 'darkslategrey', lw = 2)
    axs[2].set_title('Fey\'s function')
    plt.show()

#plot_3()

# Parser section
parser = argparse.ArgumentParser(
                    prog='Curve-Plotting',
                    description='Shows plot(s) of Polar function in Cartesian coordinates.',
                    epilog='Text at the bottom of help')
# easiest to make -graph optional argument, default is to all 
parser.add_argument('-g', '--graph', choices = ['deltoid', 'galilean', 'feys', 'all', 
                                         'Deltoid', 'Galilean', 'Feys', 'All'], default = 'all', # default is to plot all graphs
                    action = 'store', help='Must select one of the possible graph options.') # added helpful help message
parser.add_argument('-n', '--ntheta', type = int, # optional argument 
                    action = 'store', default = 1000, # default of 1000 points to generate curves
                    help='Optional argument, choose number of points to plot function') 

args = parser.parse_args()
print(args)

if args.graph == 'deltoid' or args.graph =='Deltoid':
    x,y = deltoid_curve(args.ntheta)
    plot(x,y)
    plt.title('Deltoid curve')
elif args.graph == 'galilean' or args.graph == 'Galilean':
    x,y = Galilean_spiral(args.ntheta)
    plot(x,y)
    plt.title('Galilean spiral')
elif args.graph == 'feys' or args.graph == 'Feys':
    x,y = Feys_funct(args.ntheta)
    plot(x,y)
    plt.title('Fey\'s function')
elif args.graph == 'all' or args.graph == 'All':
    plot_3(args.ntheta)
else:
    plot_3(args.ntheta)

