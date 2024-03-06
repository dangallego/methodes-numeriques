''' 5.21 DIFFERENTIATION - E field of charge distribution
Suppose we have a distribution of charges and we want to calculate the resulting electric field. 
One way to do this is to first calculate the electric potential φ and then take its gradient. 
For a point charge q at the origin, the electric potential at a distance r from the origin is ___
and the electric field is E = −∇φ.

    1. You have two charges, of ±1 C, 10 cm apart. Calculate the resulting electric potential on a 1 m x 1 m square plane surrounding the charges and passing through them. 
    Calculate the potential at 1 cm spaced points in a grid and make a visualization on the screen of the potential using a density plot.
    
    2. Now calculate the partial derivatives of the potential with respect to x and y and hence find the electric field in the xy plane. 
    Make a visualization of the field also. This is a little trickier than visualizing the potential, because the electric field has both magnitude and direction. 
    A visualization might use the arrow object from the visual package, drawing a grid of arrows with direction and length chosen to represent the field.
'''
import numpy as np 
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u
import CalcFunctions as cf

# let's place q1 at (x=45, y=50) and q2 at (x=55, y=50) 

# Electric Potential function 
def e_potential(r_x, r_y):
    ''' Calculates electric potential for already fixed q1 and q2. '''
    # constants --> using only values for each 
    e = c.e.value      # c.e is electron charge
    e0 = c.eps0.value  # vacuum permittivity

    q1 = 1*e ; q2 = -1*e    # charges
    # q1: (x=45, y=50) and q2: (x=55, y=50) 
    r1 = np.sqrt((45 - r_x)**2 + (50-r_y)**2) # for distance from reference point to q1
    r2 = np.sqrt((55-r_x)**2 + (50-r_y)**2)   # for distance to q2

    potential =(1/(4*np.pi*e0)) * (q1/r1 + q2/r2)
    return potential

# first difficult part - need to iterate over varying x and y values 

N = 100
xarr = np.arange(0,N)
yarr = np.arange(0,N)

xi = np.zeros(N)
yi = np.zeros(N)

potentials = []

# need to keep indices so we know what x and y values the potential corresponds to 
for x in xarr: 
    for y in yarr: 
        # calculate electric potential from point to q1 and q2 for each possible combination
        P = e_potential(x,y)
        # save results 
        potentials.append((x, y, P))

potentials = np.array(potentials)
#for result in potentials:
#    print("X:", potentials[0], "Y:", potentials[1], "Potential", potentials[2])
        
print(potentials.shape)

# Extract x, y, and electric potential values from the array
x_values = potentials[:, 0]
y_values = potentials[:, 1]
e_potential_values = potentials[:, 2]

# Reshape the electric potential values to match the grid - assumes x and y values are evenly spaced
grid_size = int(np.sqrt(len(x_values)))  # must be square grid 
electric_potential_grid = e_potential_values.reshape((grid_size, grid_size))

# Create x and y coordinate grid
x_grid, y_grid = np.meshgrid(np.linspace(min(x_values), max(x_values), grid_size),
                              np.linspace(min(y_values), max(y_values), grid_size))

# Set x and y limits (points in the middle so makes no sense to plot the whole grid)
plt.xlim(40, 60)  # Replace x_min and x_max with your desired limits
plt.ylim(40, 60)  # Replace y_min and y_max with your desired limits

plt.contour(x_grid, y_grid, electric_potential_grid, levels = 15)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Contour Plot of Electric Potential')
plt.colorbar(label='Electric Potential')


#plt.contour(potentials[:,0], potentials[:,1], potentials[:,2])
plt.show()