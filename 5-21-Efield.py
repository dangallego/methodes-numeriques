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

    potential = (1/(4*np.pi*e0)) * (q1/r1 + q2/r2)
    return potential

# need to iterate over varying x and y values 
N = 100
xarr = np.arange(0,N)
yarr = np.arange(0,N)

# need to keep indices so we know what x and y values the potential corresponds to -- want to keep as (100,100) 2D array 
e_potentials = np.zeros((100, 100)) # intialize 2D array 

for i, x in enumerate(xarr): 
    for j, y in enumerate(yarr):
        # calculate electric potential from point to q1 and q2 for each possible combination 
        e_potentials[i, j] = e_potential(x,y)

#print(e_potentials)


def electric_contour(xarr, yarr, electric_potentials):
    ''' Takes array of x and y values to create numpy meshgrid and plot the electric potentials. '''
    # Create a grid of x and y values
    x_grid, y_grid = np.meshgrid(xarr, yarr)
    # Set x and y limits (points in the middle so makes no sense to plot the whole grid)
    plt.xlim(40, 60)  # Replace x_min and x_max with your desired limits
    plt.ylim(40, 60)  # Replace y_min and y_max with your desired limits
    # Create a contour plot
    contour = plt.contour(x_grid, y_grid, electric_potentials, cmap='viridis', levels = 15)
    plt.clabel(contour, inline=True, fontsize=8) # to show values
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Contour Plot of Electric Potential')
    plt.show()

#### shows contour plot of electric potentials for q1 and q2
print("Plot of electric potentials and WIP of electric field")
electric_contour(xarr, yarr, e_potentials)

'''
PART 2 Calculate partial derivatives of potential w.r.t. x and y (the E-field in xy plane)
 d/dx = f(x+h/2, y) - f(x-h/2,y)
 d/dy = ...
 in layman terms need to get value of function (Potential) above and below each point where we want to calculate derivative [f(x+h) - f(x-h) / h] 
    where h is the step size (here it is a step size of 1 for the grid since we are going by 1 cm)
'''

# Assuming electric_potential is your 100x100 array of potential values

# partial derivs w.r.t x
ddx = np.zeros_like(e_potentials)
for i in range(e_potentials.shape[0]):
    ddx[i, :-1] = (e_potentials[i, 1:] - e_potentials[i, :-1])

# partial derivs w.r.t y
ddy = np.zeros_like(e_potentials)
for j in range(e_potentials.shape[1]):
    ddy[:-1, j] = (e_potentials[1:, j] - e_potentials[:-1, j])

print(ddx.shape)

# Generate grid of points
x = np.arange(100)
y = np.arange(100)
X, Y = np.meshgrid(x, y)

# Create quiver plot for both partial derivative values
plt.quiver(X, Y, ddx, ddy)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('E-field (partial derivatives of x and y)')
plt.axis('equal')

plt.show()





'''
# shape of electric potentials
original_shape = e_potentials.shape

# array to save partial derivative values
ddx = np.zeros((original_shape[1], original_shape[1] - 1))
#ddx = np.zeros((1,99))

#print('test ddx shape', ddx.shape)

# find partial derivatives using central difference theorem 
for i in range(original_shape[0]):
    for j in range(1, original_shape[1] - 1):
        ddx[i, j-1] = (e_potentials[i, j-1] - e_potentials[i, j+1]) - e_potentials[i, j-1]

# boundaries
for i in range(original_shape[0]):
    # Left boundary
    ddx[i, 0] = (e_potentials[i, 1] - e_potentials[i, 0])
    # Right boundary
    ddx[i, -1] = (e_potentials[i, -1] - e_potentials[i, -2])

# resulting array 
#print(ddx)

# array to store the partial derivative values
ddy = np.zeros((original_shape[0]-1, original_shape[1]))

# find the partial derivs using cdm
for i in range(1, original_shape[0] - 1):
    for j in range(original_shape[1]):
        ddy[i-1, j] = (e_potentials[i-1, j] - e_potentials[i+1, j]) - e_potentials[i-1, j]

# accounting for boundaries
for j in range(original_shape[1]):
    # Top boundary
    ddy[0, j] = (e_potentials[1, j] - e_potentials[0, j])
    # Bottom boundary
    ddy[-1, j] = (e_potentials[-1, j] - e_potentials[-2, j])

#####################

print("Shape of ddx:", ddx.shape)
print("Shape of ddy:", ddy.shape)


ddx2 = ddx[1:, :]
print("Shape of ddx2:", ddx2.shape)
ddy2 = ddy[:, :-1]
print("Shape of ddy2:", ddy2.shape)

# Generate grid of points
x = np.arange(100)
y = np.arange(100)
X, Y = np.meshgrid(x, y)

# Create quiver plot for both partial derivative values
plt.quiver(X, Y, ddx, ddy)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Partial Derivatives with respect to X and Y')
plt.axis('equal')

plt.show()



# Generate grid of points
x = np.arange(100)
y = np.arange(100)
X, Y = np.meshgrid(x, y)



# Create quiver plot for partial derivative with respect to x
plt.figure(figsize=(10, 5))  # Adjust figure size as needed
plt.subplot(1, 2, 1)  # Subplot for partial derivative with respect to x
plt.quiver(X[:-1, :], Y[:-1, :], ddx, np.zeros_like(ddx))
plt.xlabel('x')
plt.ylabel('y')
plt.title('E-field (x-component)')
plt.axis('equal')

# Create quiver plot for partial derivative with respect to y
plt.subplot(1, 2, 2)  # Subplot for partial derivative with respect to y
plt.quiver(X[:, :-1], Y[:, :-1], np.zeros_like(ddy), ddy)
plt.xlabel('x')
plt.ylabel('y')
plt.title('E-field (y-component)')
plt.axis('equal')

plt.tight_layout()  # Adjust layout to prevent overlap
plt.show()


# Create quiver plot for partial derivative with respect to x
plt.figure(figsize=(10, 5))  # Adjust figure size as needed
plt.subplot(1, 2, 1)  # Subplot for partial derivative with respect to x
plt.quiver(X[:-1, :], Y[:-1, :], ddx, ddy)
plt.xlabel('x')
plt.ylabel('y')
plt.title('E-field (x-component)')
plt.axis('equal')
#plt.show()
'''