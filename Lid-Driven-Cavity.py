'''
LID DRIVEN CAVITY NAVIER STOKES SOLVER / PLOTTER

Momentum:           ∂u/∂t + (u ⋅ ∇) u = - 1/rho ∇p + nu ∇²u + f

Incompressibility:  ∇ ⋅ u = 0

Begin with fluid at rest and then horizontal velocity applied on top of unit square

STEPS / EQ'ns (Chorin's splitting projection)

u = [u, v]
x = [x, y]

1. Solve for current velocity with velocity boundary condition

    ∂u/∂t + u ∂u/∂x + v ∂u/∂y = nu ∂²u/∂x² + nu ∂²u/∂y²

    ∂v/∂t + u ∂v/∂x + v ∂v/∂y = nu ∂²v/∂x² + nu ∂²v/∂y²

2. Solve pressure poisson with pressue boundary condition (same everywhere except for the top)

    ∂²p/∂x² + ∂²p/∂y² = rho/Δt (∂u/∂x + ∂v/∂y)

3. Correct velocity (again with velocity boundary conditions)

    u ← u - Δt/rho ∂p/∂x

    v ← v - Δt/rho ∂p/∂y

'''

import matplotlib.pyplot as plt
import numpy as np
import argparse
from tqdm import tqdm


# Set up the argument parser
parser = argparse.ArgumentParser(description="Plots a lid-driven cavity flow problem with user defined parameters.")

# Add arguments
parser.add_argument("--n_points", type=int, default=50, help="Number of points in each dimension of the unit square grid.")
parser.add_argument("--time_step", type=float, default=0.001, help="Time step length (h).")
parser.add_argument("--n_iterations", type=int, default=1000, help="Number of iterations for the simulation.")
parser.add_argument("--kvisc", type=float, default=0.1, help="Kinematic Viscosity.")
parser.add_argument("--top_vel", type=float, default=1.0, help="Horizontal velocity applied at top of lid.")
parser.add_argument("--n_poisson", type=int, default=100, help="Number of iterations to use when solving pressure Poisson equation.")
parser.add_argument("--rho", type=float, default=1.0, help="Density of fluid.")


# Parse arguments
args = parser.parse_args()

# Set constants and parameters
N_POINTS = args.n_points
h = args.time_step
N_ITERATIONS = args.n_iterations
KVISC = args.kvisc
VELOCITY_TOP = args.top_vel
RHO = args.rho
POISSON_ITERATIONS = args.n_poisson


# Calculate element length and create grid
element_length = 1 / (N_POINTS - 1) # space between points in grid ; denom. accounts for "picket fence" problem 
x = np.linspace(0, 1, N_POINTS) # creates array from 0 to 1 (unit square), with N_POINTS evenly spaced samples
y = np.linspace(0, 1, N_POINTS)
X, Y = np.meshgrid(x, y)

# Initialized velocity and pressure fields: fluid at rest 
u = np.zeros_like(X) # velocity in x-direction
v = np.zeros_like(Y) # velocity in y-direction
p = np.zeros_like(X) # pressure (scalar)


# Define functions for finite differences
def cdf_x(funct):
    '''
    For a given array f, calculates the first derivative with respect to the x-direction 
    using the central difference formula. 
    '''
    dfdx = np.zeros_like(funct)
    # Operate only on interior points to avoid boundary issues
    dfdx[:, 1:-1] = (funct[:, 2:] - funct[:, :-2]) / (2 * element_length)# for 2D array: [row, column] == [y,x] , which can be confusing but means we hold y constant 
    return dfdx

def cdf_y(funct):
    '''
    For a given array f, calculates the first derivative with respect to the y-direction 
    using the central difference formula. 
    '''
    dfdy = np.zeros_like(funct)
    # Operate only on interior points to avoid boundary issues
    dfdy[1:-1, :] = (funct[2:, :] - funct[:-2, :]) / (2 * element_length)
    return dfdy

def laplacian(funct):
    '''
    Laplacian operator for 2D grid, found using

    ∇²f = f(x+a, y) + f(x-a, y) + f(x, y+a) + f(x, y-a) - 4*f(x,y) / a² , where a = element_length
    '''
    laplacian = np.zeros_like(funct)
    # Calculate the Laplacian for internal points
    laplacian[1:-1, 1:-1] = (funct[1:-1, 2:] + funct[1:-1, :-2] + funct[2:, 1:-1] + funct[:-2, 1:-1] - 4 * funct[1:-1, 1:-1]) / (element_length**2)
    return laplacian


# Chorin's Projection Calculation / Steps
for _ in tqdm(range(N_ITERATIONS)): # to get progress bar
    # Calculate spatial derivatives -- first solve momentum without pressure gradient to get tent. velocity
    du_d_x = cdf_x(u)
    du_d_y = cdf_y(u)
    dv_d_x = cdf_x(v)
    dv_d_y = cdf_y(v)
    laplace_u = laplacian(u)
    laplace_v = laplacian(v)

    # Tentative velocity step (u*, v*) -- first solve momentum WITHOUT pressure gradient
    u_tent = u - h * (u * du_d_x + v * du_d_y) + KVISC * laplace_u * h 
    v_tent = v - h * (u * dv_d_x + v * dv_d_y) + KVISC * laplace_v * h

    # Apply boundary conditions everywhere except at top of cavity 
    u_tent[0, :] = 0 # bottom of cavity
    u_tent[:, 0] = 0 # left boundary
    u_tent[:, -1] = 0 # right boundary
    u_tent[-1, :] = VELOCITY_TOP # horizontal velocity at top of cavity

    v_tent[0, :] = v_tent[:, 0] = v_tent[:, -1] = v_tent[-1, :] = 0 # quickhand way of setting all to 0 

    # Compute the divergence of tentative velocity field
    divergence = cdf_x(u_tent) + cdf_y(v_tent)

    for _ in range(POISSON_ITERATIONS):
        p_next = np.zeros_like(p)
        # interior points 
        p_next[1:-1, 1:-1] = ( p[1:-1, 0:-2]+ p[0:-2, 1:-1]+ p[1:-1, 2:]+ p[2:, 1:-1] - element_length**2 
                                    * (RHO/h * divergence)[1:-1, 1:-1] )/4
        
        # Apply pressure boundary conditions 
        p_next[:, -1] = 0 # boundary value becomes that of the next interior point
        p_next[0, :] = 0
        p_next[-1, :] = 0
        p_next[:, 0] = 0
        p = p_next # update pressure 


    # Correct velocities using pressure gradient -- need derivatives of pressure field w.r.t axes 
    d_p_next_dx = cdf_x(p_next)
    d_p_next_dy = cdf_y(p_next)

    # Correct velocities using the pressure gradient
    u_next = u_tent - h / RHO * d_p_next_dx
    v_next = v_tent - h / RHO * d_p_next_dy

    # Apply velocity boundary conditions (again)
    u_next[0, :] = u_next[:, 0] = u_next[:, -1] = 0
    u_next[-1, :] = VELOCITY_TOP
    v_next[0, :] = v_next[:, 0] = v_next[:, -1] = v_next[-1, :] = 0

    # Update final 
    u = u_next
    v = v_next
    p = p_next

# Plotting results
plt.figure(figsize=(8, 8))
plt.contourf(X, Y, p_next, levels = 20, cmap = 'coolwarm')
plt.colorbar(label='Pressure')
plt.quiver(X[::2, ::2], Y[::2, ::2], u_next[::2, ::2], v_next[::2, ::2], color="k")
plt.title("Velocity and Pressure Field in Lid-Driven Cavity")
plt.xlabel("X")
plt.ylabel("Y")

# save figure 
#plt.savefig('cavity_2.png', format = 'png', dpi=300)
plt.show()




