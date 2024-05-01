'''
Exercise 8.7: Trajectory with air resistance
Many elementary mechanics problems deal with the physics of objects moving or flying through the air, but they almost always ignore friction and air resistance 
to make the equations solvable. If we're using a computer, however, we don't need solvable equations.
Consider, for instance, a spherical cannonball shot from a cannon standing on level ground. The air resistance on a moving sphere is a force in the opposite direction 
to the motion with magnitude

F = \frac{1}{2} \pi R^2\rho C v^2,

where R is the sphere's radius, rho is the density of air, v is the velocity, and C is the so-called coefficient of drag (a property of the shape of the moving object, in this case a sphere).

    Starting from Newton's second law, F = ma, show that the equations of motion for the position (x, y) of the cannonball are

    \ddot{x} = - {\pi R^2\rho C\over2m}\, \dot{x}\sqrt{\dot{x}^2+\dot{y}^2},
    \ddot{y} = - g - {\pi R^2\rho C\over2m}\, \dot{y}\sqrt{\dot{x}^2+\dot{y}^2},

    where m is the mass of the cannonball, g is the acceleration due to gravity, and \dot{x} and \ddot{x} are the first and second derivatives of x with respect to time.
    Change these two second-order equations into four first-order equations using the methods you have learned, then write a program that solves the equations for a 
    cannonball of mass 1 kg and radius 8 cm, shot at 30º to the horizontal with initial velocity 100 ms-1. The density of air is rho = 1.22 kg m-3 and the coefficient of 
    drag for a sphere is C = 0.47. Make a plot of the trajectory of the cannonball (i.e., a graph of y as a function of x).
    When one ignores air resistance, the distance traveled by a projectile does not depend on the mass of the projectile. In real life, however, mass certainly does make 
    a difference. Use your program to estimate the total distance traveled (over horizontal ground) by the cannonball above, and then experiment with the program to 
    determine whether the cannonball travels further if it is heavier or lighter. You could, for instance, plot a series of trajectories for cannonballs of different masses,
      or you could make a graph of distance traveled as a function of mass. Describe briefly what you discover.


'''
import numpy as np
import matplotlib.pyplot as plt 
import plotly.express as px
import pandas as pd

def RK4(f, y0, t0, tf, dt):
    """
    Generic Runge-Kutta 4th order (RK4) solver.
    =========================================================
    Parameters:
    f : callable
        The derivative function of the ODE (dy/dt = f(t, y)).
    y0 : float
        Initial value of the dependent variable.
    t0 : float
        Initial value of the independent variable.
    tf : float
        Final value of the independent variable.
    dt : float
        Step size (equivalent to h).

    Returns:
    t_values : list
        List of independent variable values.
    y_values : list
        List of solutions corresponding to t_values.
    """
    t_values = [t0]
    y_values = [y0]
    t = t0
    y = y0

    while t < tf:
        k1 = dt * f(t, y)
        k2 = dt * f(t + 0.5 * dt, y + 0.5 * k1)
        k3 = dt * f(t + 0.5 * dt, y + 0.5 * k2)
        k4 = dt * f(t + dt, y + k3)
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += dt
        t_values.append(t)
        y_values.append(y)

    return t_values, y_values


# set times and number of iterations, which we will use to find dt ("h")
t0 = 0 
tf = 10 # seconds
N = 1000
h = (tf - t0) / N

# Constants
g = 9.81  # gravity in m/s^2
R = 0.08  # Radius of the cannonball in meters
theta_0 = 30 * np.pi / 180  # Initial launch angle in radians (converted from degrees) -- 30 degrees to horizontal
v_0 = 100  # Initial velocity in m/s
rho = 1.22  # Density of air in kg/m^3
C = 0.47  # Drag coefficient for a sphere

t0 = 0  # Initial time in seconds
tf = 10  # Final time in seconds for this to run
N = 10000  # Number of time steps 
h = (tf - t0) / N  # Time step size ("dt" in RK4 function)

c = np.pi * R ** 2 * rho * C / 2  # Precomputed constant part of the drag force equation

def constant(m):
    """
    Calculate the constant used in the air resistance force component as a function of mass.
    -(pi R^2 rho C)/2 m

    Parameters:
        m (float): Mass of the projectile in kilograms.

    Returns:
        float: The overall constant divided by the given mass.
    """
    return c / m


def funct(r, t, m):
    """
    Calculates the derivatives for the ODE solver.

    Parameters:
        r (array): Current  vector [x, vx, y, vy] where x, y are positions and vx, vy are velocities w.r.t. x and y.
        t (float): Current time. 
        m (float): Mass of the projectile in kg.

    Returns:
        array: Array of derivatives [vx, ax, vy, ay] -- (or, x', x'', y', y'').
    """
    vx = r[1]
    vy = r[3]
    v = np.sqrt(vx ** 2 + vy ** 2)
    return np.array([vx, -constant(m) * vx * v,
                  vy, -g -constant(m) * vy * v], float)

def trajectory(m):
    """
    Calculates the trajectory of cannonball using the RK4 method.
    vx = v0 cos(theta)
    vy = v0 sin(theta)
    ============================================================================
    Parameters:
        m (float): Mass of the projectile in kilograms.

    Returns:
        tuple: Two arrays containing the x and y coordinates of the trajectory.
    """
    xpoints = []
    ypoints = []
    # Projectile begins at the origin(0,0) 
    r = np.array([0, v_0 * np.cos(theta_0), 0, v_0 * np.sin(theta_0)], float)
    for t in tpoints:
        xpoints.append(r[0])
        ypoints.append(r[2])
        k1 = h * funct(r, t, m)
        k2 = h * funct(r + 0.5 * k1, t + 0.5 * h, m)
        k3 = h * funct(r + 0.5 * k2, t + 0.5 * h, m)
        k4 = h * funct(r + k3, t + h, m)
        r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return np.array(xpoints, float), np.array(ypoints, float)

tpoints = np.arange(t0, tf, h)  # Time points for the RK4 method

# Calculate trajectories for different masses and plot them
trajectory1_x, trajectory1_y = trajectory(1)  # for 1 kg cannonball
trajectory2_x, trajectory2_y = trajectory(5)  # for 2 kg cannonball
trajectory3_x, trajectory3_y = trajectory(10)  # for 4 kg cannonball


# Calculate trajectories - requires using Pandas Dataframe for Plotly
masses = [1, 5, 10]
data = pd.DataFrame()
for m in masses:
    x, y = trajectory(m)
    df = pd.DataFrame({'x': x, 'y': y, 'Mass': f'{m} kg'})
    data = pd.concat([data, df])

# Plot using Plotly
fig = px.line(data, x='x', y='y', color='Mass', title='Trajectories of Cannonballs with Different Masses',
              labels={'x': 'Horizontal Distance (m)', 'y': 'Vertical Distance (m)'})
fig.show()


'''
NORMAL PLOTTING BELOW
plt.plot(trajectory1_x, trajectory1_y, 'k', label = '1 kg')  
plt.plot(trajectory2_x, trajectory2_y, 'b', label = '5 kg') 
plt.plot(trajectory3_x, trajectory3_y, 'g', label = '10 kg')  

plt.xlabel('x')  # Label the x-axis
plt.ylabel('y')  # Label the y-axis
plt.legend()
plt.show()  # Display the plot
'''
