'''
Exercise 10.3: Brownian motion
Brownian motion is the motion of a particle, such as a smoke or dust particle, in a gas, as it is buffeted by random collisions with gas molecules. 
Make a simple computer simulation of such a particle in two dimensions as follows. The particle is confined to a square grid or lattice L x L squares on a side, 
so that its position can be represented by two integers i, j = 0 . . . L - 1. It starts in the middle of the grid. On each step of the simulation, 
choose a random direction—up, down, left, or right—and move the particle one step in that direction. This process is called a random walk. The particle is 
not allowed to move outside the limits of the lattice—if it tries to do so, choose a new random direction to move in.
Write a program to perform a million steps of this process on a lattice with L = 101 and make an animation on the screen of the position of the particle. 
(We choose an odd length for the side of the square so that there is one lattice site exactly in the center). 

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# Grid size
L = 101  
# Number of steps to simulate: 
steps = 1_000_000
# Start at middle of the grid
start_position = [L // 2, L // 2]  # should be 50 for both 

# Initialize the position of the particle
position = np.array(start_position)

# Possible directions to move : up, down, left, right
moves = np.array([[0, 1], [0, -1], [-1, 0], [1, 0]])

# Stores path if we want to plot 
path = []

# Perform the random walk
for _ in range(steps):
    while True:
        # Choose random direction
        move = moves[np.random.randint(0, 4)]  
        # Calculate new position
        new_position = position + move  
        # Check that new position is within the grid
        if 0 <= new_position[0] < L and 0 <= new_position[1] < L:
            #Updates position
            position = new_position  
            #exits the loop to perform next step
            break  
    path.append(position.copy())  # Save position for plotting

path = np.array(path)

#print(path.shape)

x = path[:,0]
y = path[:,1]



''''
BELOW WAS IF WE WANTED TO SHOW A STATIC END PLOT THAT TRACED PATH! 
# Static path plot
plt.figure(figsize=(8, 8))
plt.plot(x, y color='slategrey', linewidth=0.05, alpha=0.4)
plt.scatter(start_position[0], start_position[1], color='green', label='Start')
plt.scatter(path[-1, 0], path[-1, 1], color='red', label='End')
plt.title("Brownian Motion")
plt.xlabel("X")
plt.ylabel("Y")
plt.legend()
plt.show()
'''

# Create animation 
fig = plt.figure(figsize=(8,8,))
# Limits should be size of grid
axs= plt.axes(xlim=(0,100),ylim=(0,100))
point=plt.Circle((0,0),radius=1,facecolor='darkslategrey')
plt.title('Brownian Motion Random Walk')
plt.xlabel('x') ; plt.ylabel('y')
axs.add_patch(point)

def init():
    point.center = (0, 0) 
    axs.add_patch(point)
    return point

def animate(i):
    x_animation = x[i] 
    y_animation = y[i]
    point.center = (x_animation, y_animation)
    return point,


anim = animation.FuncAnimation(fig, animate,frames=len(x), init_func= init)
plt.show()
