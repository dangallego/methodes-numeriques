''' EXERCISE 8.1 -- RC Circuit
Write a program to solve this equation for V_out (t) using the fourht-order Runge-Kutta method when in the input signal 
is a square wave with frequence 1 and amplitude 1
'''

import numpy as np
import matplotlib.pyplot as plt
import os

def Vin(t):
    '''
    V_in function. 
    '''
    if (np.floor(2*t) % 2) == 0:    # Even case
        return 1
    elif (np.floor(2*t) % 2) != 0:  # Odd case
        return -1 
    
def rk2(Vout, t, h, RC):
    ''' Solves the ODE using the second-order Runga-Kutta method. '''
    # Computes k1
    k1 = h * (1/RC) * ( Vin(t) - Vout )
    # Computes k2
    k2 = h * (1/RC) * ( Vin(t + 0.5*h) - (Vout + 0.5*k1) )
    # Updates Vout using k2
    return Vout + k2 

def run_rk2(RC, t_end=10, h = 0.001):
    '''
    Function to run 2nd-order Runge-Kutta over a time interval for given h steps and RC value. 
    '''
    # Array of points over which ODE is solved, accounting for "h" steps from 0 to t_end
    times = np.arange(0, t_end, h)
    Vouts = np.zeros_like(times) # fills array of 0 same shape and size as "times" array 
    # Intial condition Vout(0) = 0
    Vouts[0] = 0

    for i in range(1, len(times)): # starting at 1 b/c first element is already 0
        '''For each iteration of loop, this calculates the value of Vout at the current time step, i , (which starts at 1, not 0)
            using the value of Vout from the previous time step, i-1 .  '''  
        Vouts[i] = rk2(Vouts[i-1], times[i-1], h, RC) 
    
    return times, Vouts

times, Vouts = run_rk2(RC=1, t_end=10, h = 0.01)

plt.scatter(times, Vouts,s=2, alpha = 0.8,c='darkslategrey', label = '2nd Order Runge-Kutta')
plt.xlabel('Time (s)') ; plt.ylabel('$V_{out}(t)$')
#plt.title('2nd Order Runge-Kutta Method')

#plt.show()

#### 4th Order Runge-Kutta Methods ####

def rk4(Vout, t, h, RC):
    ''' Solves the ODE using the fourth-order Runga-Kutta method. '''
    # Computes k1
    k1 = h * (1/RC) * ( Vin(t) - Vout )
    # Computes k2
    k2 = h * (1/RC) * ( Vin(t + 0.5*h) - (Vout + 0.5*k1) )
    # Computes k3
    k3 = h * (1/RC) * ( Vin(t + 0.5*h) - (Vout + 0.5*k2) )
    # Computes k4 (notice formula is different that above)
    k4 = h * (1/RC) * ( Vin(t + h) - (Vout + k3) )
    # Updates Vout using k1, k2, k3, k4
    return Vout + (1/6) * ( k1 + 2*k2 + 2*k3 + k4 )


def run_rk4(RC, t_end=10, h = 0.001):
    '''
    Function to run 4th-order Runge-Kutta over a time interval for given h steps and RC value. 
    '''
    # Array of points over which ODE is solved, accounting for "h" steps from 0 to t_end
    times = np.arange(0, t_end, h)
    Vouts = np.zeros_like(times) # fills array of 0 same shape and size as "times" array 
    # Intial condition Vout(0) = 0
    Vouts[0] = 0
    # same loop as for 2nd-order RK, because function above was all that had to change
    for i in range(1, len(times)): # starting at 1 b/c first element is already 0
        '''For each iteration of loop, this calculates the value of Vout at the current time step, i , (which starts at 1, not 0)
            using the value of Vout from the previous time step, i-1 .  '''  
        Vouts[i] = rk4(Vouts[i-1], times[i-1], h, RC) 
    
    return times, Vouts


times4, Vouts4 = run_rk4(RC=1, t_end=10, h = 0.01)

plt.scatter(times4, Vouts4, s=4,alpha=0.5, marker='v', c='salmon',label = '4th Order Runge-Kutta')
plt.xlabel('Time (s)') ; plt.ylabel('$V_{out}(t)$')
#plt.title('4th Order Runge-Kutta Method')
plt.legend()

# Save fig
save_dir = "/Users/Daniel/Downloads"
filename = "Circuit_RK2.png"
full_path = os.path.join(save_dir, filename)
#plt.savefig(full_path, dpi=300)  

plt.show()
