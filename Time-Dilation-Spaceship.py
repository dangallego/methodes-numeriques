
### EXERCISE 2.4 -- SPECIAL RELATIVITY SPACESHIP
'''A spaceship travels from Earth in a straight line at relativistic speed v to another planet x light years away. 
    Write a program to ask the user for the value of x and the speed v as a fraction of the speed of light c, then 
    print out the time in years that the spaceship takes to reach its destination (a) in the rest frame of an observer on Earth
    and (b) as perceived by a passenger on board the ship. Use your program to calculate the answers for a planet 10 light years 
    away with v = 0.99c'''

import numpy as np
import argparse

def time_dilation_calculator(x, v):
    '''Calculates the time (in years) that it takes for a spaceship to get from Earth to a planet 
        X light years away. Returns the time it takes from the perpsective of an observer on Earth (del_t)
        and the time from the perpsective of someone on board the ship (in a moving frame of reference, t0). 
        The distance X and speed of the spaceship (in lightyears) should be given by the user. Speed should
        be given in terms of a fraction of the speed of light (eg. 0.1, 0.99).'''
    c = 299_792_458 # meters/s; the speed of light 
    vfrac = v*c # velocity is now a fraction the speed of light 
    # time calculation of moving clock; ie. time as observed by outsided stationary observer on Earth:
    del_t = x/vfrac # standard eq. for: time = change in position / velocity 
    gamma = 1 / np.sqrt(1 - (vfrac/c)**2) # calculates gamma factory for eq: del_t = gamma * t0

    t0 = del_t / gamma # stationary clock time calculation for observer in the spaceship
    text1 = f"t0 (time observed from inside spaceship) is {t0:.2e} years, "  # formats text to only show 2 decimal points in scientific notation
    text2 = f"and del_t (time observed from Earth) is {del_t:.2e}"
    return text1 + "\n" + text2 # expected: t0 <= del_t (del_t will be larger b/c more time elapsed relative to t0)

#argparse section
parser = argparse.ArgumentParser(
                    prog='Time-Dilation-Spaceship',
                    description='Prints the relative and "normal" time for a spaceship travlling from Earth to another planet, \
                    as observed from Earth (del_t) and from inside the spaceship (t0).',
                    epilog='Text at the bottom of help')

parser.add_argument('-x', '--distance', type = float,  # option that takes a distance, must be a float
                    action='store', help='Must enter the distance as a float/decimal. If no distance indicated \
                    in command line as argument, program will ask for user input.')  #added help value
parser.add_argument('-v', '--velocity', type = float,  # option that takes a velocity
                    action='store', help='Must enter speed of spaceship as a fraction of the speed of light \
                        If no velocity is indicated in command line, program will ask for user input.') # default value set to very fast (99% speed of light)  

args = parser.parse_args()
#print(args.velocity, args.distance)

# conditional statement to ask user for their input if optional arguments not passed in command line
if args.distance is None or args.velocity is None: 
    x = float(input("Enter the distance from Earth to the planet of travel \n"))
    v = float(input("Enter the speed of the spaceship, as a fraction of the speed of light \n"))
    print(time_dilation_calculator(x,v))
else: 
    print(time_dilation_calculator(args.distance, args.velocity))



