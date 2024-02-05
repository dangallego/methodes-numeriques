## Exercise 2.1: Another ball dropped from a tower
''' A ball is dropped from a tower of height h with initial velocity zero. 
Write a program that asks the user to enter the height in meters of the tower and then 
calculates and prints the time the ball takes until it hits the ground, ignoring air resistance. 
Use your program to calculate the time for a ball dropped from a 100 m high tower. '''

import numpy as np
import sys

g = 10 
# if using argv - no need to input 
#h = float(input("Enter the height of the tower that the ball is being dropped from (in meters): \n"))
#h = float(sys.argv[1])

if len(sys.argv) != 1: #checks if length of arguments returned is not 1 (in which case we entered the height in the terminal after the file)
    h = float(sys.argv[1])
else: #if we just run the file, this condition prompts us to enter the distance the ball is dropped
    h = float(input("Enter the height of the tower that the ball is being dropped from (in meters): \n"))

#if type(h) != float or type(h) != int: ### THIS WON'T WORK bc type is always originally string (then converted) --> argparse works better here
#    print("You must type a number") ### issues here actually arise from not being able to convert (ie. user typing '75 meters' instead of '75')

t = np.sqrt((2*h)/(g)) #equation to find time 

#print("The total time the ball takes to hit the ground is: ", t, " seconds")
text = "The total time the ball takes to hit the ground is {:.2f} seconds!".format(t) # formats text to print time to only show last two digits 
print(text, "\n", sys.argv, "length:", len(sys.argv))



