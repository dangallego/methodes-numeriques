''' Exercise 7.1: Fourier transforms of simple functions

Write Python programs to calculate the coefficients in the discrete Fourier transforms of the following periodic functions sampled at 
N = 1000 evenly spaced points, and make plots of their amplitudes:

a) A single cycle of a square-wave with amplitude 1
b) The sawtooth wave yn = n
c) The modulated sine wave yn = sin(πn/N) sin(20πn/N)

'''
import numpy as np
import matplotlib.pyplot as plt

def dft(y): 
    '''
    Computes and returns coefficient values from Discrete Fourier Transform (DFT). 
    -------------------------------------------------------------------------------
    Parameters: 
    y (array_like): input values on which we want to perform the DFT. 

    Returns: 
    numpy.ndarray: Array of coefficients from the DFT. 
    '''
    N = len(y)
    c = np.zeros(N//2+1, complex)

    for k in range(N//2+1):
        for n in range(N):
            c[k] += y[n]*np.exp(-2j*np.pi*k*n/N)
    return c

def dft_plots(x, y, percent_coeff=100, function_title='Function'):
    '''
    Generates original plot of a function and uses defined Discrete Fourier Transform (DFT)
    function to calculate and plot the Power Spectrum of DFT coefficients. Also plots the 
    Inverse Fast Fourier Transform using Numpy's function to check that original function can be retrieved. 
    ---------------------------------------------------------------------------------------------------------
    Parameters: 
    x (array_like): range of values used to evaluate function we are calculating DFT over. 
    y (array_like): input values used to calculate DFT. 
    percent_coeff (int): percentage of beginning coefficients desired to be displayed for the Power Spectrum. 
    function_title (Str): title of function for first plot, which is the function we are performing DFT and IFFT on. 

    Returns: 
    3 plots
    '''
    # Calculate DFT of y values
    coeff = dft(y)

    # Calculate |magnitude| of coefficients from DFT
    mod_coeff = np.abs(coeff)
    # FOR SOME REASON THE BELOW DOES NOT WORK (though it should)
    #coeff = np.conjugate(dft(y)) * dft(y)

    # If we want to keep only low frequency (first 20% of coefficients)
    mod_coeff = mod_coeff[:len(mod_coeff) * percent_coeff // 100]
    #mod_coeff[(len(mod_coeff) // 20):] = 0 # this code changes those values to zero 

    # Now use Numpy Inverse FFT function to check that we can get back original function
    inv_fft = np.fft.irfft(coeff)

    # Create a figure and a grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))

    # Plot original function on top left
    axs[0, 0].plot(x,y)
    axs[0, 0].set_title(function_title)

    # Plot Power Spectrum on top right
    axs[0, 1].plot(mod_coeff)
    axs[0, 1].set_title('Power Spectrum')

    # Plot Inverse FFT on the third subplot 
    axs[1, 0].plot(inv_fft, c='orangered')
    axs[1, 0].set_title('Inverse FFT')

    # Hide the unused subplot (bottom right)
    axs[1, 1].axis('off')

    plt.tight_layout()
    plt.show()

#### a) Single cycle of a square-wave with amplitude 1
def square_wave(x):
    ''' 
    Generates a square wave with amplitude 1. 
    -------------------------------------------------------------
    x (array_like): array of time points at which wave is evaluated
    '''
    return np.sign(np.sin(2 * np.pi * x))

# # Generate x values and calculate y values for Square Wave funct
N = 1000
x = np.linspace(0, 1, N)  # Generate N points in time from 0 to 1 x values
square_y = square_wave(x)  # calculate function values outside of plotting function

# Print all 3 plots
dft_plots(x, square_y, percent_coeff=20, function_title='Square Wave Function')


### Part b) Sawtooth wave
def sawtooth_wave(length, num_periods=5):
    '''
    Generates a sawtooth wave. 
    --------------------------------------------
    Parameters: 
    length (int): Total desired number of points in wave function. 
    num_periods (int): Number of desired periods for sawtooth wave to generate. 

    Returns: 
    numpy.ndarray: An array of output values of the sawtooth function. 
    '''
    n = np.arange(length) # creates array of desired length
    # Period is how many sample points there are per period of sawtooth (integer division because this number must be an integer)
    period = length // num_periods
    # Following creates repeating pattern
    sawtooth = n % period
    # Normalize to values between 0 and 1 
    sawtooth = sawtooth / period

    return sawtooth

# Generate x values and calculate y values for Sawtooth Wave funct
N = 1000
x = np.arange(0, N) 
sawtooth_y = sawtooth_wave(N)

# Print all 3 plots
dft_plots(x, sawtooth_y, percent_coeff=50, function_title='Sawtooth Wave Function')

### Part c) The modulated sine wave yn = sin(πn/N) sin(20πn/N)
import numpy as np

def sine_wave(n, N=10):
    '''
    Generates a modulated sine wave: y_n = sin(πn/N) sin(20πn/N).
    -------------------------------------------------------------------------------
    Parameters:
    n (array_like): The values of points in the wave. 
    N (int): The total number of points in the wave. 

    Returns:
    numpy.ndarray: An array containing the values of the generate sine wave.
    '''
    y_n = np.sin(np.pi * n / N) * np.sin(20 * np.pi * n / N)
    return y_n


N = 1000  # Total number of points
x = np.arange(N)  
sine_y = sine_wave(n=x, N=10)

# Print all 3 plots
dft_plots(x, sine_y, percent_coeff=100, function_title='Sine Wave Function')

print('The following are plots of three different functions and their Discrete Fourier Transforms (DFT), as well as',
      'their corresponding Power Spectra and Inverse Fast Fourier Transform, to show that original function can be recovered from DFT. ')


