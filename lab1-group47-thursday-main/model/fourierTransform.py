#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
import cmath, math
import sys
from scipy import signal

def plotSpectrum(x, Fs, type = 'FFT'):

    n = len(x)             # length of the signal
    df = Fs/n              # frequency increment (width of freq bin)

    # compute Fourier transform, its magnitude and normalize it before plotting
    if type == 'FFT':
        Xfreq = np.fft.fft(x)
    elif type == 'IFFT':
        Xfreq = np.fft.ifft(x)
    elif type == 'DFT':
        Xfreq = my_DFT(x)
    elif type == 'IDFT':
        Xfreq = my_IDFT(x)
    XMag = abs(Xfreq)/n

    # Note: because x is real, we keep only the positive half of the spectrum
    # Note also: half of the energy is in the negative half (not plotted)
    XMag = XMag[0:int(n/2)]

    # freq vector up to Nyquist freq (half of the sample rate)
    freq = np.arange(0, Fs/2, df)

    fig, ax = plt.subplots()
    ax.plot(freq, XMag)
    ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
        title='Frequency domain plot')
    # fig.savefig("freq.png")
    plt.show()

def plotTime(x, time):

    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

    dt = 1.0/Fs                          # sampling period (increment in time)
    time = np.arange(0, interval, dt)    # time vector over interval

    # generate the sin signal
    x = amplitude*np.sin(2*math.pi*frequency*time+phase)

    return time, x

def cli_error_msg():

    # error message to provide the correct command line interface (CLI) arguments
    print('Valid arguments:')
    print('\trc:  reference code')
    print('\til1: in-lab 1')
    print('\til2: in-lab 2')
    print('\til3: in-lab 3')
    print('\tth:  take-home')
    sys.exit()

def my_IFFT (x):
    Xf = np.fft.ifft(x)
    Xf = X
    return Xf

def my_DFT (x):
    n=len(x)
    X = np.zeros(n, dtype=np.complex_)
    for m in range(n):
        X[m] = 0
        for k in range(n):
            X[m] = X[m] + x[k]*cmath.exp(-2j*np.pi*((k*m)/n))
    Xf = X
    return Xf

def my_IDFT (x):
    n = len(x)
    Xf = np.zeros(n, dtype=np.complex)
    for m in range(n):
        for k in range(n):
            Xf[m] = Xf[m] + x[k]*cmath.exp(2j*np.pi*k*m/n)
        Xf[m] = Xf[m]/n
    return Xf

def signalEn (x):
    energy = 0
    for i in range (len(x)):
        energy = energy + abs(x[i])**2
    return energy

def multiTone(Fs, interval, frequency, amplitude, phase):
    time, tone1 = generateSin (Fs, interval, 5, 11, 0.0)
    time, tone2 = generateSin (Fs, interval, 8, 4, 2.0)
    time, tone3 = generateSin (Fs, interval, 10, 16, 3.0)
    x = tone1 + tone2 + tone3

    return time, x

#def squareWave (Fs, interval, frequency , duty):
#    dt = 1.0/Fs
#    time = np.arange(0, interval, dt)
#   rest of function
#    return time, x

if __name__ == "__main__":

    if len(sys.argv[0:]) != 2:
        cli_error_msg()

    Fs = 100.0          # sampling rate
    interval = 1.0      # set up to one full second

    if (sys.argv[1] == 'rc'): # runs the reference code (rc)

        print('Reference code for the Fourier transform')

        # generate the user-defined sin function
        time, x = generateSin(Fs, interval)
        # plot the signal in time domain
        plotTime(x, time)
        # plot the signal in frequency domain
        plotSpectrum(x, Fs, type = 'FFT')

    elif (sys.argv[1] == 'il1'):

        print('In-lab experiment 1 for the Fourier transform')
        # compute the spectrum with your own DFT
        # you can use cmath.exp() for complex exponentials
        # plotSpectrum(x, Fs, type = 'your DFT name')
        # confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
        # for further details, if any, check the lab document
        time, x = generateSin(Fs, interval)
        plotTime(x, time)
        #comment out 1 by 1 to check IFT and DFT and to check IDF
        plotSpectrum(x, Fs, type = 'DFT')
        plotSpectrum(x, Fs, type = 'FFT')
        plotSpectrum(x, Fs, type = 'IDFT')
        plotSpectrum(x, Fs, type = 'IFFT')
        DFT1 = my_DFT(x)
        IDFT1 = my_IDFT(x)

        if np.allclose(np.real(DFT1), np.real(IDFT1), rtol = 1.e-3, atol=1.e-3, equal_nan = False): #All Close
            print('DFT and IDFT is correctness.')

    elif (sys.argv[1] == 'il2'):

        print('In-lab experiment 2 for the Fourier transform')

        # use np.random.randn() for randomization
        # we can owverwrie the default values

        Fs =1000.0
        interval = 1.0
        frequency = 8.0                      # frequency of the signal
        amplitude = (20)*np.random.rand()-10 # amplitude of the signal
        phase = 1.0                           # phase of the signal
        time, x = generateSin(Fs, interval, frequency, amplitude, phase)

        plotTime(x, time)
        plotSpectrum(x, Fs, type = 'DFT')
        x = 10*np.random.normal(size=1000)
        print(signalEn(x), signalEn(my_DFT(x))/len(x), signalEn(my_IDFT(my_DFT(x))))
        # You should also numerically check if the signal energy
        # in time and frequency domains is identical

        # for further details, if any, check the lab document

    elif (sys.argv[1] == 'il3'):

        print('In-lab experiment 3 for the Fourier transform')

        # generate randomized multi-tone signals
        # plot them in both time and frequency domain
        Fs =1000.0
        interval = 1.0
        frequency = 8.0                      # frequency of the signal
        amplitude = (20)*np.random.rand()-10 # amplitude of the signal
        phase = 1.0
        time, x = multiTone(Fs, interval, frequency, amplitude, phase)
        plotTime(x, time)
        plotSpectrum(x, Fs, type = 'FFT')

        # for further details, if any, check the lab document

    elif (sys.argv[1] == 'th'):

        print('Take-home exercise for the Fourier transform')

        #time, x = squareWave(Fs, interval, frequency = 2, duty = 0.5)
        t = np.linspace(0, 1, 500, endpoint = False)
        plt.plot(t, signal.square(2 * np.pi * 1 *t, duty = 0.25))
        plt.ylim(-2, 2)
        plotSpectrum(t, 1000, type = 'FFT')
        #plotTime(x, time)
        #plotSpectrum(x, Fs, type = 'FFT')

    else:

        cli_error_msg()

    plt.show()
