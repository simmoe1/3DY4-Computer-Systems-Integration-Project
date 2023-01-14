#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import math
import sys

# use generateSin/plotTime from the fourierTransform module
from fourierTransform import generateSin, plotTime, multiTone, plotSpectrum

def lowPassImpulse (Fc, Fs, Ntaps):

	normCut = Fc / (Fs/2) #equation
	h = np.zeros(Ntaps)

	for i in range (Ntaps):
		if i == ((Ntaps-1)/2):
			h[i] = normCut
		else:
			variable = math.pi * normCut * (i - ((Ntaps-1)/2))
			h[i] = normCut * (np.sin(variable) / variable)
		h[i] = h[i] * (np.sin((i*math.pi)/Ntaps)) * (np.sin((i*math.pi)/Ntaps))

	return h

def freqzPlot(coeff, Fs, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
	w, h = signal.freqz(coeff)

	# Reminder: np.pi rad/sample is actually the Nyquist frequency
	w = w * Fs/(2*np.pi) # needed to draw the frequency on the X axis

	# plots the magnitude response where the x axis is normalized in rad/sample
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color='b')
	ax1.set_xlabel('Frequency [Hz]')

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response is for some applications,
	# it is not critical at this stage because we expect a linear phase in the passband

	# ax2 = ax1.twinx()
	# angles = np.unwrap(np.angle(h))
	# ax2.plot(w, angles, 'g')
	# ax2.set_ylabel('Angle (radians)', color='g')

def filterMTSignal (Fs, Fc, coeff):

	phase = [1.0, 2.0, 3.0]
	frequency = [5, 10, 15]
	amplitude = [1, 2, 3]
	time, x = multiTone(Fs, 1.0, frequency, amplitude, phase)
	plotTime(x, time)
	plotSpectrum(x, Fs, type = 'FFT')

	fx = signal.lfilter(coeff, 1.0, x)
	plotTime(fx, time)
	plotSpectrum(fx, Fs, type = 'FFT')


def filterSin(Fs, Fc, coeff):

	# we can control the frequency relative to the filter cutoff
	time, x = generateSin(Fs, interval = 1.0, frequency = Fc * 0.4)
	plotTime(x, time)

	# use lfilter from SciPy for FIR filtering:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
	fx = signal.lfilter(coeff, 1.0, x)

	# you should clearly observe the effects (attenuation, delay) introduced by the filter
	plotTime(fx, time)


def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\tth:  take-home')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0           # sampling rate
	Fc = 15.0            # cutoff frequency
	N_taps = 41          # number of taps for the FIR


	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for the digital filter design')

		# derive filter coefficients using firwin from Scipy:
		# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
		# second argument is the normalized cutoff frequency, i.e., the
		# cutoff frequency divided by Nyquist frequency (half of sampling rate)
		firwin_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))

		# plot the frequency response obtained through freqz
		freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the digital filter design')

		# implement your own method for finding the coefficients for a low pass filter
		# my_own_coeff = ... provide the following arguments: Fc, Fs and N_taps
		# compare through visual inspection the frequency response against firwin
		# freqzPlot(my_own_coeff, Fs, 'my own FIR design with ' + str(N_taps) + ' taps')
		my_own_coeff = lowPassImpulse(Fc, Fs, N_taps)
		freqzPlot(my_own_coeff, Fs, 'my own FIR design with ' + str(N_taps) + ' taps')


	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the digital filter design')

		# you can confirm that a single tone has been filtered
		my_own_coeff = lowPassImpulse(Fc, Fs, N_taps)
		#filterSin(Fs, Fc, my_own_coeff)
		filterMTSignal(Fs, Fc, my_own_coeff)
		#filterMTSignal(Fs, Fc, firwin_coeff)


	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for the digital filter design')

		# for specific details check the lab document
		#implement band pass filter
		#plot frequency response and look at filtered data in time domain
		#pass only middle tone of 3 tone signal
		firwin_coeff = signal.firwin(N_taps, [5.0, 15.0], pass_zero=False, window = 'hann', fs=Fs)
		freqzPlot(firwin_coeff, Fs, 'firwin for ' + str('5hz and 15hz') + ' Hz cutoff with ' + str(N_taps) + ' taps')
		filterMTSignal(Fs, Fc, firwin_coeff)


	else:

		cli_error_msg()

	plt.show()
