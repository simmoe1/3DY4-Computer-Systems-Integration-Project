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
from scipy.io import wavfile
from scipy import signal
import sys
from fourierTransform import my_DFT, my_IDFT
from filterDesign import lowPassImpulse

def convolution (x, h, a, zi):
    x = x / a
    h = h / a

    x_len = len(x)
    h_len = len(h)
    z_len = len(zi)

    length1 = h_len + x_len - 1
    #intalize y
    y = np.zeros(x_len)
    x_rev = x[::-1]
    h = np.pad(h, (x_len, x_len), 'constant')

    for n in range (len(y)):
        for k in range (len(h)):
            if n-k >= 0:
                y[n] += (h[k] * x_rev[n-k])

    y[:z_len] = y[:z_len] + zi


    return y[:h_len], y[h_len:]

def my_filter_BP (audio_data, \
					block_size, \
					audio_Fc, \
					audio_Fs, \
					N_taps):

	my_coeff =  lowPassImpulse(audio_Fc, audio_Fs, N_taps)

    #assume the data is stereo, as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
    #start at first block
	position = 0

	lengthCoeff = len(my_coeff)

	l_zi, r_zi = np.zeros(lengthCoeff - 1), np.zeros(lengthCoeff - 1)

	while True:

		filtered_data[position:position+block_size, 0], l_zi = convolution(my_coeff, audio_data[position:position+block_size, 0], 1.0, l_zi)

		filtered_data[position:position+block_size, 0], r_zi = convolution(my_coeff, audio_data[position:position+block_size, 0], 1.0, r_zi)

		position = position + block_size
		if position > len(audio_data):
			break

	return filtered_data



def filter_block_processing(audio_data, \
							block_size, \
							audio_Fc, \
							audio_Fs, \
							N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps,
								audio_Fc/(audio_Fs/2),
								window=('hann'))

	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)

	# start at the first block (with relative position zero)
	position = 0

	# intiial filter state - state is the size of the impulse response minus 1
	# we need to channels for the state (one for each audio channel)
	filter_state = np.zeros(shape = (len(firwin_coeff)-1, 2))

	r_zi = signal.lfilter_zi(firwin_coeff, 1.0)
	l_zi = signal.lfilter_zi(firwin_coeff, 1.0)

	while True:

		filtered_data[position:position+block_size, 0], l_zi = \
			signal.lfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 0], zi=l_zi)
		filtered_data[position:position+block_size, 1], r_zi = \
			signal.lfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 1], zi=r_zi)

		position += block_size

		# the last incomplete block is ignored
		if position > len(audio_data):
			break

	return filtered_data

def filter_single_pass(audio_data, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps,
								audio_Fc/(audio_Fs/2),
								window=('hann'))

	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# filter left channel
	filtered_data[:,0] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,1])

	return filtered_data

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\tth:  take-home')
	sys.exit()

# audio test file from: https://www.videvo.net/royalty-free-music/
if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	# use use wavfile from scipy.io for handling .wav files
	print('Opening audio file (.wav format)')
	audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	print(' Audio sample rate = {0:d} \
		\n Number of channels = {1:d} \
		\n Numbef of samples = {2:d}' \
		.format(int(audio_Fs), audio_data.ndim, len(audio_data)))

	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for processing streams divided in blocks')

		# you can control the cutoff frequency and number of taps
		single_pass_data = filter_single_pass(audio_data, \
							audio_Fc = 10e3, \
							audio_Fs = audio_Fs, \
							N_taps = 51)

		# write filtered data back to a .wav file
		wavfile.write("../data/single_pass_filtered.wav", \
					audio_Fs, \
					single_pass_data.astype(np.int16))

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for processing streams divided in blocks')

		# you can control also the block size
		block_processing_data = filter_block_processing(audio_data, \
							block_size = 51, \
							audio_Fc = 10e3, \
							audio_Fs = audio_Fs, \
							N_taps = 51)

		wavfile.write("../data/block_processing_filtered.wav", \
					audio_Fs, \
					block_processing_data.astype(np.int16))

	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for processing streams divided in blocks')

		block_processing_data = my_filter_BP(audio_data, \
							block_size = 51, \
							audio_Fc = 10e3, \
							audio_Fs = audio_Fs, \
							N_taps = 51)

		wavfile.write("../data/block_processing_filtered.wav", \
					audio_Fs, \
					block_processing_data.astype(np.int16))


		start = 0
		# number_of_samples = 5 * 1000

		fig, (ax1) = plt.subplots(1)
		ax1.plot(block_processing_data[start:start+200])
		ax1.set(title='block')

		plt.show()



	else:

		cli_error_msg()

	# plt.show()
