#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

from fmSupportLib import fmDemodArctan, fmPlotPSD
from fmSupportLib import secondaryFmDemodArctan
from fmSupportLib import pcoeffBPIR
from fmPll import fmPll

#variables
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

# add other settings for audio, like filter taps, ...

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    in_fname = "../src/stereo_l0_r9.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')

	#print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (raw_data - 128.0) / 128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc / (rf_Fs / 2), window='hann')

    # coefficients for the filter to extract mono audio
    audio_coeff = signal.firwin(audio_taps, audio_Fc / ((rf_Fs / rf_decim) / 2), window='hann')

	#stereo channel recovery and extraction filtering
    audio_coeffRec  = pcoeffBPIR(18.5e3, 19.5e3, 240e3, audio_taps)
    audio_coeffEtr  = pcoeffBPIR(22e3, 54e3, 240e3, audio_taps)

	# set up the subfigures for plotting
    subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
    plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
    block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps - 1)
    state_q_lpf_100k = np.zeros(rf_taps - 1)
    state_lpf_16k = np.zeros(audio_taps - 1)
    state_lpf_16k_s = np.zeros(audio_taps - 1)

    # states for stereo - recover and extraction
    state_monoRec   = np.zeros(audio_taps - 1)
    state_monoEtr  = np.zeros(audio_taps - 1)
    stereo_filt = np.zeros(audio_taps - 1)

	# used for left/right audio blocks
    audio_left = np.array([])
    audio_right = np.array([])

    # beginning states for pll
    integrator = 0.0
    phaseEst = 0.0
    feedbackI = 1.0
    feedbackQ = 0.0
    trigOffset = 0
    ncoOut = np.zeros(1)

    # audio buffer that stores all the audio blocks
    stereo_mixed  = np.zeros(int(block_size / 2 / rf_decim))
    mono_audio_data = np.array([])
    Stereo_data = np.array([]) # used to concatenate filtered blocks (audio data)

    # if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
    while (block_count + 1) * block_size < len(iq_data):
        # if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
        print('Processing block ' + str(block_count))

        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size:(block_count + 1) * block_size:2],
                                                  zi=state_i_lpf_100k)
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, iq_data[block_count * block_size + 1:(block_count + 1) * block_size:2],
                                                  zi=state_q_lpf_100k)

        # downsample the I/Q data from the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]

        # FM demodulator custom demod
        fm_demod, prev_i, prev_q = secondaryFmDemodArctan(i_ds, q_ds)

        # stereo - extraction first then recovery
        stereo_extract, state_monoEtr  = signal.lfilter(audio_coeffEtr , 1.0, fm_demod, zi=state_monoEtr )
        stereo_recovery, state_monoRec   = signal.lfilter(audio_coeffRec , 1.0, fm_demod, zi=state_monoRec  )

        ncoOut, integrator, phaseEst, feedbackI, feedbackQe, trigOffset_ = fmPll(stereo_recovery, 19e3, 240e3, integrator, phaseEst, feedbackI, feedbackQ, trigOffset, ncoOut, 2.0)

        # stereo-processing mixer
        stereo_mixed  = 2 * stereo_extract * ncoOut[:-1]

		# digital filtering and downsampling
		# audio_block = ... change as needed
        stereo_audio_filt, stereo_filt = signal.lfilter(audio_coeff, 1.0, stereo_mixed , zi=stereo_filt)

        #mono audio
        mono_audio_filt , state_lpf_16k = signal.lfilter(audio_coeff, 1.0, fm_demod, zi=state_lpf_16k)

        # downsample audio data
        stereo_audio_block = stereo_audio_filt[::audio_decim]
        mono_audio_block = mono_audio_filt [::audio_decim]

        audio_left_block = (mono_audio_block + stereo_audio_block) / 2
        audio_right_block =  (mono_audio_block - stereo_audio_block) / 2

        audio_left = np.concatenate((audio_left, audio_left_block))
        audio_right = np.concatenate((audio_right, audio_right_block))

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
        Stereo_data = np.concatenate((Stereo_data, stereo_audio_block))

        # to save runtime select the range of blocks to log data
		#this includes both saving binary files as well plotting PSD
		#below we assume we want to plot for graphs for blocks 10 and 11
        if block_count >= 10 and block_count < 12:
            # plot PSD of selected block after FM demodulation
            ax0.clear()
            fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')

            # output binary file name (where samples are written from Python)
            fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
            fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after downsampling mono audio
			# ... change as needed

			# save figure to file
            fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

        block_count += 1
    print('Finished processing all the blocks from the recorded I/Q samples')

    audio_data = np.vstack((audio_left, audio_right))
    audio_data = audio_data.transpose()

	# write audio data to file
    out_fname = "../data/fmStereo.wav"
    wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

    print(rf_coeff)

	# uncomment assuming you wish to show some plots
	# plt.show()
