/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (unsigned int m = 0; m < Xf.size(); m++) {
		for (unsigned int k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
  for (unsigned int i = 0; i < Xf.size(); i++) {
    Xmag[i] = std::abs(Xf[i])/Xf.size();
  }
}

// add your own code to estimate the PSD
void estimatePSD(const std::vector<float> &samples, double Fs, std::vector<float> &freq, std::vector<float> &psd_estimate)
{
	const auto freq_bins = NFFT;
	const auto df = Fs / freq_bins;

	freq.resize(freq_bins / 2 + (freq_bins % 2 != 0), 0);//ceiling for integer division

	std::vector<float> hann(freq_bins);
	for (auto i = 0; i < hann.size(); i++) hann[i] = std::pow(std::sin((i * PI) / freq_bins), 2); //design hann window

	//create an empty list where the PSD for each segment is computed
	std::vector<float> psd_list;

	//samples should be a multiple of frequency bins, so
	//the number of segments used for estimation is an integer
	//note: for this to work you must provide an argument for the
	//number of frequency bins not greater than the number of samples!
	auto no_segments = int(floor(samples.size()) / float(freq_bins));

	//iterate through all the segments
	for (auto k = 0; k < no_segments; k++) {
		std::vector<float> windowed_samples(freq_bins);
		for (auto i = 0; i < freq_bins; i++) {
			windowed_samples[i] = samples[k * freq_bins + i] * hann[i];
		}
		//compute fourier transform
		std::vector<std::complex<float>> Xf;
		DFT(windowed_samples, Xf);

		//only keep positive half of the spectrum.to have more accurate PSD estimate when plotting
		Xf = std::vector<std::complex<float>>(Xf.begin(), Xf.begin() + (freq_bins/2));

		float psd_seg;
		for (auto i = 0; i < Xf.size(); i++) {
			psd_seg = 1.0 / (Fs * freq_bins / 2) * std::pow(std::abs(Xf[i]), 2);
			psd_seg = 2.0 * psd_seg;
			psd_seg = 10.0 * std::log10(psd_seg);
			psd_list.push_back(psd_seg);
		}
	}

	//compute estmate to be returned by function through averaging
	psd_estimate.resize(int(freq_bins/2), 0);

	//iterate through all the frequency bins, average them
	for (auto k = 0; k < (freq_bins/2); k++) {
		for (auto l = 0; l < no_segments; l++) {
			psd_estimate[k] += psd_list[k + l * (freq_bins/2)];
		}
		psd_estimate[k] = psd_estimate[k] / no_segments;
	}
}
