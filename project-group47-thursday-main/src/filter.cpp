/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, const int &gain, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear();
	h.resize(num_taps, 0.0);

	float normCutoff = Fc/(Fs/2);

	for (int i = 0; i < num_taps; i++) {
		if (i == (num_taps-1)/2){
			h[i] = normCutoff;
		}
		else{
			h[i] = normCutoff * (std::sin(PI*normCutoff*(i-(num_taps-1)/2)))/(PI*normCutoff*(i-(num_taps-1)/2));
		}
		//state saving for next block as well as incorporating gain
		h[i] = h[i]*std::pow(std::sin((i*PI)/num_taps), 2);
		h[i] = h[i]*gain;

	}

}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
//also downsamples and upsamples based on the fast implementation
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &down, const int &up)
{
		// allocate memory for the impulse response
		y.clear();
		y.resize(x.size()*up/down, 0.0);
		int phase = 0;
		int fix = 0;
		int temp = 0;

		for (int i = 0; i < y.size(); i++){
			// compute phase
			temp = i * down;
			phase = (down * i)%up;

			// compute fix
			fix = ((down * i)-(phase)) / up;

			for (int k = phase; k < h.size();k+=up){
				if (temp-k >= 0){
					y[i] += h[k] * x[fix];
				}
				else{
					y[i] += h[k] * zi[zi.size()+fix];
				}
				fix--;
			}
		}
		//state saving for next block
		int index = x.size() - h.size()/up + 1;
		zi = std::vector<float>(x.begin() + index, x.end());
}

void convolveOLD(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &step)
{
	// allocate memory for the impulse response
	y.clear();
	y.resize(x.size(), 0.0);
	int xIndex;
	int taps = h.size();

	for (int i = 0; i < x.size(); i+=step){
		y[i] = 0;
		xIndex = i - taps + 1;
		for (int k = 0; k < taps; k++){
			if (k + xIndex < 0) y[i] += zi[k+xIndex+taps]*h[taps-k-1];
			else y[i] += x[k+xIndex]*h[taps-k-1];
		}
	}
	//state saving
	for (int j = 0; j < taps; j++){
		zi[j] = x[x.size() + j - taps];
	}
}

//bandpass impulse response filter, taken from python bandpass
void impulseResponseBPF(const float Fs, const float Fb, const float Fe, const unsigned short int num_taps, std::vector<float> &h)
{
	//allocating memory for new output
	h.clear();
	h.resize(num_taps, 0.0);

 	float temp;

	float normCenter = ((Fe + Fb)/2)/(Fs/2);
	float normPass = (Fe-Fb)/(Fs/2);

	for (int i = 0; i < num_taps; i++){
		if (i == (num_taps-1)/2){
			h[i] = normPass;
		}
		else{
			temp = PI*(normPass/2)*(i-(num_taps-1)/2);
			h[i] = normPass*((std::sin(temp))/(temp));
		}
		h[i] = h[i]*std::cos(i*PI*normCenter);
		h[i] = h[i]*std::pow(std::sin((i*PI)/num_taps), 2);
	}
}
