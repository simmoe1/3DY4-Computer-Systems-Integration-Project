/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "math.h"
#include "numeric"
#include "algorithm"

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float normCut = Fc / (Fs/2);
	for (auto i = 0; i < num_taps; i++){
		if (i==(num_taps-1)/2){
			h[i]=normCut;
		}
		else {
			float variable = PI * normCut * (i-(num_taps-1)/2);
			h[i]=normCut * (std::sin(variable)/variable);
		}
		h[i]=h[i]* std::sin((i * PI) / num_taps) * std::sin((i * PI) / num_taps);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for (auto i = 0; i < y.size(); i++) {
			y[i] = 0;
			for (auto m = 0; m < h.size(); m++) {
				if (m <= i){
					y[i] = y[i] + x[i-m] * h[m];
				}
				else{
					y[i] = 0;
			}
		}

	}
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
}
