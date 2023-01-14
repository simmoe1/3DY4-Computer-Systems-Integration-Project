
#include "dy4.h"
#include "supportLib.h"

//function to split our array of inputs according to block size
void cutArray(const int position1, const int position2, const std::vector<float> &filtered_data, std::vector<float> &block_signal){
	int size = position2 - position1;
	block_signal.clear();
	block_signal.resize(size);

	for (int i = 0; i < size; i++){
		block_signal[i] = filtered_data[position1 + i];
	}

}

//function to split a block of data into its I and Q values
void splitIQ(const std::vector<float> &inData, std::vector<float> &I, std::vector<float> &Q){

  I.clear();
	I.resize(inData.size()/2);

  Q.clear();
	Q.resize(inData.size()/2);

	int index = 0;
  for (int i=0; i<inData.size(); i+=2) {
    I[index] = inData[i];
		Q[index] = inData[i+1];
		index += 1;
  }

}

//demod function, adapted from python
void fmDemodulation(const std::vector<float> &I, const std::vector<float> &Q, float &prev_I, float &prev_Q, std::vector<float> &fm_demod) {
	//allocating memory for output
	fm_demod.clear();
	fm_demod.resize(I.size());

	for (int k = 0; k < I.size(); k++) {
      if (I[k] == 0 or Q[k] == 0.0) {
          fm_demod[k] = 0.0;
      } else {
          fm_demod[k] = (I[k] * (Q[k] - prev_Q) - (Q[k] * (I[k] - prev_I))) / ((I[k] * I[k]) + (Q[k] * Q[k]));
			}
			//state saving vars
      prev_I = I[k];
      prev_Q = Q[k];
  }

}

//simple downsample function
void downsample(const std::vector<float> &in, std::vector<float> &out, int &val){
	//determing size needed after downsampling so we can allocate proper memory
	std::vector<float>::size_type s = in.size()/val;
	//allocating new memory for our output
	out.clear();
	out.resize(s);

	int x = 0;
	for (int i=0; i<in.size(); i+=val) {
		out[x] = in[i];
		x += 1;
	}

}

//simple upsampling function
void upsample(const std::vector<float> &in, std::vector<float> &out, int val){

	int x = 0;
	out.clear();
	out.resize(in.size()*val, 0.0);
	for (int i=0; i<in.size()*val; i+=val) {
		out[i] = in[x];
		x += 1;
	}
}

//PLL function taken from lecture and edited to include state saving
void fmPLL(std::vector<float> pllIn, const float freq, const float Fs, std::vector<float> &ncoOut, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &trigOffset, const float ncoScale, const float phaseAdjust, const float normBand){
	//scale factors for proportional/integrator terms
	float Cp = 2.666;
	float Ci = 3.555;
	float errI, errQ, errD, trigArg;
	//gain for the proportional term
	float Kp = normBand*Cp;
	//gain for the integrator term
	float Ki = (normBand*normBand)*Ci;

	for(int k=0; k < pllIn.size(); k++){
		//phase detector
		errI = pllIn[k]*(+feedbackI);
		errQ = pllIn[k]*(-feedbackQ);


		//arctan function for phase err detection
		errD = std::atan2(errQ,errI);

		//loop filter and updating phase estimate
		integrator = integrator + Ki*errD;
		phaseEst += (Kp*errD) + integrator;

		//oscillator
		trigOffset += 1;
		trigArg = 2*PI*(freq/Fs)*(trigOffset)+phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos((trigArg*ncoScale));

	}
	ncoOut[0] = ncoOut[ncoOut.size()-1];
}
