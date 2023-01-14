
#ifndef DY4_SUPPORTLIB_H
#define DY4_SUPPORTLIB_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

void cutArray(const int position1, const int position2, const std::vector<float> &filtered_data, std::vector<float> &block_signal);
void splitIQ(const std::vector<float> &inData, std::vector<float> &I, std::vector<float> &Q);
void fmDemodulation(const std::vector<float> &I, const std::vector<float> &Q, float &prev_I, float &prev_Q, std::vector<float> &fm_demod);
void downsample(const std::vector<float> &in, std::vector<float> &out, int &val);
void upsample(const std::vector<float> &in, std::vector<float> &out, int val);
void fmPLL(std::vector<float> pllIn, const float freq, const float Fs, std::vector<float> &ncoOut, float &integrator, float &phaseEst, float &feedbackI, float &feedbackQ, float &trigOffset, const float ncoScale, const float phaseAdjust, const float normBand);

#endif // DY4_LOGFUNC_H
