/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, const int &gain, std::vector<float> &h);
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &down, const int &up);
void convolveOLD(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &zi, const int &step);
void impulseResponseBPF(const float Fs, const float Fb, const float Fe, const unsigned short int num_taps, std::vector<float> &h);

#endif // DY4_FILTER_H
