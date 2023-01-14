/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <valarray>

#define PI 3.14159265358979323846


// function for DFT (reused from previous experiment)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf)
{
	Xf.clear();
	Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (unsigned int m = 0; m < Xf.size(); m++) {
		for (unsigned int k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to generate a sine with N samples per second over interval
void generateSin(std::vector<float> &t, std::vector<float> &x, float Fs, float interval, float frequency = 7.0, float amplitude = 5.0, float phase = 0.0)
{
	// we do NOT allocate memory space explicitly
	// for the time (t) vector and sample (x) vector
	t.clear(); x.clear();
	float dt = 1/Fs;
	for (float i = 0.0; i < interval; i += dt) {
		// vector size increases when pushing new elements into it
		t.push_back(i);
		x.push_back(amplitude*std::sin(2*PI*frequency*i+phase));
	}
}

// function to add an array of sines
void addSin(const std::vector<std::vector<float>> &sv, std::vector<float> &added)
{
	// assumes at least one sine passed to this function
	// assumes all input sines are of the same size
	for (unsigned int i = 0.0; i < sv[0].size(); i ++) {
		float addval = 0.0;
		// note: sv.size() returns the number of sines (or rows in 2D repr)
		// sv[0].size() returns the number of samples in a sine (or cols in 2D repr)
		for (unsigned int k = 0; k < sv.size(); k++)
			addval += sv[k][i];
		added.push_back(addval);
	}
}

// function to print a real vector (reused from previous experiment)
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (unsigned int i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

// function to print a complex vector (reused from previous experiment)
void printComplexVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (unsigned int i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

// function to record data in a format to be read by GNU plot
// the arguments are VERY specific to this usage in this experiment
// we have the time vector (t), a vector of sines (sv),
// input samples (x, i.e., added sines for this experiment),
// output samples (y, filtered input samples)
// frequency vectors (both Xf and Yf)
// note: the reference code does NOT do filtering, hence y and Yf are zeros by default
void plotaddedSinesSpectrum(const std::vector<float> &t, const std::vector<std::vector<float>> &sv, const std::vector<float> &x, const std::vector<float> &y, const std::vector<std::complex<float>> &Xf, const std::vector<std::complex<float>> &Yf)
{
	// write data in text format to be parsed by GNU plot
	const std::string filename = "../data/example.dat";
	std::fstream fd;  // file descriptor
	fd.open(filename, std::ios::out);
	fd << "#\tindex\tsine(0)\tsine(1)\tsine(2)\tdata in\tdata out\tspectrum in\tspectrum out\n";
	for (unsigned int i = 0; i < t.size(); i++) {
		fd << "\t " << i << "\t";
		for (unsigned int k = 0; k < sv.size(); k++)
			fd << std::fixed << std::setprecision(3) << sv[k][i] << "\t ";
		fd << x[i] << "\t "<< y[i] << "\t\t ";
		fd << std::abs(Xf[i])/Xf.size() << "\t\t " << std::abs(Yf[i])/Yf.size() <<"\n";
	}
	std::cout << "Generated " << filename << " to be used by GNU plot\n";
	fd.close();
}

// function to compute the impulse response "h" based on the sinc function
// see pseudocode from previous lab that was implemented in Python
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
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (auto i = 0; i < y.size(); i++) {
		y[i] = 0;
		for (auto m = 0; m < h.size(); m++) {
			if (m <= i){
				y[i] = y[i] + x[i-m] * h[m];
			}else{
				y[i] = 0;
			}
		}
	}
}

//FUNCTION FOR SQUARE WAVE
void squareWave (float Fs, float interval, int frequency, int duty, std::vector<float> &squareWaveX)
{
	float dT = 1.0/Fs;

	//time = np.arange(0, interval, dT)

	std::vector<std::complex<float>> timeX;
	timeX.clear(); timeX.resize(interval, 0.0);
	int j=0;
	for (int i = 0; i < interval; i = i + dT){
		timeX[i] += i;
		j+=1;
	}

	float period = 1/frequency;
	float periodSamples = period/dT;
	float positiveFloatSamples = periodSamples*duty;

	//std::vector<std::complex<float>> x;
	squareWaveX.clear(); squareWaveX.resize(timeX.size(), 0.0);

	while (int i = 0 < int(timeX.size()/periodSamples)){
		for (int j = (i+1)*periodSamples-positiveFloatSamples; j < (i+1)*periodSamples; j+=1){
			squareWaveX[j] = 1;
			j+=1;
		}

	}
}

int main()
{

	float Fs = 1024.0;                   // samples per second
	float interval = 1.0;                // number of seconds
	unsigned short int num_taps = 101;   // number of filter taps
	float Fc = 30.0;                     // cutoff frequency (in Hz)

	// declare a vector of vectors for multiple sines
	std::vector<std::vector<float>> sv;
	// declare time and sine vectors
	std::vector<float> t, sine;
	// note: there is no explicit memory allocation through vector resizing
	// vector memory space will increase via the push_back method

	// generate and store the first tone
	// check the function to understand the order of arguments
	//generateSin(t, sine, Fs, interval, 10.0, 5.0, 0.0);
	//sv.push_back(sine);
	// generate and store the second tone
	//generateSin(t, sine, Fs, interval, 60.0, 2.0, 0.0);
	//sv.push_back(sine);
	// generate and store the third tone
	//generateSin(t, sine, Fs, interval, 80.0, 3.0, 0.0);
	//sv.push_back(sine);

	std::vector<float> x;
	std::vector<float> squareWaveX;
	squareWave(Fs, interval,rand() % 15 + 5, 0.5, squareWaveX);
	sv.push_back(squareWaveX);
	sv.push_back(squareWaveX);
	sv.push_back(squareWaveX);

	// declare the added sine vector and add the three tones
	//addSin(sv, x);
	// printRealVector(x);

	// declare a vector of complex values for DFT; no memory is allocated for it
	std::vector<std::complex<float>> Xf;
	DFT(x, Xf);
	// printComplexVector(Xf);

	// generate the impulse response h
	// convolve it with the input data x
	// in order to produce the output data y
	std::vector<float> h;              // impulse response
	impulseResponseLPF(Fs, Fc, num_taps, h);
	std::vector<float> y;              // filter out
	convolveFIR(y, x, h);

	// compute DFT of the filtered data
	std::vector<std::complex<float>> Yf;
	DFT(y, Yf);




	// prepare the data for GNU plot
	plotaddedSinesSpectrum(t, sv, x, y, Xf, Yf);


	// naturally, you can comment the line below once you are comfortable to run GNU plot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
