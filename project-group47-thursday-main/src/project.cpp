/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/


#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "supportLib.h"


void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for (int k=0; k<(int)num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}


int main(int argc, char* argv[])
{

	// Default mode 0
	int mode = 0;

	// Mode Selection
	if (argc<2){
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if (argc==2){
		mode=atoi(argv[1]);
		if (mode>3){
			std::cerr << "Wrong mode " << mode << std::endl;
			exit(1);
		}
	} else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0-3" << argv[0] << std::endl;
		exit(1);
	}

	std::cerr << "Operating in mode " << mode << std::endl;

	// Defining mode based variables
	float Fs;
	int rf_decim;
	int audio_decim;
	int audio_up;
	int block_size;

	//Setting parameters based on mode of operation
	if (mode == 0){
		Fs = 2400000;
		rf_decim = 10;
		audio_decim = 5;
		audio_up = 1;
		block_size = 1024*rf_decim * ((float)audio_decim/audio_up) * 2;
	}
	else if (mode == 1){
		Fs = 1152000;
		rf_decim = 4;
		audio_decim = 6;
		audio_up = 1;
		block_size = 1024*rf_decim * ((float)audio_decim/audio_up) * 2;
	}
	else if (mode == 2){
		Fs = 2400000;
		rf_decim = 10;
		audio_up = 147;
		audio_decim = 800;
		block_size = 1029*rf_decim * ((float)audio_decim/audio_up) * 2;
	}
	else{
		Fs = 960000;
		rf_decim = 3;
		audio_up = 441;
		audio_decim = 3200;
		block_size = 1323*rf_decim * ((float)audio_decim/audio_up) * 2;
	}

	// GENERAL VARIABLES
	unsigned short int N_taps = 101; // Number of Taps for all
	std::vector<float> block_data(block_size);


	// FRONT END VARIABLES
	float rf_Fc = 100000;                    // Cutoff Freq for Front end
	std::vector<float> dataI;                // I and Q after splitting input data block
	std::vector<float> dataQ;
	std::vector<float> downI;                // I and Q after downsampling
	std::vector<float> downQ;
	std::vector<float> h;                    // h vector for front end and mono coeff
	std::vector<float> filteredI;            // I and Q after filtering
	std::vector<float> filteredQ;
	std::vector<float> filtIz(N_taps, 0.0);  // I and Q state saving when filtering
	std::vector<float> filtQz(N_taps, 0.0);
	float lastI = 0;  						           // State saving for demodulation in Front End
	float lastQ = 0;
	std::vector<float> interF;               // Intermediate Frequency


	// MONO VARIABLES
	float mono_Fc = 16000;                     // Cutoff Freq for Mono
	std::vector<float> IFfilt;                 // Filtered IF
	std::vector<float> filtMonoz(N_taps, 0.0); // Mono state saving when filtering
	std::vector<float> monoOut;                // Mono Output


	// STEREO VARIABLES
	float carrierFb;
	float carrierFe;
	std::vector<float> carrierH;
	std::vector<float> carrierState(N_taps, 0.0);
	std::vector<float> carrierFiltered;

	float extrFb;
	float extrFe;
	std::vector<float> extrH;
	std::vector<float> extrState(N_taps, 0.0);
	std::vector<float> extrFiltered;

	float PLLlockFreq = 19000;

	std::vector<float> ncoOut(block_size/(2*rf_decim), 1.0);
	float integrator = 0.0;
	float phaseEst = 0.0;
	float feedbackI = 1.0;
	float feedbackQ = 0.0;
	float trigOffset = 0;
	float ncoScale = 1;
	float phaseAdjust = 0;
	float normBand = 0.01;

	std::vector<float> stereoMixed(block_size/(2*rf_decim), 0.0);

	std::vector<float> stereoConverted;
	std::vector<float> stereoH;
	std::vector<float> stereoConvState(N_taps, 0.0);
	std::vector<float> strOut;

	std::vector<float> stereoLR(block_size/(rf_decim*((float)audio_decim/audio_up)), 0.0);

	for (unsigned int block_id = 0; ;block_id++){


		// ----- READING -----
		readStdinBlockData(block_size, block_id, block_data);
		if (std::cin.rdstate()!=0){
			std::cerr << "End of Input Stream Reached" << std::endl;
			exit(1);
		}

		std::cerr << "Read block" << block_id << std::endl;



		// ----- RF FRONT END -----
		//std::cerr << "block_size: " << block_size << std::endl;
		//std::cerr << "block_data SIZE: " << block_data.size() << std::endl;

		// Cut Block into I and Q
		splitIQ(block_data, dataI, dataQ);

		//std::cerr << "dataI SIZE: " << dataI.size() << std::endl;

		// Pass Through Low Pass filter
		impulseResponseLPF(Fs, rf_Fc, N_taps, 1, h);
		convolveOLD(filteredI, dataI, h, filtIz, rf_decim);
		convolveOLD(filteredQ, dataQ, h, filtQz, rf_decim);

		//std::cerr << "filteredI SIZE: " << filteredI.size() << std::endl;

		// Downsample
		downsample(filteredI, downI, rf_decim);
		downsample(filteredQ, downQ, rf_decim);

		//std::cerr << "downI SIZE: " << downI.size() << std::endl;

		// Demodulate
		fmDemodulation(downI, downQ, lastI, lastQ, interF);

		//std::cerr << "interF SIZE: " << interF.size() << std::endl;



		// ----- STEREO -----

		// Carrier BP Filter
		carrierFb = 18500;
		carrierFe = 19500;
		impulseResponseBPF(Fs/rf_decim, carrierFb, carrierFe, N_taps, carrierH);
		convolveOLD(carrierFiltered, interF, carrierH, carrierState, 1);

		//std::cerr << "carrierFiltered SIZE: " << carrierFiltered.size() << std::endl;

		// Extraction BP Filter
		extrFb = 22000;
		extrFe = 54000;
		impulseResponseBPF(Fs/rf_decim, extrFb, extrFe, N_taps, extrH);
		convolveOLD(extrFiltered, interF, extrH, extrState, 1);

		//std::cerr << "extrFiltered SIZE: " << extrFiltered.size() << std::endl;

		// PLL and NCO
		fmPLL(carrierFiltered, PLLlockFreq, Fs/rf_decim, ncoOut, integrator, phaseEst, feedbackI, feedbackQ, trigOffset, ncoScale, phaseAdjust, normBand);

		//std::cerr << "ncoOut SIZE: " << ncoOut.size() << std::endl;

		// Mixing
		for (int i = 0; i < extrFiltered.size(); i++) {
				stereoMixed[i] = extrFiltered[i]*ncoOut[i]*2;
		}
		//std::cerr << "STEREO MIXED SIZE: " << stereoMixed.size() << std::endl;



		// ----- MONO PATH -----
		if (mode == 2 or mode == 3){
			// MONO
			impulseResponseLPF((Fs/rf_decim)*audio_up, mono_Fc, audio_up * N_taps, audio_up, h);
			convolveFIR(IFfilt, interF, h, filtMonoz, audio_decim, audio_up);
			monoOut = IFfilt;

			// STEREO
			impulseResponseLPF((Fs/rf_decim)*audio_up, mono_Fc, audio_up * N_taps, audio_up, stereoH);
			convolveFIR(stereoConverted, stereoMixed, stereoH, stereoConvState, audio_decim, audio_up);
			strOut = stereoConverted;
		}
		else{
			// MONO
			impulseResponseLPF(Fs/rf_decim, mono_Fc, N_taps, 1, h);
			convolveOLD(IFfilt, interF, h, filtMonoz, audio_decim);
			downsample(IFfilt, monoOut, audio_decim);

			// STEREO
			impulseResponseLPF(Fs/rf_decim, mono_Fc, N_taps, 1, stereoH);
			convolveOLD(stereoConverted, stereoMixed, stereoH, stereoConvState, audio_decim);
			downsample(stereoConverted, strOut, audio_decim);
		}
		//std::cerr << "STR OUT SIZE: " << strOut.size() << std::endl;
		//std::cerr << "MONO OUT SIZE: " << monoOut.size() << std::endl;


		// STEREO COMBINATION
		for (int j = 0; j < strOut.size(); j++){
			stereoLR[2*j] = (monoOut[j] + strOut[j])/2;
			stereoLR[2*j+1] = (monoOut[j] - strOut[j])/2;
		}
		//std::cerr << "stereoLR SIZE: " << stereoLR.size() << std::endl;


		// // ----- OUTPUT FOR STEREO -----
		// std::vector<short int> audio_data(stereoLR.size());
		// for (unsigned int k=0; k<stereoLR.size(); k++){
		// 	if (std::isnan(stereoLR[k])) audio_data[k] = 0;
		// 	else audio_data[k] = static_cast<short int> (stereoLR[k]*16384);
		// }
		//
		// fwrite(&audio_data[0], sizeof(short int), audio_data.size()-1, stdout);

		// ----- OUTPUT FOR MONO -----
		std::vector<short int> audio_data(monoOut.size());
		for (unsigned int k=0; k<monoOut.size(); k++){
			if (std::isnan(monoOut[k])) audio_data[k] = 0;
			else audio_data[k] = static_cast<short int> (monoOut[k]*16384);
		}

		fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);

	}

	return 0;
}
