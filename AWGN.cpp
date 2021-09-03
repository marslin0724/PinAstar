
#include "PinAstar/AWGN.h"
using namespace std;

static unsigned int Num = 0;
// AWGN Channel

unsigned int rnd()
{
	unsigned int x;
	x = rand();
	Num++;
	if (Num >= RAND_MAX){
		time_t t;
		t = time(NULL);
		srand((unsigned)(t % RAND_MAX));
		Num = 0;
	}
	return x;
} // end rnd()

double Uniform()
{
	return (double)(rnd() + 1) / (RAND_MAX + 2);
} //end Uniform()

double Gaussian_Variable_Original(double mean, double variance)
{
	double UnifromVR1((double)Uniform()), UnifromVR2((double)Uniform());
	return (double)(mean + (sqrt(-2. * log(UnifromVR1)) * cos(2. * PI * UnifromVR2) * sqrt(variance)));
} // end Gaussian_Variable()

double Gaussian_Variable_New(double mean, double variance)
{
	std::random_device rd;
	std::mt19937 generator(rd());
	std::normal_distribution<double> norm(mean, sqrt(variance));
	return norm(generator);
}

void AWGN_Channel(
	double mean, 
	double snr_dB, 
	double code_rate, 
	vector<double> &tx_signal_seq , 
	vector<double> &rx_signal_seq ,
	bool channel_info_flag)
{
	//code_rate = 1;
	//cout << code_rate << endl;
	double variance((double)((Es / (2.*code_rate)) * pow(10., (-(snr_dB / 10.)))));		// SNR [dB] ¡÷ SNR ¡÷ N_0/2
	//cout << variance <<endl;
	//size_t sequence_length(tx_signal_seq.size());
	//for (size_t i(0); i < sequence_length; ++i)
	//	rx_signal_seq.at(i) = tx_signal_seq.at(i) + Gaussian_Variable_New(mean, variance);
		//Gaussian_Variable_Original(mean, variance); //Original methods
	
	std::random_device rd;
	std::mt19937 generator(rd());
	std::normal_distribution<double> norm(mean, sqrt(variance));


	// generator a gaussian random vector
	auto gen = [&norm, &generator]() {
		return norm(generator);
	};
	vector<double> Gussian_vec(tx_signal_seq);
	generate(begin(Gussian_vec), end(Gussian_vec), gen);


	// addition of vector
	std::transform(
		tx_signal_seq.begin(),	//input1
		tx_signal_seq.end(),	
		Gussian_vec.begin(),	//input2
		rx_signal_seq.begin(),	//output
		std::plus<double>());	//operation type



	//Show AWGN channel information
	if (channel_info_flag == TRUE){
		cout
			<< "\n\n --------- Channel Information ---------\n\n"
			<< " Code Rate : " << code_rate << "\n"
			<< " Es : " << Es << " [ joule / per symbol ] \n"
			<< " Eb : " << Es * (1. / code_rate) << " [ joule / per info. bit ] \n"
			<< " SNR_dB : " << snr_dB << " dB \n"
			<< " Gaussian Mean : " << mean << "\n"
			<< " Gaussian Variance : " << variance << "\n"
			<< "\n --------- Channel information ---------\n\n";
	}
} // end AWGN_Channel