#pragma once
#ifndef _AWGN_
#define _AWGN_

#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <random>
#include "DefineParameter.h"
#include <algorithm> // for transform
#include <functional> // for plus

using namespace std;

unsigned int rnd(void);
double Uniform(void);
double Gaussian_Variable_Original(double mean, double variance);

double Gaussian_Variable_New(double mean, double variance);
void AWGN_Channel(
	double mean,
	double snr_dB,
	double code_rate,
	vector<double> &tx_signal_seq,
	vector<double> &rx_signal_seq,
	bool channel_info_flag);
#endif // !_AWGN_