#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <algorithm> 
#include "DefineClass.h"
#include "DefineParameter.h"
#include "AstarDecode.h"
#include "LinearBlockCodes.h"

using namespace std;

// Polar Code Encoder
void Polar_Code_G_Generator(MATRIX<__int8> &G_, vector<int>&frozen, vector<int>&non_frozen, vector<double> &Channel_parameter);

// Polar + Astar
void Hybrid_Polar_Astar_Encoder(
	MATRIX<__int8> &G, 
	MATRIX<__int8> &G_Inner, 
	vector<__int8> &message,
	vector<__int8> &codeword,
	DECODING_INFO &decoding_info);

void Hybrid_Polar_Astar_Decoder(
	MATRIX<__int8> &G,
	MATRIX<__int8> &G_Inner,
	DECODING_INFO &decoding_info);

//Polar Code Decoder
void Polar_code_BP_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void Inner_Polar_BP_Decoder_Outer_Astar_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// Some algorithms
double min_sum(double a, double b);

// RM Decoding Algorithms
void Majority_Soft_decoding(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void Majority_Soft_decoding(MATRIX<__int8> &G, vector<double> &Rx, double var);
void Majority_Soft_decoding_ver2(MATRIX<__int8> &G, MATRIX<__int8> &G_, vector<double> &Rx, double var);
void Majority_Soft_decoding_ver3(MATRIX<__int8> &G, MATRIX<__int8> &G_, vector<double> &Rx, double var);

void Majority_Hard_decoding(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void Majority_Hard_decoding(MATRIX<__int8> &G, MATRIX<__int8> &G_, DECODING_INFO &decoding_info);

double RM_LLR_order2(double L1, double L2, double L3, double L4);
double RM_LLR_order1(double L1, double L2);


/*
class POLAR_FUNC {

public:
	vector<int> frozen;
	vector<int> non_frozen;
	//friend CODE;
};*/