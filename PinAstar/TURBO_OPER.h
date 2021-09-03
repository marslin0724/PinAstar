#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <algorithm> 
#include <math.h>
#include "DefineClass.h"
#include "DefineParameter.h"
#include "AstarDecode.h"
#include "LinearBlockCodes.h"
#include "LDPC_OPER.h"

void Turbo_LDPC_encoder(MATRIX<__int8> &G_, MATRIX<__int8> &G, vector<__int8> &message_seq, vector<__int8> &output_codeword_seq);
void Turbo_LDPC_decoder(MATRIX<__int8> &H_, MATRIX<__int8> &H, MATRIX<__int8> &G, DECODING_INFO &decoding_info, vector<int> &permutation_seq_, vector<int> &permutation_seq);
void Turbo_LDPC_Repitition_Method(MATRIX<__int8> &H_,
	MATRIX<__int8> &H, MATRIX<__int8> &G,
	DECODING_INFO &decoding_info,
	vector<int> &permutation_seq_,
	vector<int> &permutation_seq);
void Turbo_LDPC_decoder_punctured(MATRIX<__int8> &H_, MATRIX<__int8> &H, DECODING_INFO &decoding_info, vector<int> &permutation_seq_, vector<int> &permutation_seq);
// Use Two LDPC Code To SISO
void Turbo_Encoder_Punctured(vector<__int8> &message_seq, vector<__int8> &output_codeword_seq, int row, int col);  // Rate =1 /2
void Turbo_LDPC_encoder_Punctured_Version(MATRIX<__int8> &G, vector<__int8> &message_seq, vector<__int8> &output_codeword_seq);

// Some Tools
void Remove_Even_Index(vector<__int8> &Vec);
void Remove_Odd_Index(vector<__int8> &Vec);
void Zero_Padding_Even(vector<double> &Vec);
void Zero_Padding_Odd(vector<double> &Vec);

void Block_Interleaver(vector<__int8> &input);
void Block_Interleaver(vector<double> &input);
void Block_DeInterleaver(vector<__int8> &input);
void Block_DeInterleaver(vector<double> &input);