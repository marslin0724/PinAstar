#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <fstream>
#include "DefineParameter.h"

#ifndef _CLASS_
#define _CLASS_
using namespace std;
void Present_Large_Number(__int64 input_number);

// Class
template<class T >
class MATRIX
{
private:
	void Delete_Matrix(void);
public:
	std::vector< std::vector<T> > _matrix;
	std::vector< std::vector<T> > _matrix_inner;
	std::vector< std::vector<T> > _matrix_outer;
	std::vector< std::vector<T> > _H;
	std::vector< std::vector<T> > subG;
	std::vector<int> interleaver;
	size_t Col_number = 0, Row_number = 0;
	
	// Functions for class MATRIX
	void Building_Empty_Matrix(size_t row, size_t col);
	void Show_matrix(size_t x);

	// Constructor
	MATRIX(size_t row, size_t col);
	MATRIX();

	// Destructor
	~MATRIX();

	//PoHan Parity Check Matrix
	size_t Num_Low_Density_Row = 0;
	vector <size_t> Low_Density_Row_Index;
};

// Constructors
template<class T>
MATRIX<T>::MATRIX(size_t row, size_t col)
{
	Building_Empty_Matrix(row, col);
}

template<class T>
MATRIX<T>::MATRIX()
{
	// Empty constructor
}

// Destructor
template<class T>
MATRIX<T>::~MATRIX()
{
	Delete_Matrix();
}

template<class T>
void MATRIX<T>::Building_Empty_Matrix(size_t row, size_t col)
{
	Row_number = row;
	Col_number = col;
	std::vector<T> v(col, 0);
	_matrix.assign(row, v);
	v.~vector();
}	

template<class T>
void MATRIX<T>::Delete_Matrix(void)
{
	_matrix.~vector();
}	

template<class T>
void MATRIX<T>::Show_matrix(size_t BlankSpace)
{
	for (size_t i(0); i < Row_number; ++i){
		for (size_t j(0); j < Col_number; ++j)
			std::cout << std::setw(BlankSpace) << (int)_matrix[i][j];
		std::cout << std::endl;
	}
}


class NODE_PATH
{
private:
public:
	double metric = 0, heuristic = 0, fano_metric = 0;
	int parity_index = 0, non_updated_num = 0;
	int base = 0;
	//segment
	int temp_idx = 0, idx_end = 0, segment = 0;
	double segment_metric = 0;
	vector <__int8>
		// corresponding message bits
		message_bits;
		// This is used to store the positions of bits which differ from hard-decision z
		// for speeding up the rate of algorithms. We can save the number of encoding during
		// exeution of tree search through this way.
	vector <__int16>
		Diff_Index,
		Same_Index;

	// D_z  : corresponding hamming distance compare with hard-decision z in MRIP
	// D_zp : corresponding hamming distance compare with hard-decision z in CB
	size_t level = 0, D_z = 0, D_zp = 0;
	
	NODE_PATH(size_t message_length) { message_bits.assign(message_length, 0); }
	NODE_PATH(void) {}

	double auxiliary(void){
		return (metric + heuristic);
	}
};

class CODE
{
private:
	void ChooseCode();

public:
	double Code_Rate;
	size_t Code_number;
	size_t Row_Number, Col_Number;
	size_t Inner_Row_Number, Inner_Col_Number;
	size_t Short_Inner_Row_Number, Short_Inner_Col_Number;
	string Title;
	MATRIX<__int8> G,G_,G_Inner,H,H_;
	
	//LDPC
	vector<int> Permutation_Seq, Permutation_Seq_;
	int Punctured_Number;

	//Polar
	vector<int> frozen;
	vector<int> non_frozen;
	vector<double> Channel_Parameter;

	// Yin: Part of CRC
	vector<__int8> CRC;
	void CRC_Desicion(vector<__int8>& CRC);
	vector<__int8> CRC_Generator(vector<__int8> Message, vector<__int8> variant); // Return the parity bits added to the Message as CRC
	bool CRC_Examination(vector<__int8> Message_and_CRC, vector<__int8> variant);         // Return True if "Message" fits the CRC condition ("variant") 
	// Yin
	

	CODE() { 
		if (CRC_Check == TRUE) { CRC_Desicion(CRC); }
		ChooseCode(); 
	}
	//friend POLAR_FUNC;
};

class OPER_PARA
{
public:
	
	double
		SNR_dB_start = 0,
		SNR_dB_end	 = 0,
		SNR_dB_step  = 0,
		Alpha        = 0, //  for A*-Parity-f 
		OSC_Alpha    = 0,
		Worst_metric_ratio = 0,
		Constraint_i_ratio = 0,
		// Fano metric
		Fano_Metric_Parameter = 0,
		New_OSC_Alpha = 0; // for OSCII
	size_t
		// Es or Eb/N0
		Eb_N0_1_or_Es_N0_2 = 0,

		StackSize = 0,
		ErrorBlock_Thr = 0,
		BlockStep = 0,
		Decoder_Version = 0,
		//SPA
		SPA_I = 0,
		// Memory Reduction
		Constraint_i = 0,	// MRIP constraint

		// Deleting Mechanism
		// DM-I 
		Constraint_j = 0,	// ControlBand constraint 
		Control_Level = 0,

		// DM-II
		DM_StackSize = 0,
		Check_Level = 0,

		//Multi-Stack 
		Multi_Stack_Size = 0,

		// Parity Check
		Parity_Check_bits = 0,

		// CBC - Filpping Bit
		CBC_FlippingBit = 0,

		// List Decoder
		SCL_L = 0,
		Parameter_ReadTxt_Flag = FALSE,
		//section
		section1_i = 0,
		section2_i = 0,
		N_section;
			

	string Algorithm;
	size_t block_number = 0;

	OPER_PARA(CODE Code);
	void ShowDecoder();
	void DecodingParameterKeyIn();
	void DecodingParameterShow();
};

/////////////// Worst case
class DECODING_INFO
{
public:
	//Yin
	/*
	double Dz_1_Max_K = 0, Dz_2_Max_K = 0, Dz_3_Max_K = 0, Dz_4_Max_K = 0, Dz_5_Max_K = 0;
	double Dz_1_Tot_K = 0, Dz_2_Tot_K = 0, Dz_3_Tot_K = 0, Dz_4_Tot_K = 0, Dz_5_Tot_K = 0;
	int Dz_1_Count = 0, Dz_2_Count = 0, Dz_3_Count = 0, Dz_4_Count = 0, Dz_5_Count = 0;
	int Dz_Break = 0;
	*/
	// Turbo
	vector<double> LLR_1_priori, LLR_2_priori, LLR_1_extrinsic, LLR_2_extrinsic;
	int Turbo_Index;
	MATRIX<double> Iterative_Decoding_Candidates;
	MATRIX<double> Amplitude_Flipping_number;
	MATRIX<__int8> Iter_Decode_Candidates_HardDecision, Iter_Decode_Candidates_Est;
	vector<int> Error_Accumulation;
	bool EarlyTerminate;
	vector<double> Candidate_Metric;
	//Pin test
	size_t err_count[128];

	long double Ave_LLR;
	int Turbo_Small_Value_Counter;
	double ML_metric;
	double Best_metric;
	vector<__int8> Inner_message;
	//Polar 
	vector<int> frozen;
	vector<int> non_frozen;
	vector<double> Channel_Parameter;

	//LDPC
	bool Return_Est_1_or_LLR_0;

	int Code_Number;

	int Adaptive_i_Iteration;
	int Cancelled_Candidate_i;
	int Dz_0_Number = 0, Dz_1_Number = 0, Dz_2_Number = 0, Dz_3_Number = 0, Dz_4_Number = 0, Dz_5_Number = 0;
	int Dz_6_Number = 0, Dz_7_Number = 0, Dz_8_Number = 0, Dz_9_Number = 0, Dz_10_Number = 0, Dz_11_Number = 0;
	int Dz_12_Number = 0, Dz_13_Number = 0, Dz_14_Number = 0, Dz_15_Number = 0, Dz_16_Number = 0, Dz_17_Number = 0;
	int Dz_18_Number = 0, Dz_19_Number = 0, Dz_20_Number = 0, Dz_21_Number = 0, Dz_22_Number = 0, Dz_23_Number = 0;
	int D1 = 96, D2 = 96, D3 = 96, D4 = 96, D5 = 96;
	
	double Dz_1_Max_OSCr = 0, Dz_2_Max_OSCr = 0, Dz_3_Max_OSCr = 0, Dz_4_Max_OSCr = 0, Dz_5_Max_OSCr = 0;
	size_t num_Best_Node_Update = 0, num_Deleted_Candidate_in_Stack = 0, temp_num_Deleted_Candidate_in_Stack = 0;
	double OSC_Ratio = 0;
	bool DoubleDecoder;
	bool StackImformation;
	unsigned long long int CBC_Candidate = 0;
	unsigned long long int CBC_STE = 0;
	unsigned long long int DM_STE = 0;
	int CBC_length;
	//vector <size_t> Location_Record;
	//vector <__int8> Unsorted_Est;
	vector <__int8> CRC;
	//vector <__int8> CodeWordNoCRC;
	MATRIX<__int8> G_crc;   // ( k*(k+CRC length) )
	MATRIX<__int8> G_total; // G_crc*G (k*n)
	
	vector <__int8> message_seq,code_seq,code_seq_MRIP;
	vector<int> code_seq_permutation;
	double snr;
	double var;
	int Counter;
	unsigned long long int TotalCounter = 0;
	vector<int> Accumulate_Break;
	vector<double> Sorted_R;
	vector<double> Sorted_R_Base2;
	//vector<double> OSC_Ratio;
	//vector<double> OSC_Ratio_Max;
	//vector<double> OSC_Ratio_Min;
	int Frame_Counter;
		
	//Yin
	vector <__int8> estimated_codeword_66;
	vector <__int8> estimated_codeword;
	vector <double> rx_signal_seq;
	vector <int> initial_vector;
	//int FlipRecorder[20][20] = { 0 };
	size_t
		StackSize = 0,
		Constraint_i = 0,
		section1_i = 0,
		section2_i = 0,
		Constraint_j = 0,
		Control_Level = 0,
		DM_II_StackSize = 0,
		Check_Level = 0,
		CBC_FlippingBit = 0,
		Parity_Check_bits = 0,
		Multi_Stack_Size = 0,
		//
		MaxUsed_Stack = 0,
		Reached_Goal_Node = 0,
		Worst_Case_Candidate = 0,
		Number_of_the_last_symbols = 0,
		//
		SPA_I = 0,
		N_section;
		
		
	double
		Alpha = 0, // which is used in A*-Parity-f 
		OSC_Alpha = 0,
		Worst_Case_STE = 0,
		Worst_Case_COM = 0,
		//DM_STE = 0,
		DM_COM = 0,
		CandidateCodeWord = 0,
		STE = 0,
		phase = 0,
		STE_1 = 0,
		STE_2 = 0,
		STE_3 = 0,
		COM = 0,
		COM_1 = 0,
		COM_2 = 0,
		Binary_STE = 0,
		SNR = 0,
		Code_Rate = 0,
		First_nonzero_metric = 0,  // PoHan
		Constraint_i_ratio = 0,
		New_OSC_Alpha = 0,
		Del_Multi_Stack = 0,
		Fano_Metric_Parameter = 0,
		Worst_OSC_Ratio = 0,
		Worst_Metric_Ratio = 0;
	bool
		skip_Flag = FALSE;
	// Count the number of node which is deleted by DM
	size_t DeletedNode_counter;
	DECODING_INFO(size_t codeword_size, OPER_PARA oper_para)
	{
		StackSize = oper_para.StackSize;
		estimated_codeword.resize(codeword_size, 0);
		rx_signal_seq.resize(codeword_size, 0);
		initial_vector.resize(codeword_size, 0);
		{
			Alpha = oper_para.Alpha;
			OSC_Alpha = oper_para.OSC_Alpha;
			New_OSC_Alpha = oper_para.New_OSC_Alpha;
			Worst_Metric_Ratio = oper_para.Worst_metric_ratio;
			CBC_FlippingBit = oper_para.CBC_FlippingBit;
			Constraint_i = oper_para.Constraint_i;
			Constraint_j = oper_para.Constraint_j;
			Check_Level = oper_para.Check_Level;
			Control_Level = oper_para.Control_Level;
			Parity_Check_bits = oper_para.Parity_Check_bits;
			DM_II_StackSize = oper_para.DM_StackSize;
			Multi_Stack_Size = oper_para.Multi_Stack_Size;
			Fano_Metric_Parameter = oper_para.Fano_Metric_Parameter;
			Constraint_i_ratio = oper_para.Constraint_i_ratio;
			//Pin
			section1_i = oper_para.section1_i;
			section2_i = oper_para.section2_i;
			N_section = oper_para.N_section;
			SPA_I = oper_para.SPA_I;
			//Yin
			//CBC_length = Control_Level-
		}
		if (CRC_Check == TRUE) {
			if (CRC_length == 16) CRC = CRC_16;
			else if (CRC_length == 8) CRC = CRC_8;
			else if (CRC_length == 10) CRC = CRC_10;
			else if (CRC_length == 12) CRC = CRC_12;
			else if (CRC_length == 14) CRC = CRC_14;
			else if (CRC_length == 24) CRC = CRC_24;
			else cout << "Wrong CRC length!" << endl;
		}
	}

	void initial_para(){
		//DM_STE = 0;
		DM_COM = 0;
		DeletedNode_counter = 0;
		STE = 0;
		STE_1 = 0;
		STE_2 = 0;
		STE_3 = 0;
		COM = 0;
		COM_1 = 0;
		COM_2 = 0;
		New_OSC_Alpha = 0;
		Binary_STE = 0;
		CandidateCodeWord = 0;
		estimated_codeword.assign(initial_vector.begin(), initial_vector.end());
		rx_signal_seq.assign(initial_vector.begin(), initial_vector.end());
		num_Best_Node_Update = 0, num_Deleted_Candidate_in_Stack = 0, temp_num_Deleted_Candidate_in_Stack = 0;
	}

	void WorstSTE(size_t STE) {
		if (STE > Worst_Case_STE) Worst_Case_STE = STE;
	}

	//Yin
	double Combination_C(int i,int j){ // Mathematic Combination
		double ReturnNumber = 1;
		while (j > 0) {
			ReturnNumber = ReturnNumber * (i--) / (j--);
		}
		return ReturnNumber;
	}

	bool CRC_Examination(vector<__int8> Message_and_CRC, vector<__int8> variant) {
		int digit = 0;
		int Boundary = Message_and_CRC.size() - variant.size();
		vector<__int8> Dividend = Message_and_CRC;
		int Variant_Size = variant.size();
		//std::cout << Dividend.size() << endl;
		while (digit <= Boundary) {
			//std::cout <<endl<< digit << ":";
			if (Dividend.at(digit) == 1) {
				for (int i = 0; i < Variant_Size; ++i) {
					//std::cout << digit+i << ",";
					Dividend.at(i + digit) ^= variant.at(i);
				}
			}
			//for (int j = 0; j < Dividend.size(); ++j) std::cout << (int)Dividend.at(j) << " ";
			//std::cout << std::endl;
			++digit;
		}
		//std::cout << "AAAA";
		--digit;
		bool Result = TRUE;
		//std::cout << "o";
		while (digit < Message_and_CRC.size()) {
			if (Dividend.at(digit) == 1) {
				Result = FALSE;
				//std::cout << "x";
				break;
			}
			//for (int j = 0; j < Dividend.size(); ++j) std::cout << (int)Dividend.at(j) << " ";
			//std::cout << std::endl;
			++digit;
		}
		return Result;
	}

	void CRC_GMatrix_Genenerator(int Mes_Length, vector<__int8> CRC, MATRIX<__int8> &G);
	vector<__int8> Mes_To_MesCrc(vector<__int8> Message, MATRIX<__int8> G);
	MATRIX <__int8> G_times_G( MATRIX <__int8> G1, MATRIX <__int8> G2 );

};

class Adaptive_I_Decode_Info {

public:
	vector<__int8>
		Hard_RX,
		Hard_RX_Base1,
		Hard_RX_Base2,
		MRIP_codeword,
		MRIP_codeword_Base1,
		MRIP_codeword_Base2;
	MATRIX<__int8>
		Sorted_G,
		Sorted_G_Base1,
		Sorted_G_Base2;
	MATRIX<double>
		Metric_Table,
		Fano_Metric_Table,
		Metric_Table_Base1,
		Metric_Table_Base2;
	NODE_PATH Best_Goal;
	double OSC_metric_thr;

};

class NODE_COMB {
public:
	double metric = 0;
	vector<__int8> codeword_bits;
	NODE_COMB(int codeword_len) {
		codeword_bits.assign(codeword_len, 0);
	}
};
class MD_NODE {
public:
	vector<__int8> message_bits;
	int level = 0;
	size_t D_z = 0;
	MD_NODE(size_t len) {
		message_bits.assign(len, 0);
	}
};

#endif // !_CLASS_