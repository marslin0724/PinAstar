#pragma once
#define NULL 0

#include <math.h>
#include <iostream>
#include <vector>
#include <list>
#include<forward_list>
#include<deque>
#include <algorithm>
#include <random>
#include <functional>   // std::bit_xor
#include <numeric>
#include <time.h>


#include "MatrixOperation.h"
#include "DefineParameter.h"
#include "DefineClass.h"
#include "LinearBlockCodes.h"


using namespace std;

// Some opration functions for A*
/**********************************************************************/
void Generate_Message(size_t seq_length, vector<__int8> &message_seq);
void BPSK_Modulation(vector<__int8> &input_message_seq, vector<double> &output_signal_seq);


void Determine_Permutation( 
	vector<double>			&rx_signal, 
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&permutation_seq
);

void Determine_Permutation_MultiBase(
	vector<double>			&Rx_signal,
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&permutation_seq,
	int						Total_Basis,
	int						Current_Basis);

void Sort_Function( 
	vector<double>			&rx_signal,  
	vector<size_t>			&permutation_seq,  
	vector<double>			&sorted_rx_signal
);

void Sort_Function(
	vector<__int8>		&Rx_signal,
	vector<size_t>		&permutation_seq,
	vector<__int8>		&sorted_rx_signal
);

void Desort_Function( 
	vector<size_t>			&permutation_seq,  
	vector<__int8>			&estimate_sorted_codeword,
	vector<__int8>			&estimate_codeword
);

void Build_Metric_Table( 
	vector<double>			&sorted_rx_signal_seq, 
	MATRIX<double>			&metric_table
);

void Build_Fano_Metric_Table(
	vector<double>			&sorted_rx_signal,
	MATRIX<double>			&fano_metric_table,
	DECODING_INFO			&decoding_info
);			//PoHan

void Place_Node(
	vector<NODE_PATH>		&stack, 
	NODE_PATH				&child_node, 
	DECODING_INFO			&decoding_info
);

//Yin
void Place_Node_ver2(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info, vector<__int8> &Record);

void Place_Node_Fano(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info);

void Place_Node_multistack(
	vector<NODE_PATH>		&Stack,
	NODE_PATH				&child_node,
	DECODING_INFO &decoding_info
);

void Pre_Procedure(
	vector<double>			&rx_signal, 
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&Sorted_G,
	vector<size_t>			&permutation_seq, 
	MATRIX<double>			&metric_table
);
//Yin

void ParityCheck_Permutation(
	vector<double>			&ThirdParityBitLocation,
	MATRIX<__int16>			&ParityCheckbits,
	MATRIX<__int16>			&sorted_ParityCheckbits,
	vector<size_t>			&permutation_seq
);

void Sort_Parity_Length(vector<size_t>	&input_seq, vector<size_t>	&permutation_seq, vector<size_t> &output_seq);

void Pre_Procedure(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	DECODING_INFO       &decoding_info);

void Pre_Procedure(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	MATRIX<double>		&fano_metric_table,
	DECODING_INFO       &decoding_info);	//PoHan

void Pre_Procedure_Fano(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&fano_metric_table,
	DECODING_INFO       &decoding_info);

void Pre_Procedure_MultiBase(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	int					Total_Basis,
	int					Current_Basis,
	DECODING_INFO       &decoding_info);	//PoHan

void Pre_Procedure_MultiBase(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	MATRIX<double>		&fano_metric_table,
	int					Total_Basis,
	int					Current_Basis,
	DECODING_INFO       &decoding_info);

void Extend_Node_Procedure(
	NODE_PATH				&pointer, 
	NODE_PATH				&child_node, 
	MATRIX<double>			&metric_table, 
	__int8					&new_bit
);

void Extend_Node_Procedure_Fano(
	NODE_PATH		&Pointer,
	NODE_PATH		&child_node,
	MATRIX<double>	&metric_table,
	MATRIX<double>	&fano_metric_table,
	__int8			&new_bit
);

 bool Update_Best_Goal_Procedure(
	 NODE_PATH				&child_node, 
	 NODE_PATH				&best_goal, 
	 vector<NODE_PATH>		&stack
 );
 bool Update_Best_Goal_Procedure(
	 NODE_PATH				&Child_Node,
	 NODE_PATH				&Best_Goal,
	 vector<NODE_PATH>		&Stack,
	 DECODING_INFO			&decoding_info
 );
 bool Update_Best_Goal_Procedure(
	 NODE_PATH				&Child_Node,
	 NODE_PATH				&Best_Goal,
	 vector<NODE_PATH>		&Stack,
	 vector<NODE_PATH>		&Multi_Stack
 );

 void swap_PoHan(NODE_PATH &a, NODE_PATH &b);
 int Partition(vector<NODE_PATH> &Stack_CBC1, int front, int end, DECODING_INFO &decoding_info);
 void QuickSort(vector<NODE_PATH> &Stack_CBC1, int front, int end, DECODING_INFO &decoding_info);

// Original
/**********************************************************************/
// Hard Decision
void Hard_decision_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

/* A*I */
void A_star_I(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);


// Heuristic Function
/**********************************************************************/
/* A*Parity */
void A_star_Parity(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
__int8 Bit_Encoder(MATRIX<__int8> &G, size_t parity_index, NODE_PATH &node);
void Q_Function(MATRIX<__int8> &G, MATRIX<size_t>& Q_table);

/* A*Parity-f */
void Place_Node_f(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info, double alpha);
void A_star_Parity_f(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);

// Stack Reduction
/**********************************************************************/
void A_star_1_stack(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_3_stack(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_4_stack(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);

void A_star_PC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_PC_out(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);

void A_star_2_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_3_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_4_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);

void A_star_1_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_3_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);

void A_star_Parity_PC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_Parity_PC_out(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);

void A_star_Parity_2_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_Parity_3_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_Parity_4_stack(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);

void A_star_2_multiple_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// NEW Stack
/**********************************************************************/
void Place_Node_Qucik_Sorting(vector<NODE_PATH> &stack, NODE_PATH child_node, DECODING_INFO &decoding_info);
void QuickSort(vector<NODE_PATH> &arr, size_t front, size_t end, double &ComNN_temp);
size_t Partition(vector<NODE_PATH> &arr, size_t front, size_t end, double &ComNN_temp);
void swap(NODE_PATH &a, NODE_PATH &b);

void A_star_I_QuickSort(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_I_New_Comparsion(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// DM-I
/**********************************************************************/
void ParityPath_Checking(MATRIX<__int8> &G, DECODING_INFO &decoding_info, MATRIX<double> &Metric_Table, NODE_PATH &Child_Node, NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack, vector <__int8> &Hard_RX);
void A_star_PC_out_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

void A_star_1_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_4_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

void A_star_1_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_4_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// DM-II
/**********************************************************************/
void A_star_PC_DM_II(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_PCout_DM_II(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
bool Deleting_Mechanisim(MATRIX<__int8> &G, NODE_PATH &Check_Node, vector<__int8> &Hard_RX, DECODING_INFO &decoding_info);

void A_star_1_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_2_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_3_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);

void A_star_2_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_4_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

void A_star_BMA(MATRIX<__int8> &G, DECODING_INFO &decoding_info);


// OSC
/**********************************************************************/
void A_star_PC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_1_stack_OSC(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);
void A_star_PC_out_OSC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_2_stack_OSC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_3_stack_OSC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_4_stack_OSC(MATRIX<__int8> &G_matrix, DECODING_INFO &decoding_info);
void A_star_2_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_4_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_multiple_stack_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_stack_DM_I_Parity_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
double OSC_early_termination(MATRIX<double> &Metric_Table, size_t d_min, vector<__int8> &codeword, vector<__int8> &h_seq);
double OSC_early_termination2(
	MATRIX<double>		&Metric_Table,
	size_t				d_min,
	vector<__int8>		&codeword,
	vector<__int8>		&z_seq,
	NODE_PATH			&candidate,
	size_t				message_length);

void A_star_1_stack_PC(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info);


// 20190629 CBC
/**********************************************************************/
NODE_PATH Control_Band_Check_1bit(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info);

NODE_PATH Control_Band_Check_1bit(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table, 
	vector <__int8>				&Input_codeword_seq, 
	vector <__int8>				&Hard_RX, 
	NODE_PATH					&Input_Node, 
	NODE_PATH					&Best_Goal_Node, 
	double						OSC_metric_thr, 
	double						Worst_metric, 
	DECODING_INFO				&Decoding_info); 

NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info);

NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	double						OSC_metric_thr,
	DECODING_INFO				&Decoding_info);

NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
    int                          CBC1);

NODE_PATH Control_Band_Check_3bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
	double                       limited_ratio);

// Yin
NODE_PATH Control_Band_Check_2bits_Limited(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
	size_t limited_CBC_length);

void A_star_1_stack_CBC_OSC(
	MATRIX<__int8>			&G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Node,
	NODE_PATH				&Pre_Best_Goal,
	size_t					pc_i,
	double					metric_thr,
	DECODING_INFO			&decoding_info);

void A_star_2_stack_CBC_OSC(
	MATRIX<__int8>			&G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Node,
	NODE_PATH				&Pre_Best_Goal,
	size_t					pc_i,
	double					metric_thr,
	DECODING_INFO			&decoding_info);

void A_star_PC_out_CBC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_stack_CBC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);


void A_star_PC_out_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_stack_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_3_stack_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// Yin
void A_star_PC_out_CBC_OSC_Verified(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Verified_Function(
	MATRIX<__int8> &G, 
	DECODING_INFO &decoding_info,
	vector<double> &Rx,
	vector<__int8> &outoputcodeword,
	long double &metric,
	double &Average_LLR); // for A* in function

void A_star_PC_out_CBC_OSC_MultiDecoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Block(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info);   // 給上一行的function使用

void A_star_PC_out_CBC_OSC_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_Limited_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Dynamic_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

void A_star_PC_out_CBC_OSC_Test(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Test2(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

// PoHan
void A_star_PC_out_CBC_OSC_Adaptive_i_NEW(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_OSC_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_OSC_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_III(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_multiple_basis_CBC_OSC(
	MATRIX<__int8>			&Sorted_G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Best_Goal,
	vector<NODE_PATH>		&Stack_CBC,
	double					OSC_metric_thr,
	int						Basis,
	DECODING_INFO			&decoding_info);
void A_star_Fano_metric(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Adaptive_i_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Block_Fano_First(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1);
void A_star_PC_out_CBC_OSC_Block_Fano_Last(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1);
void A_star_2_Base_PC_out_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Latest(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_multiple_stack(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info, vector<NODE_PATH> &Multi_Stack);
void A_star_2_Base_PC_out_CBC_OSC_Parallel(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_Block_First(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1);
void A_star_PC_out_CBC_OSC_Block_New(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_Current, vector<NODE_PATH> &Stack_Next);
void A_star_PC_out_CBC_OSC_Block_Last(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1);
void A_star_RSS_Test(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient_Condition(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient_Condition(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info);
void A_star_PC_out_CBC_OSC_WorstMetric(MATRIX<__int8> &G, DECODING_INFO &decoding_info); 
void A_star_PC_out_Dynamic_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info);

//Pin
void A_star_Segment(MATRIX<__int8>& G, DECODING_INFO& decoding_info);
void A_star_Segment_ver2(MATRIX<__int8>& G, DECODING_INFO& decoding_info);
inline void Pre_Procedure_Segment(DECODING_INFO& decoding_info,
	MATRIX<double>& Metric_Table,
	vector<NODE_PATH>& S_Stack, vector<NODE_PATH>& C_Stack, size_t segment_length, size_t offset);
inline void Combine_Segment(NODE_PATH& Update_node, vector<NODE_PATH>& C_Stack_Update, vector<NODE_PATH>& C_Stack, vector<NODE_PATH>& S_Stack1, vector<NODE_PATH>& S_Stack2,
	NODE_PATH& Best_Goal, size_t segment_length,
	DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table);
inline void Combine_Segment(NODE_PATH &Update_node, vector<NODE_COMB>& C_Stack_Update, vector<NODE_COMB> &C_Stack, vector<NODE_PATH>& S_Stack1, vector<NODE_PATH>& S_Stack2,
	NODE_PATH &Best_Goal, size_t segment_length,
	DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table);
inline void Combine_Segment_ver2(NODE_PATH &Update_node, vector<NODE_PATH> &C_Stack, forward_list<NODE_PATH>& O_Stack,
	NODE_PATH &Best_Goal, size_t segment_length,
	DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table);

inline void Place_C_Stack(vector<NODE_PATH>& Stack, NODE_PATH& child_node, DECODING_INFO& decoding_info);
inline void Place_C_Stack(vector<NODE_COMB>& Stack, NODE_COMB& child_node, DECODING_INFO& decoding_info);
inline void Update_Stack( NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack);
inline void Update_Stack(NODE_PATH &Best_Goal, deque<NODE_PATH> &Stack);
inline void Update_Stack(NODE_PATH &Best_Goal, vector<NODE_COMB> &Stack);
inline bool Update_Best_Goal_Segments(NODE_PATH &Child_Node, NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack);
inline void Place_O_Stack(forward_list<NODE_PATH>& Stack, NODE_PATH& child_node, DECODING_INFO& decoding_info);
inline void next_extend(NODE_PATH& Node, vector<NODE_PATH> &C_Stack, forward_list<NODE_PATH>& O_Stack, NODE_PATH &Best_Goal,
	int segment_length, int message_length, DECODING_INFO& decoding_info);
inline bool meg_equal(vector<__int8> meg1, vector<__int8> meg2, int len);
inline void test(NODE_PATH& node, MATRIX<double>& Metric_Table);
//Pin test

void Hard_decision_test(MATRIX<__int8> &G, DECODING_INFO &decoding_info);


