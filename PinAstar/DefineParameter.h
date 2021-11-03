#pragma once
#ifndef _PARAMETER_
#define _PARAMETER_

#define PI (double)acos(-1.)
#define	Es 1.
#define AWGN_Mean 0.
#define FALSE	0
#define FLASE   0
#define TURE    1
#define TRUE	1

#define display_channel_info_flag	FALSE
#define decodeing_detail_flag		FALSE

// Iterative Decoding
#define Diff_Candidate_Test         FALSE
// TRUE: SumProduct, FALSE:MinSum
#define SumProduct                  TRUE
#define Modified_LLR_v2c            FALSE
#define Modified_LLR_v2c_factor     1.0
#define Modified_LLR_c2v            FALSE
#define Modified_LLR_c2v_factor     1.0
#define EarlyTermination            FALSE
#define Candidate_Check_Period      1
#define Interleaver_On              FALSE

//LDPC Code
#define SPA_Iter                    50
#define LDPC_H_Check_Equation       TRUE
#define Rx_diversity                FALSE

// Turbo Code
#define Iteration_Number            120
#define Local_Iteration             1
#define Accumulated_LLR_Calculate   FALSE
#define Acc_LLR_Scaling_Factor      1.0

#define Find_Reliable_Index         FALSE
#define Combine_Two_CodeWords       FALSE
#define CTC_factor                  1.0

#define Astar_after_TurboDecoding   FALSE
#define Syndrome_Check              FALSE
#define EXIT_Chart_Record           FALSE
#define LLR_MAX                     30
#define Message_Equals_0            FALSE
#define Bit_Flipping                FALSE
#define Test_parameter              1.0
#define ML_Likely_Method            FALSE
#define Repetition_Method           FALSE
#define Renew_One_Bit               FALSE

//Polar Code
#define BEC_Error                   0.5
#define BP_Decode_Iter              50
#define Polar_Astar_Outer_K         128
#define Polar_Astar_Outer_N         192   // = Polar_Astar_Inner_K
#define Polar_Astar_Inner_N         256
#define Polar_with_Astar            FALSE

#define StartPeriod                 1000
#define PeriodLength                1000

// Hybrid Polar + Astar 
#define Outer_Astar_Length          128

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))
//#define Eb_N0_1_or_Es_N0_2          2       // 1: Eb/N0
//                                            // 2: Es/N0

//Yin
//#define Adaptive_i_Method           FALSE
#define Adaptive_i_Parameter        1.2
#define Adaptive_i_Parameter_2		1.3
#define Adaptive_i_Decoder1_i       3       // Constraint-i of First Decoder 
#define Adaptive_i_Decoder2_i       5       // Constraint-i of Second Decoder 
#define Adaptive_i_Decoder3_i       0       // Constraint-i of Third Decoder 
#define Adaptive_i_Decoder4_i       0       // Constraint-i of Fourth Decoder 

#define Dynamic_I_Constraint        FALSE
#define Error_Floor_constant        8
#define CRC_Check                   FALSE   // Use CRC to increase the speed of tree searching
#define CRC_length                  14       // 8,16
#define CRC_8                       { 1,0,0,1,1,0,0,0,1 }          //  Original: { 1,0,0,1,1,0,0,0,1 } , bluetooth: { 1,1,0,1,0,0,1,1,1 }
#define CRC_10                      {1,1,0,0,0,1,1,0,0,1,1}
#define CRC_12                      {1,1,0,0,0,0,0,0,0,1,1,1,1}
#define CRC_14                      {1,0,0,0,0,1,1,0,1,1,1,0,1,1,1}
#define CRC_16                      {1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1}
#define CRC_24                      {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1} // WCDMA

// PoHan
//define Multi_Bases
#define Sorted_Base_1				1
#define Sorted_Base_2				2
#define Sorted_Base_3				3
#define Multiple_Basis_Bits			16
// define Sufficient Condition	
#define dmin						34


// RM Code (16,11)  (r,m) = (2,4)
#define RM_16_v0                       {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_16_v4                       {0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1}
#define RM_16_v3                       {0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1}
#define RM_16_v2                       {0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1}
#define RM_16_v1                       {0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1}
#define RM_16_v3v4                     {0,0,0,0,0,0,0,0, 0,0,0,0,1,1,1,1}
#define RM_16_v2v4                     {0,0,0,0,0,0,0,0, 0,0,1,1,0,0,1,1}
#define RM_16_v1v4                     {0,0,0,0,0,0,0,0, 0,1,0,1,0,1,0,1}
#define RM_16_v2v3                     {0,0,0,0,0,0,1,1, 0,0,0,0,0,0,1,1}
#define RM_16_v1v3                     {0,0,0,0,0,1,0,1, 0,0,0,0,0,1,0,1}
#define RM_16_v1v2                     {0,0,0,1,0,0,0,1, 0,0,0,1,0,0,0,1}

// RM (32,6)    (r,m)=(1,5)
#define RM_32_v0                       {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_32_v5                       {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_32_v4                       {0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1}
#define RM_32_v3                       {0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1}
#define RM_32_v2                       {0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1}
#define RM_32_v1                       {0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1}

// RM (64,42)    (r,m)=(3,6)
#define RM_64_v0                       {1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_64_v6                       {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_64_v5                       {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1}
#define RM_64_v4                       {0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0,0, 1,1,1,1,1,1,1,1}
#define RM_64_v3                       {0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1, 0,0,0,0,1,1,1,1}
#define RM_64_v2                       {0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1, 0,0,1,1,0,0,1,1}
#define RM_64_v1                       {0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1}



// Basic Versions
#define Hard_Decision_Decoder           0
#define A_Star_Original					1	// Áú¥Ã²»
#define A_Star_Han						2	// Áú¥Ã²»
#define A_Star_Parity					3	// ±i®a»²
#define A_Star_Parity_f					4	// ±i®a»² 	

// Stack Reduction
#define A_Star_PC						5	// ±i®a»²
#define A_Star_PC_out					6	// µ¾µÏ¡B§Ó¤ß
#define A_Star_2_Stack					7	// µ¾µÏ¡B§Ó¤ß
#define A_Star_3_Stack					8	// µ¾µÏ¡B§Ó¤ß
#define A_Star_4_Stack					9	// µ¾µÏ¡B§Ó¤ß

// Stack Reduction + Parity
#define A_Star_Parity_PC				10	
#define A_Star_Parity_PC_out			11	
#define A_Star_2_Parity_Stack			12
#define A_Star_3_Parity_Stack			13
#define A_Star_4_Parity_Stack			14

// DM-I
#define A_Star_PC_out_DM_I				15
#define A_Star_2_Stack_DM_I				16
#define A_Star_3_Stack_DM_I				17
#define A_Star_4_Stack_DM_I				18	

// DM-I + Parity
#define A_Star_PC_out_DM_I_Parity		19
#define A_Star_2_Stack_DM_I_Parity		20
#define A_Star_3_Stack_DM_I_Parity		21
#define A_Star_4_Stack_DM_I_Parity		22

// DM-II
#define A_Star_PC_DM_II					23
#define A_Star_PCout_DM_II				24
#define A_Star_2_Stack_DM_II			25
#define A_Star_3_Stack_DM_II			26
#define A_Star_4_Stack_DM_II			27
#define A_Star_BMA						28

// OSC
#define A_Star_PC_OSC					29
#define A_Star_PCout_OSC				30
#define A_Star_2_Stack_OSC				31
#define A_Star_3_Stack_OSC				32
#define A_Star_4_Stack_OSC				33

// CBC
#define A_Star_PCout_CBC				34
#define A_Star_2_Stack_CBC				35
#define A_Star_3_Stack_CBC				36

// CBC + OSC
#define A_Star_PCout_CBC_OSC			37
#define A_Star_2_Stack_CBC_OSC			38
#define A_Star_3_Stack_CBC_OSC			39
#define A_Star_1_Stack_CBC_OSC			40 // Chou Yin Debug not finished

// Yin
#define A_Star_PC_out_Limited_CBC_OSC        41
#define A_Star_PC_out_CBC_OSC_Dynamic_I      42
#define A_Star_PC_out_CBC_OSC_MultiDecoder   43
#define A_Star_PC_out_CBC_OSC_Verified       44

#define A_Star_PC_out_CBC_OSC_Adaptive_i     66

//PoHan
#define A_Star_PC_OSC_II													45
#define A_Star_PC_out_Dynamic_CBC_OSC										46
#define A_Star_PC_out_CBC_OSC_II											47
#define A_Star_PC_out_CBC_OSC_III											48
#define A_Star_RSS_test														49
#define A_Star_PC_out_CBC_OSC_WorstMetric									50
#define A_Star_2_Multi_Stack												55
#define A_Star_2_Multi_Stack_OSC											56
#define A_Star_PC_out_CBC_OSC_Adaptive_i_New								67
#define A_Star_Fano															70
#define A_Star_PC_out_CBC_OSC_Fano											71
#define A_Star_PC_out_CBC_OSC_Adaptive_i_Fano								72
#define A_Star_2_Base_PC_out_CBC_OSC										101
#define A_Star_2_Base_PC_out_CBC_OSC_Latest									102
#define A_Star_2_Base_PC_out_CBC_OSC_Parallel								103
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i								104
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel					105
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient				108
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient	109
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano						106
#define A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano				107

//Pin
#define A_Star_Segment_Orignal			110
#define A_Star_Segment_Orignal_ver2		111
#define Hard_test						112
#define MinD_test						113
#define A_Star_Section_PC				114
#define A_Star_Section_PC_out			115
#define SPA_A_Star						116
#define SPA_A_Star_ver2					117
#define A_Star_N_Section				118
#define A_Star_5_Stack					119
#define A_Star_6_Stack					120
#define A_Star_7_Stack					121
#define SPA_AvgLLR						122
#define A_Star_Rep						123



// Other Decoding Algorithms
#define Polar_Code_BP_Decoder                80
#define LDPC_Sum_Product_Decoder             81
#define Majority_Soft_Decoding               82
#define Majority_Hard_Decoding               83
#define Majority_Soft_Decoding_ver2          84
//#define Turbo_LDPC_Decoder                   85
#define LDPC_VC_RBP_Algorithm                85
#define LDPC_Modified_RBP_Algorithm          86

// Test
#define A_Star_PC_out_CBC_OSC_Test           99
#define A_Star_PC_out_CBC_OSC_Test2         100

#endif //!_PARAMETER.H_