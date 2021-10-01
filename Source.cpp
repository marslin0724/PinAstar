#include "PinAstar/Source.h"
// 常犯錯誤: Eb/No的CodeRate !!!!!

// FALSE -> transmit zero message 
// TRUE  -> Randomly generate transmitted message

#define RandomMessage_flag TRUE

// some ideas with which teacher comes up
// but the effect is not obvious.
#define Number_of_Last_Symbols 10

// Global Class
// Select a linear block code for simulation.
CODE LinearBlockCode;								

// Key in decoding parameters for simulation.
OPER_PARA Operation_Parameter(LinearBlockCode);

// Store the configuration of decoder and decoding statistic
DECODING_INFO Decoding_info(LinearBlockCode.Col_Number, Operation_Parameter);

using namespace std;

// Global Parameters
double start_time(clock());

const size_t
	Message_Length(LinearBlockCode.Row_Number),
	Codeword_Length(LinearBlockCode.Col_Number);

long double
	// Block_Counter(0),
	Temp_ErrorBit_Counter(0),
	ErrorBit_Counter(0),
	ErrorBit_Counter_1(0),
	ErrorBlock_Counter(0),
	ML_LB_ErrorBit_Counter,
	ML_LB_ErrorBlock_Counter(0);

double
BER(0),				// for error rate
BER_1(0),
BLER(0),
ML_LB_BER(0),
ML_LB_BLER(0),
// A*,for the number of searching tree edges
Total_Avg_STE(0),	// Total
Avg_STE(0),			// Average per message bit
Total_Avg_STE_1(0),
Avg_STE_1(0),
Total_Avg_STE_2(0),
Avg_STE_2(0),
Total_Avg_STE_3(0),
Avg_STE_3(0),
// A* ,for the number of conparisons 
Total_Avg_Com(0),	// Total
Avg_Com(0),			// Average per message bit
Total_Avg_Com_1(0),
Avg_Com_1(0),
Total_Avg_Com_2(0),
Avg_Com_2(0),


Total_Avg_CandidateCodeWord(0),
Avg_CandidateCodeWord(0),

Total_Avg_Alpha(0),	//PoHan
Avg_Binary_STE(0),
Total_Avg_Binary_STE(0),
Total_Updated_Best_Node(0),
Total_Deleted_Candidate(0),
Avg_Deleted_Candidate(0),

Avg_Alpha(0),
// DM,for the number of searching tree edges
Total_Avg_DM_STE(0),// Total
Avg_DM_STE(0),		// Average per message bit
					// DM ,for the number of conparisons 
Total_Avg_DM_Com(0),// Total
Avg_DM_Com(0),		// Average per message bit
						// DM ,for the number of deleting nodes
Total_Avg_DN(0),	// Total
Avg_DN(0),			// Average per decoding
Transmitted_Block(0);

	
vector<__int8>
	Message_Seq(Message_Length, 0),
	Codeword_Seq(Codeword_Length, 0),
 //	Error_Seq(Codeword_Length, 0);
    Error_Seq(Message_Length, 0);

vector<double> Tx_signal_seq(Codeword_Length, 1);
// File Name
string
Realtime_Data = "RealTime_Decoding_Statistic.txt",
Complete_Data = "Complete_Decoding_Statistic.txt",
OSC_Data = "OSC_Statistic.txt",
Best_Goal_OSC_Data = "Best_Goal_OSC_Statistic.txt",
RSS_Data = "RSS_Data.txt";



// LDPC part
//LDPC_FUNC LDPC;
//string H_file = "H_1024_768_z64_196.txt";
//MATRIX<__int8> LDPC_H;
//MATRIX<__int8> LDPC_G;
//vector <size_t> LDPC_Permutation_Seq;
//#define LDPC_CodeRate   0.75
//#define System_CodeRate 0.53   // Eb/No (LDPC: 0.58 / Product: 0.53 ) Es/No 1
//#define NumberOfAstar   6
//#define SPA_Iter        10

void Sort_function(
	vector<__int8>		&Rx_signal,
	vector<size_t>		&permutation_seq,
	vector<__int8>		&sorted_rx_signal);

void main()
{
	//cout << sqrt(13) << endl;
	/*
	vector<int> Test = { 1,2,3,4,5,6,7,8,9,10 };
	int size = 10;
	for (int i = 0; i < size; ++i) {
		Test.insert(Test.begin() + 2*i ,0);
	}
	//Test.push_back(0);
	for (int i = 0; i < Test.size(); ++i) {
		cout << Test.at(i) << " ";
	}
	system("pause");
	*/

	ClearFile(Complete_Data);
	//Pin test
	//ClearFile("Hard_test.txt");
	WriteFile_Parameters(Complete_Data, LinearBlockCode, Operation_Parameter);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0, 1.0);
	srand(NULL);
	//srand((unsigned)time(NULL));

	//printf(" ** Flip bits starting from the last. \n");

	if (RandomMessage_flag)
		printf(" ** Randomly generate message sequence in equal likely. \n\n");
	else
		printf(" ** Message is all zero sequence. \n\n");

	/**************************************************************/

	// Pointer of function
	// used to choose the decoding algorithm we want to simulate. 
	
	void(*Decoder)(MATRIX<__int8> &, DECODING_INFO &) = NULL;
	
	switch (Operation_Parameter.Decoder_Version) {
		
	//printf("\n\nHI\n\n");
	//system("pause");
	// Basic A*
	case Hard_Decision_Decoder:
		Decoder = Hard_decision_Decoder;
		break;
	case A_Star_Original:
		Decoder = A_star_I;
		
		// Other methods of basic A* arranging the stored nodes
		// --(1)-- Arrange the elements at stack by Quick Sorting 
		// Decoder = A_star_I_Decoder_QuickSort;		
		// --(2)-- store the position of previous inserted node, 
		// reduceing the number of comparison
		// Decoder = A_star_I_Decoder_New_Comparsion;	
		break;
	case A_Star_Parity:
		Decoder = A_star_Parity;
		break;
	case A_Star_Parity_f:
		Decoder = A_star_Parity_f;
		break;
	
	// Stack Reductions
	case A_Star_PC:
		Decoder = A_star_PC;
		break;
	case A_Star_PC_out:
		Decoder = A_star_PC_out;
		break;
	case A_Star_2_Stack:
		Decoder = A_star_2_stack;
		break;
	case A_Star_3_Stack:
		Decoder = A_star_3_stack;
		break;
	case A_Star_4_Stack:
		Decoder = A_star_4_stack;
		break;
	
	// Stack Reduction + Parity 
	case A_Star_Parity_PC:
		Decoder = A_star_Parity_PC;
		break;
	case A_Star_Parity_PC_out:
		Decoder = A_star_Parity_PC_out;
		break;
	case A_Star_2_Parity_Stack:
		Decoder = A_star_Parity_2_stack;
		break;
	case A_Star_3_Parity_Stack:
		Decoder = A_star_Parity_3_stack;
		break;
	case A_Star_4_Parity_Stack:
		Decoder = A_star_Parity_4_stack;
		break;
	
	// Deleting Mechanism I
	case  A_Star_PC_out_DM_I:
		Decoder = A_star_PC_out_DM_I;
		break;
	case A_Star_2_Stack_DM_I:
		Decoder = A_star_2_stack_DM_I;
		break;
	case A_Star_3_Stack_DM_I:
		Decoder = A_star_3_stack_DM_I;
		break;
	case A_Star_4_Stack_DM_I:
		Decoder = A_star_4_stack_DM_I;
		break;
	

	// Deleting Mechanism I + Parity 
	case A_Star_PC_out_DM_I_Parity:
		Decoder = A_star_PC_out_DM_I_Parity;
		break;
	case A_Star_2_Stack_DM_I_Parity:
		Decoder = A_star_2_stack_DM_I_Parity;
		break;
	case A_Star_3_Stack_DM_I_Parity:
		Decoder = A_star_3_stack_DM_I_Parity;
		break;
	case A_Star_4_Stack_DM_I_Parity:
		Decoder = A_star_4_stack_DM_I_Parity;
		break;
	

	//	Deleting Mechanism II 
	case A_Star_PC_DM_II:
		Decoder = A_star_PC_DM_II;
		break;
	case A_Star_PCout_DM_II:
		Decoder = A_star_PCout_DM_II;
		break;
	case A_Star_2_Stack_DM_II:
		Decoder = A_star_2_Stack_DM_II;
		break;
	case A_Star_3_Stack_DM_II:
		Decoder = A_star_3_Stack_DM_II;
		break;
	case A_Star_4_Stack_DM_II:
		Decoder = A_star_4_Stack_DM_II;
		break;
	case A_Star_BMA:
		Decoder = A_star_BMA;
		break;
	
	
	// OSC
	case A_Star_PC_OSC:
		Decoder = A_star_PC_OSC;
		break;
	case A_Star_PCout_OSC:
		Decoder = A_star_PC_out_OSC;
		break;
	case A_Star_2_Stack_OSC:
		Decoder = A_star_2_stack_OSC;
		break;
	case A_Star_3_Stack_OSC:
		Decoder = A_star_3_stack_OSC;
		break;
	case A_Star_4_Stack_OSC:
		Decoder = A_star_4_stack_OSC;
		break;


	// CBC
	case A_Star_PCout_CBC:
		Decoder = A_star_PC_out_CBC;
		break;
	case A_Star_2_Stack_CBC:
		Decoder = A_star_2_stack_CBC;
		break;
	

	// CBC + OSC
	case A_Star_PCout_CBC_OSC:
		Decoder = A_star_PC_out_CBC_OSC;
		break;
	case A_Star_2_Stack_CBC_OSC:
		Decoder = A_star_2_stack_CBC_OSC;
		break;
	case A_Star_3_Stack_CBC_OSC:
		Decoder = A_star_3_stack_CBC_OSC;
		break;
	case A_Star_1_Stack_CBC_OSC:
		Decoder = A_star_PC_out_CBC_OSC;
		break;

		// Chou-Yin
	case A_Star_PC_out_Limited_CBC_OSC:
		Decoder = A_star_PC_out_Limited_CBC_OSC;
		break;
	case A_Star_PC_out_CBC_OSC_Dynamic_I:
		Decoder = A_star_PC_out_CBC_OSC_Dynamic_I;
		break;
	case A_Star_PC_out_CBC_OSC_MultiDecoder:
		Decoder = A_star_PC_out_CBC_OSC_MultiDecoder;
		break;

	case A_Star_PC_out_CBC_OSC_Verified:  // 最新版的 "A_Star_PCout_CBC_OSC"(37)
		Decoder = A_star_PC_out_CBC_OSC_Verified;
		break;
	case A_Star_PC_out_CBC_OSC_Adaptive_i:
		Decoder = A_star_PC_out_CBC_OSC_Adaptive_i;
		break;

		// PoHan
	case A_Star_PC_OSC_II:
		Decoder = A_star_PC_OSC_II;
		break;
	case A_Star_PC_out_Dynamic_CBC_OSC:
		Decoder = A_star_PC_out_Dynamic_CBC_OSC;
		break;
	case A_Star_PC_out_CBC_OSC_II:
		Decoder = A_star_PC_out_CBC_OSC_II;
		break;
	case A_Star_PC_out_CBC_OSC_III:
		Decoder = A_star_PC_out_CBC_OSC_III;
		break;
	case A_Star_RSS_test:
		Decoder = A_star_RSS_Test;
		break;
	case A_Star_PC_out_CBC_OSC_WorstMetric:
		Decoder = A_star_PC_out_CBC_OSC_WorstMetric;
		break;
	case A_Star_2_Multi_Stack:
		Decoder = A_star_2_multiple_stack;
		break;
	case A_Star_2_Multi_Stack_OSC:
		Decoder = A_star_2_multiple_stack_OSC;
		break;
	case A_Star_PC_out_CBC_OSC_Adaptive_i_New:
		Decoder = A_star_PC_out_CBC_OSC_Adaptive_i_NEW;
		break;
	case A_Star_Fano:
		Decoder = A_star_Fano_metric;
		break;
	case A_Star_PC_out_CBC_OSC_Fano:
		Decoder = A_star_PC_out_CBC_OSC_Fano;
		break;
	case A_Star_PC_out_CBC_OSC_Adaptive_i_Fano:
		Decoder = A_star_PC_out_CBC_OSC_Adaptive_i_Fano;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC:
		Decoder = A_star_2_Base_PC_out_CBC_OSC;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Latest:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Latest;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Parallel:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Parallel;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient_Condition;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient_Condition;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano:
		Decoder = A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano;
		break;
		//
	case Majority_Soft_Decoding:
		Decoder = Majority_Soft_decoding;
		break;
	case Majority_Hard_Decoding:
		Decoder = Majority_Hard_decoding;
		break;
	case Majority_Soft_Decoding_ver2:
		//Decoder = Majority_Soft_decoding_ver2;
		break;

		// Other Decoding Algorithm
	case Polar_Code_BP_Decoder:
		Decoder = Polar_code_BP_Decoder;
		break;
	case LDPC_Sum_Product_Decoder:
		Decoder = SPA;
		break;
	case LDPC_VC_RBP_Algorithm:
		Decoder = LDPC_VC_RBP;
		break;
	case LDPC_Modified_RBP_Algorithm:
		Decoder = LDPC_Modified_RBP;
		break;

		//Pin
	case A_Star_Segment_Orignal:
		Decoder = A_star_Segment;
		break;
	case A_Star_Segment_Orignal_ver2:
		Decoder = A_star_Segment_ver2;
		break;
	case Hard_test:
		Decoder = Hard_decision_test;
		break;
	case MinD_test:
		Decoder = MD_test;
		break;
	case A_Star_Section_PC:
		Decoder = A_star_section_PC;
		break;
	case A_Star_Section_PC_out:
		Decoder = A_star_section_PC_out;
		break;

		//Test
	case A_Star_PC_out_CBC_OSC_Test:
		Decoder = A_star_PC_out_CBC_OSC_Test;
		break;
	case A_Star_PC_out_CBC_OSC_Test2:
		Decoder = A_star_PC_out_CBC_OSC_Test2;
		break;

	default:
		cout
			<< "\n Something happens error! "
			<< "\n Please restart the simulation procedure.";
		break;
	}
	

	// ************************************
	// *                                  *
	// *        Start Simulation          *
	// *                                  *
	// ************************************

	for (double snr_dB(Operation_Parameter.SNR_dB_start);
		snr_dB <= Operation_Parameter.SNR_dB_end;
		snr_dB += Operation_Parameter.SNR_dB_step) {

		double start_time_for_SNR(clock());
		
		// Initialize Decoding Statistical Parameters
		{
			//
			Total_Avg_Com = 0;
			Avg_Com = 0;
			Total_Avg_Com_1 = 0;
			Avg_Com_1 = 0;
			Total_Avg_Com_2 = 0;
			Avg_Com_2 = 0;
			
			//
			Total_Avg_STE = 0;
			Avg_STE = 0;
			Total_Avg_STE_1 = 0;
			Avg_STE_1 = 0;
			Total_Avg_STE_2 = 0;
			Avg_STE_2 = 0;
			Total_Avg_STE_3 = 0;
			Avg_STE_3 = 0;
			
			//
			Total_Avg_Alpha = 0;	//PoHan
			Avg_Alpha = 0;			//PoHan
			Total_Avg_Binary_STE = 0;
			Avg_Binary_STE = 0;
			//
			ErrorBit_Counter = 0;
			ErrorBit_Counter_1 = 0;
			ErrorBlock_Counter = 0;
			ML_LB_ErrorBit_Counter = 0;
			ML_LB_ErrorBlock_Counter = 0;
			BER = 0;
			BLER = 0;
			ML_LB_BER = 0;
			ML_LB_BLER = 0;
			//
			//Block_Counter = 0;
			//
			Decoding_info.MaxUsed_Stack = 0;
			Decoding_info.Worst_Case_STE = 0;
			Decoding_info.Worst_Case_COM = 0;
			Decoding_info.DeletedNode_counter = 0;
			Decoding_info.SNR = snr_dB;
			//
			Total_Avg_DM_STE = 0;
			Avg_DM_STE = 0;
			//
			Total_Avg_Com = 0;
			Avg_Com = 0;
			//
			Total_Avg_DN = 0;
			Avg_DN = 0;
			//
			Total_Updated_Best_Node = 0;
			Total_Deleted_Candidate = 0;
			Avg_Deleted_Candidate = 0;
			//
			Decoding_info.Worst_OSC_Ratio = 0;
			//
			Decoding_info.Dz_0_Number = 0;
			Decoding_info.Dz_1_Number = 0;
			Decoding_info.Dz_2_Number = 0;
			Decoding_info.Dz_3_Number = 0;
			Decoding_info.Dz_4_Number = 0;
			Decoding_info.Dz_5_Number = 0;
			Decoding_info.Dz_6_Number = 0;
			Decoding_info.Dz_7_Number = 0;
			Decoding_info.Dz_8_Number = 0;
			Decoding_info.Dz_9_Number = 0;
			Decoding_info.Dz_10_Number = 0;
			Decoding_info.Dz_11_Number = 0;
			Decoding_info.Dz_12_Number = 0;
			Decoding_info.Dz_13_Number = 0;
			Decoding_info.Dz_14_Number = 0;
			Decoding_info.Dz_15_Number = 0;
			Decoding_info.Dz_16_Number = 0;
			Decoding_info.Dz_17_Number = 0;
			Decoding_info.Dz_18_Number = 0;
			Decoding_info.Dz_19_Number = 0;
			Decoding_info.Dz_20_Number = 0;
			Decoding_info.Dz_21_Number = 0;
			Decoding_info.Dz_22_Number = 0;
			Decoding_info.Dz_23_Number = 0;
			/*
			//Pin test
			for (int i = 0; i < 128; i++) {
				Decoding_info.err_count[i] = 0;
			}*/
		}

		Show_Current_Time();
		
		// Show the current (E_b/N_0) [dB] and (E_s/N_0) [dB]
		cout
			<< "\n SNR  [dB] = " << snr_dB << " [dB]"
			<< " → Es/N0 [dB] = " << snr_dB + 10 * log10(LinearBlockCode.Code_Rate) << " [dB] \n"
			<< endl;
		// simulation at specific SNR 
		cout << "SPA Iter: " << SPA_Iter << ", Polar BP Iter: " << BP_Decode_Iter << ", Turbo Iter: " << Iteration_Number << endl << endl;
		// if Es/N0 change the variance of noise
		if (Operation_Parameter.Eb_N0_1_or_Es_N0_2 == 2) LinearBlockCode.Code_Rate = 1;
		Decoding_info.Code_Number = LinearBlockCode.Code_number;



		// Calculate Error in different index
		Decoding_info.Error_Accumulation.assign(192, 0);


		//LinearBlockCode.Code_number = 88;
		//      ************************************              Some Concatenation               ********************************************* 
		if (LinearBlockCode.Code_number > 89) {

			Decoding_info.CBC_length = Decoding_info.Control_Level - Message_Length;
			Decoding_info.var = (Es / (2.*LinearBlockCode.Code_Rate)) * pow(10., (-(snr_dB / 10.)));

			LDPC_FUNC LDPC;
			bool This_is_LDPC = FALSE;
			if (LinearBlockCode.Code_number >= 90 && LinearBlockCode.Code_number < 100) {
				This_is_LDPC = TRUE;
				Decoding_info.Return_Est_1_or_LLR_0 = 1;
			}
			if (LinearBlockCode.Code_number == 102) {
				This_is_LDPC = TRUE;
				Decoding_info.Return_Est_1_or_LLR_0 = 0;
			}
			//if (Adaptive_i_Decoder1_i != 0)Decoding_info.Adaptive_i_Iteration = 2;
			//else Decoding_info.Adaptive_i_Iteration = 1;
			//if (Adaptive_i_Decoder2_i != 0)Decoding_info.Adaptive_i_Iteration = 3;

			if (CRC_Check) {
				Decoding_info.CRC_GMatrix_Genenerator(LinearBlockCode.G.Row_number - CRC_length, Decoding_info.CRC, Decoding_info.G_crc); // 產生G_crc
				//cout << Decoding_info.G_crc.Row_number << Decoding_info.G_crc.Col_number << endl;
				//cout << LinearBlockCode.G.Row_number << LinearBlockCode.G.Col_number << endl;
				//Decoding_info.G_total = Decoding_info.G_times_G(Decoding_info.G_crc, LinearBlockCode.G);
			}
			if (Operation_Parameter.Decoder_Version == Polar_Code_BP_Decoder || Polar_with_Astar_middle_index) {
				Decoding_info.frozen = LinearBlockCode.frozen;
				Decoding_info.non_frozen = LinearBlockCode.non_frozen;
				//vector<int> non_frozen;
			}
			if (Dynamic_I_Constraint) {
				double Snr;
				if (LinearBlockCode.Col_Number == 192) {
					if (snr_dB == 2.0) Snr = 0.000596;
					else if (snr_dB == 2.5) Snr = 6.65e-05;
					else if (snr_dB == 3.0) Snr = 3.33e-06;
					else if (snr_dB == 3.5) Snr = 1.52e-07;
					else Snr = 7.71e-09;

				}
				else if (LinearBlockCode.Col_Number == 128)
				{
					if (snr_dB == 2.0) Snr = 0.001798;
					else if (snr_dB == 2.5) Snr = 0.000268;
					else if (snr_dB == 3.0) Snr = 2.34e-05;
					else if (snr_dB == 3.5) Snr = 1.17e-06;
					else Snr = 5.01e-08;
				}
				else {
					cout << "Adapting Constraint I Error!" << endl;
					break;
				}
				Snr *= Error_Floor_constant;
			
			}
			Decoding_info.Ave_LLR = 0;
			long double AveLLR_Before = 0;

			if (LinearBlockCode.Inner_Col_Number != 0) {
				Tx_signal_seq.resize(LinearBlockCode.Inner_Col_Number);
				Decoding_info.frozen = LinearBlockCode.frozen;
				Decoding_info.non_frozen = LinearBlockCode.non_frozen;
				Decoding_info.Channel_Parameter = LinearBlockCode.Channel_Parameter;
			}
			if (LinearBlockCode.Short_Inner_Col_Number != 0)  Tx_signal_seq.resize(LinearBlockCode.Col_Number*LinearBlockCode.Short_Inner_Col_Number / LinearBlockCode.Short_Inner_Row_Number);

			vector<double> LLR_Distribution_0(2000, 0);
			vector<double> LLR_Distribution_1(2000, 0);
			if (EXIT_Chart_Record) {
				vector<double> LLR_all_zero(150000, 0);
				Decoding_info.LLR_1_priori = LLR_all_zero, Decoding_info.LLR_2_priori = LLR_all_zero,
					Decoding_info.LLR_1_extrinsic = LLR_all_zero, Decoding_info.LLR_2_extrinsic = LLR_all_zero;
				Decoding_info.Turbo_Index = 0;
			}
			if (Diff_Candidate_Test) {
				Decoding_info.Iterative_Decoding_Candidates.Building_Empty_Matrix(Iteration_Number,LinearBlockCode.H.Col_number);
				Decoding_info.Iter_Decode_Candidates_HardDecision.Building_Empty_Matrix(Iteration_Number, LinearBlockCode.H.Col_number);
				Decoding_info.Iter_Decode_Candidates_Est.Building_Empty_Matrix(Iteration_Number, LinearBlockCode.H.Col_number);
			}
			Decoding_info.Amplitude_Flipping_number.Building_Empty_Matrix(12, 14);

			if (Combine_Two_CodeWords || ML_Likely_Method) {
				LinearBlockCode.G_Inner = LinearBlockCode.G_;
				for (int row = 0; row < LinearBlockCode.G_.Row_number; row++) {
					LinearBlockCode.G_Inner._matrix.at(row).insert(
						LinearBlockCode.G_Inner._matrix.at(row).end(),
						LinearBlockCode.G._matrix.at(row).begin()+ LinearBlockCode.G_.Row_number,
						LinearBlockCode.G._matrix.at(row).end());
				}
				LinearBlockCode.G_Inner.Col_number = 256;

			}

			double r0, r1, r2;
			r0 = 0;
			r1 = 0;
			r2 = 0;
			double Inner_Error = 0;
			long double LLR_u0 = 0;
			cout << endl << "CodeRate: " << LinearBlockCode.Code_Rate << endl;

			for (
				Transmitted_Block = 1;
				Transmitted_Block <= Operation_Parameter.block_number;
				++Transmitted_Block) {
				//printf("#%.0f: \n", Transmitted_Block);
				//cout << "GG";
				Decoding_info.initial_para();
				//
				//++Block_Counter;

				// Tx , Random message or all zero message 
				// RandomMessage_flag = 0 -> Random message 
				//					  = 1 -> All zero message

				if (RandomMessage_flag) {
					Generate_Message(Message_Length, Message_Seq);
					
					if (Message_Equals_0) {
						for (int i = 0; i < Message_Seq.size(); ++i) {
							Message_Seq.at(i) = 0;
						}
						//cout << "J";
					}
					if (Operation_Parameter.Decoder_Version == Polar_Code_BP_Decoder || (LinearBlockCode.Code_number >= 90 && LinearBlockCode.Code_number < 100)) {
						for (int i = 0; i < Codeword_Length; ++i) {
							int temp = 0;
							for (int j = 0; j < Message_Length; j++) {
								temp ^= (Message_Seq.at(j)&LinearBlockCode.G._matrix[j][i]);
								//temp += Message_Seq.at(j)*LinearBlockCode.G._matrix[j][i];
							}
							Codeword_Seq.at(i) = temp;
						}
						//cout << "J";
					}
					else if (LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC 
						|| LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_ver2
						|| LinearBlockCode.Code_number == Turbo_320_128_with_2LDPC
						|| LinearBlockCode.Code_number == Turbo_256_128_RC_LDPC_192_128)  // Turbo encoder
					{
						//cout << "A";
						Turbo_LDPC_encoder(LinearBlockCode.G_, LinearBlockCode.G, Message_Seq, Codeword_Seq);
						Tx_signal_seq.resize(Codeword_Seq.size());
						Decoding_info.rx_signal_seq.resize(Codeword_Seq.size());
					}
					else if (LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_Punc) {
						Turbo_LDPC_encoder_Punctured_Version(LinearBlockCode.G, Message_Seq, Codeword_Seq);
					}
					else if (LinearBlockCode.Code_number == Polar_with_Astar_middle_index) {
						Hybrid_Polar_Astar_Encoder(
							LinearBlockCode.G,
							LinearBlockCode.G_Inner,
							Message_Seq,
							Codeword_Seq,
							Decoding_info);
					}
					else Systematic_Linear_Block_Code_Encoder(LinearBlockCode.G, Message_Seq, Codeword_Seq);
					
					/*
					cout << "(" << LinearBlockCode.G._matrix.size() << "," << LinearBlockCode.G._matrix.at(0).size() << ")" << endl;
					for (int i = 0; i < 40; ++i) {
						cout << (int)Message_Seq.at(i);
					}
					cout << endl;
					for (int i = 0+64; i < 40+64; ++i) {
						cout << (int)Codeword_Seq.at(i);
					}*/
					/*
					for (int i = 0; i < LinearBlockCode.G.Row_number; ++i) {
						for (int j = 0; j < LinearBlockCode.G.Col_number-LinearBlockCode.G.Row_number; ++j) {
							cout << (int)LinearBlockCode.G._matrix[i][j+ LinearBlockCode.G.Row_number];
						}
						cout << endl;
					}*/
					/*
					for (int i = 0; i < LinearBlockCode.H.Row_number; ++i) {
						for (int j = 0; j < LinearBlockCode.H.Col_number/2; ++j) {
							cout << (int)LinearBlockCode.H._matrix[i][j+(LinearBlockCode.H.Col_number / 2)];
						}
						cout << endl;
					}*/
					//system("pause");
					//cout << Codeword_Length << "!!"<<endl;
					//system("pause");
					Decoding_info.message_seq = Message_Seq;
					Decoding_info.code_seq = Codeword_Seq;
					Decoding_info.snr = snr_dB;
					if (LinearBlockCode.Inner_Col_Number != 0) {   // P con & L con
						Decoding_info.Inner_message =(Codeword_Seq);
						Codeword_Seq.resize(LinearBlockCode.Inner_Col_Number);
						for (int i = 0; i < LinearBlockCode.Inner_Col_Number; ++i) {
							int temp = 0;
							for (int j = 0; j < Codeword_Length; j++) {
								temp ^= (Decoding_info.Inner_message.at(j)&LinearBlockCode.G_Inner._matrix[j][i]);
								//temp += Message_Seq.at(j)*LinearBlockCode.G._matrix[j][i];
							}
							Codeword_Seq.at(i) = temp;
						}
						Decoding_info.rx_signal_seq.resize(LinearBlockCode.Inner_Col_Number);
					}
					
					if (LinearBlockCode.Short_Inner_Col_Number != 0) {   // Hamming 255 con & RM 256 con & RM 256 con
						int NumberOfCodeBlock = (int)LinearBlockCode.Col_Number / LinearBlockCode.Short_Inner_Row_Number;
						//cout << NumberOfCodeBlock << endl;
						int NewCodeLength = LinearBlockCode.Col_Number*LinearBlockCode.Short_Inner_Col_Number / LinearBlockCode.Short_Inner_Row_Number;
						Decoding_info.Inner_message = (Codeword_Seq);
						//int NewCodeLength = LinearBlockCode.Col_Number*LinearBlockCode.Short_Inner_Col_Number / LinearBlockCode.Short_Inner_Row_Number;
						Codeword_Seq.resize(NewCodeLength);
						vector<__int8> M_Temp(LinearBlockCode.Short_Inner_Row_Number,0), C_Temp(LinearBlockCode.Short_Inner_Col_Number, 0);
						for (int i = 0; i < NumberOfCodeBlock; ++i) {
							for (int mesbit = 0; mesbit < LinearBlockCode.Short_Inner_Row_Number; mesbit++) {
								M_Temp.at(mesbit) = Decoding_info.Inner_message.at(i*LinearBlockCode.Short_Inner_Row_Number + mesbit);
							}
							if (LinearBlockCode.Code_number == Hamming_255_187_Astar_187_128) {
								Systematic_Linear_Block_Code_Encoder(LinearBlockCode.G_Inner, M_Temp, C_Temp);
							}
							else if (LinearBlockCode.Code_number == RM_256_192_Astar_192_128 || LinearBlockCode.Code_number == RM_256_168_Astar_168_128) {
								//cout << LinearBlockCode.G_Inner.Row_number << "," << LinearBlockCode.G_Inner.Col_number << endl;
								//system("pause");
								for (int mm = 0; mm < LinearBlockCode.G_Inner.Col_number; ++mm) {
									int temp = 0;
									for (int nn = 0; nn < LinearBlockCode.G_Inner.Row_number; nn++) {
										temp ^= (M_Temp.at(nn)&LinearBlockCode.G_Inner._matrix[nn][mm]);
										//temp += Message_Seq.at(j)*LinearBlockCode.G._matrix[j][i];
									}
									C_Temp.at(mm) = temp;
								}
							}
							for (int codebit = 0; codebit < LinearBlockCode.Short_Inner_Col_Number; codebit++) {
								Codeword_Seq.at(codebit + i * LinearBlockCode.Short_Inner_Col_Number) = C_Temp.at(codebit);
							}
						}
						Decoding_info.rx_signal_seq.resize(NewCodeLength);
					}
					/*
					for (int i = 0; i < NumberOfCodeBlock; ++i) {
						for (int j = 0; j < LinearBlockCode.Short_Inner_Col_Number; ++j) {
							cout << (int)Codeword_Seq.at(i*LinearBlockCode.Short_Inner_Col_Number + j);
						}
						cout << endl;
					}
					cout << endl << endl;
					system("pause");*/
					//cout << LinearBlockCode.Inner_Col_Number << endl;
					//system("pause");
					BPSK_Modulation(Codeword_Seq, Tx_signal_seq);
					//cout << Decoding_info.rx_signal_seq.size() << endl;
					Decoding_info.code_seq = Codeword_Seq;
					Decoding_info.message_seq = Message_Seq;
					//system("pause");
				}
			
				AWGN_Channel(
					AWGN_Mean,
					snr_dB,  // snr_dB
					LinearBlockCode.Code_Rate, // LinearBlockCode.Code_Rate System_CodeRate
					Tx_signal_seq,
					Decoding_info.rx_signal_seq,
					display_channel_info_flag);
				/*
				cout << endl;
				cout << "R':" << endl;
				for (int j = 0; j < 10; ++j) {
					cout << Decoding_info.rx_signal_seq.at(j) << "   ";
				}*/

				//cout << "A" << endl;

				// Error-free Test
				/*
				for (int i = 0; i < Decoding_info.rx_signal_seq.size(); ++i) {
					Decoding_info.rx_signal_seq.at(i) = Tx_signal_seq.at(i);
				}*/

				/*
				vector<__int8> Inter_M = Message_Seq;
				Block_Interleaver(Inter_M);
				cout << "M:" << endl;
				for (int j = 0; j < 5; ++j) {
					cout << (int)Message_Seq.at(j) << "   ";
				}
				cout << endl;
				cout << "M':" << endl;
				for (int j = 0; j < 5; ++j) {
					cout << (int)Inter_M.at(j) << "   ";
				}
				cout << endl;
				*/

				//cout << LinearBlockCode.Code_Rate << endl;
				//for(int i=0;i<CRow)
				//cout << Message_Seq.size() << endl;
				//cout << LinearBlockCode.Code_Rate << endl;
				//cout << Decoding_info.Constraint_i << "," << Decoding_info.CBC_FlippingBit << endl;

				if (LinearBlockCode.Punctured_Number != 0) {  // L  注入ghost error
					Decoding_info.estimated_codeword.resize(Codeword_Length + LinearBlockCode.Punctured_Number);
					Decoding_info.rx_signal_seq.insert(Decoding_info.rx_signal_seq.begin(), LinearBlockCode.Punctured_Number, 0);
				}
				if (LinearBlockCode.Code_number == 101) {   // P con
					//cout << "A";
					Inner_Polar_BP_Decoder_Outer_Astar_Decoder(LinearBlockCode.G_Inner, Decoding_info);
					
				}
				if (LinearBlockCode.Short_Inner_Col_Number != 0) {  // Hamming
					vector<double> R_Total(LinearBlockCode.Col_Number, 0);
					int NumberOfCodeBlock = LinearBlockCode.Col_Number / LinearBlockCode.Short_Inner_Row_Number;
					vector<double> R_Temp_Mes(LinearBlockCode.Short_Inner_Row_Number, 0),R_Temp_Code(LinearBlockCode.Short_Inner_Col_Number, 0);
					vector<__int8> Temp_CodeWord(LinearBlockCode.Short_Inner_Col_Number, 0);
				//	cout << Decoding_info.var << endl;
					double LLR_count = 0;
					double LLR_Temp = 0;
					for (int block = 0; block < NumberOfCodeBlock; block++) {
						int Bit;

						vector<double> RRR = R_Temp_Code;

						for (int j = 0; j < LinearBlockCode.Short_Inner_Col_Number; ++j) {
							R_Temp_Code.at(j) = Decoding_info.rx_signal_seq.at(j + block * LinearBlockCode.Short_Inner_Col_Number);
							LLR_count += abs(R_Temp_Code.at(j));
						}
						///
						if (LinearBlockCode.Code_number == Hamming_255_187_Astar_187_128) {
							SPA(LinearBlockCode.H, R_Temp_Code, Decoding_info.var);
							/*
							vector<int> Est_C(R_Temp_Code.size(), 0);
							for (int k = 0; k < R_Temp_Code.size(); ++k) {
								if (R_Temp_Code.at(k) > 0) Est_C.at(k) = 0;
								else Est_C.at(k) = 1;
							}
							int Fit = 1;
							for (int i = 0; i < LinearBlockCode.H.Row_number; ++i) {
								int temp = 0;
								for (int j = 0; j < LinearBlockCode.H.Col_number; ++j) {
									temp ^= (Est_C.at(j)&LinearBlockCode.H._matrix[i][j]);
								}
								if (temp == 1) {
									Fit = 0;
									break;
								}
							}
							if (Fit = 1) {
								for (int k = 0; k < R_Temp_Code.size(); ++k) {
									R_Temp_Code.at(k) *= 1.5;
								}
							}*/
						}
						else if (LinearBlockCode.Code_number == RM_256_192_Astar_192_128) {
							Majority_Soft_decoding(LinearBlockCode.G_Inner, R_Temp_Code, Decoding_info.var);
						}//RM_256_168_Astar_168_128
						else if (LinearBlockCode.Code_number == RM_256_168_Astar_168_128) {
							Majority_Soft_decoding_ver3(LinearBlockCode.G_Inner, LinearBlockCode.G_, R_Temp_Code, Decoding_info.var);
						}

						for (int i = 0; i < R_Temp_Code.size(); ++i)
						{
							R_Temp_Code.at(i) *= (Decoding_info.var / 2);
							/*
							if (R_Temp_Code.at(i) > 1000 || R_Temp_Code.at(i) < -1000) {
								cout << i << ", " << R_Temp_Code.at(i) << endl;
								for (int j = 0; j < RRR.size(); ++j) {
									cout << RRR.at(j) << " ";
								}
								cout << endl;
								for (int j = 0; j < RRR.size(); ++j) {
									cout << R_Temp_Code.at(j) << " ";
								}
								cout << endl;
							}*/
							//if (i == 0 && LinearBlockCode.Code_number == RM_256_192_Astar_192_128) R_Temp_Code.at(i) /= 3;
							//else if (i < 5 && i >= 1) R_Temp_Code.at(i)/=3;
							//R_Temp_Code.at(i) = atan(R_Temp_Code.at(i) / 2);
							/*
							if (i == 0) r0 += abs(R_Temp_Code.at(i));
							else if (i < 5 && i >= 1) r1 += abs(R_Temp_Code.at(i));
							else if (i < 11 && i >= 5)r2 += abs(R_Temp_Code.at(i));*/
						}
						/*
						int De_temp;
						cout << "#block " << block << ":            ";
						for (int i = 0; i < LinearBlockCode.Short_Inner_Row_Number; ++i) {
							if (R_Temp_Code.at(i) < 0) De_temp = 1;
							else De_temp = 0;
							if (De_temp != Decoding_info.Inner_message.at(i + block * LinearBlockCode.Short_Inner_Row_Number)) {
								++Inner_Error;
								cout << i << ",";
							}
						}
						cout << endl;*/
						/*
						cout << "C:" << endl;
						for (int i = 0; i < LinearBlockCode.Short_Inner_Row_Number; ++i) {
							cout << (int)Decoding_info.Inner_message.at(i + block * LinearBlockCode.Short_Inner_Row_Number) << " ";
						}
						cout << endl;
						cout << "Result:" << endl;
						for (int i = 0; i < LinearBlockCode.Short_Inner_Row_Number; ++i) {
							if (R_Temp_Code.at(i) < 0) cout << "1 ";
							else cout << "0 ";
						}
						cout << endl <<endl;*/
					
						for (int j = 0; j < LinearBlockCode.Short_Inner_Row_Number; ++j) {
							//cout << abs(R_Temp_Code.at(j)) << ", ";
							LLR_Temp += abs(R_Temp_Code.at(j));
							R_Total.at(j + block * LinearBlockCode.Short_Inner_Row_Number) = R_Temp_Code.at(j);// > 0 ? log10(R_Temp_Code.at(j)) : -log10(abs(R_Temp_Code.at(j)));
						}
						//cout << endl;
					}

					//cout << LLR_Temp << endl;
					AveLLR_Before += (2 * LLR_count / LinearBlockCode.Short_Inner_Col_Number / NumberOfCodeBlock / Decoding_info.var);
					//cout << (2 * LLR_Temp / LinearBlockCode.Short_Inner_Row_Number / NumberOfCodeBlock / Decoding_info.var) << endl;
					//Decoding_info.Ave_LLR += (2 * LLR_Temp / LinearBlockCode.Short_Inner_Row_Number / NumberOfCodeBlock / Decoding_info.var);
					//cout << Decoding_info.Ave_LLR << endl;
					Decoding_info.rx_signal_seq = R_Total;
					Decoding_info.rx_signal_seq.resize(LinearBlockCode.Col_Number);
					int counter = 0;
					for (int k = 0; k < Decoding_info.rx_signal_seq.size(); ++k) {
						if (Decoding_info.Inner_message.at(k) == 0) {
							LLR_Temp += Decoding_info.rx_signal_seq.at(k);
							++counter;
						}
					}
					LLR_Temp /= counter;
					Decoding_info.Ave_LLR += (2 * LLR_Temp / LinearBlockCode.Short_Inner_Row_Number / NumberOfCodeBlock / Decoding_info.var);
				}

				//   ******SPA回傳 Binary 或 LLR值 選擇********
				if (LinearBlockCode.Code_number > 90 && LinearBlockCode.Code_number < 100) Decoding_info.Return_Est_1_or_LLR_0 = 1;
				//   ******************************************

				if (This_is_LDPC) {
					LDPC.Sort_Sequence_Backward(Decoding_info.rx_signal_seq, LinearBlockCode.Permutation_Seq);
					vector<double> Rx = Decoding_info.rx_signal_seq;
					//if (Decoder == LDPC_VC_RBP) 
					LDPC_VC_RBP(LinearBlockCode.H, Decoding_info);
					//else SPA(LinearBlockCode.H, Decoding_info);
					
					if (Rx_diversity) {
						for (int i = 0; i < Decoding_info.rx_signal_seq.size(); ++i) {
							if (Rx.at(i)*Decoding_info.rx_signal_seq.at(i) > 0) Decoding_info.rx_signal_seq.at(i) += Rx.at(i);
						}
					}
					//Rx = Decoding_info.rx_signal_seq;
					//sort(Rx.begin(), Rx.end());
					//double threshold = abs(Rx.at(192));
					//for (int i = 0; i < Decoding_info.rx_signal_seq.size(); ++i) {
						//if (abs(Decoding_info.rx_signal_seq.at(i)) > threshold) Decoding_info.rx_signal_seq.at(i) > 0 ? 1000 : -1000;
					//}

					//Determine_Permutation(rx_signal, G, Sorted_G, Permutation);

					//SPA_Partial(LinearBlockCode.H, Decoding_info);
					if (Decoding_info.Return_Est_1_or_LLR_0)
					{
						LDPC.Sort_Sequence_Forward(Decoding_info.estimated_codeword, LinearBlockCode.Permutation_Seq);
						if (LinearBlockCode.Code_number == 93) {
							for (int i = 0; i < Message_Length; ++i) {
								Decoding_info.estimated_codeword.at(i) = Decoding_info.estimated_codeword.at(i + Message_Length);
							}
						}
					}
					else {
						LDPC.Sort_Sequence_Forward(Decoding_info.rx_signal_seq, LinearBlockCode.Permutation_Seq);
						Decoder(LinearBlockCode.G, Decoding_info); 
						if (LinearBlockCode.Code_number == 93) {
							for (int i = 0; i < Message_Length; ++i) {
								Decoding_info.estimated_codeword.at(i) = Decoding_info.estimated_codeword.at(i + Message_Length);
							}
						}
						/*
						cout << endl << "R(afterDecoding):" << endl;
						for (int i = 0; i < 10; ++i) {
							cout << Decoding_info.rx_signal_seq.at(i) << " ";
						}
						cout << endl << endl;
						*/
						if (LinearBlockCode.Code_number == 92)
							Decoding_info.estimated_codeword.erase(Decoding_info.estimated_codeword.begin(), Decoding_info.estimated_codeword.begin() + LinearBlockCode.H.Row_number);
						
					}
				}
				else if (LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC
					|| LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_ver2
					|| LinearBlockCode.Code_number == Turbo_320_128_with_2LDPC
					|| LinearBlockCode.Code_number == Turbo_256_128_RC_LDPC_192_128
					|| LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_Punc) {
					
					double LLRTemp = 0;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						LLRTemp += abs(Decoding_info.rx_signal_seq.at(i) * 2 / Decoding_info.var);
					}
					LLRTemp /= Message_Seq.size();
					AveLLR_Before += LLRTemp;
					//cout << "B:" << LLRTemp << ",";
					vector<double> Rx(Decoding_info.rx_signal_seq);

					//cout << This_is_LDPC << endl;
					if(LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_Punc) 
						Turbo_LDPC_decoder_punctured(LinearBlockCode.H_, 
							LinearBlockCode.H,
							Decoding_info, 
							LinearBlockCode.Permutation_Seq_, 
							LinearBlockCode.Permutation_Seq);
					else if (Repetition_Method) {
						Turbo_LDPC_Repitition_Method(LinearBlockCode.H_,
							LinearBlockCode.H,
							LinearBlockCode.G_Inner,
							Decoding_info,
							LinearBlockCode.Permutation_Seq_,
							LinearBlockCode.Permutation_Seq);
					}
					else {
						
						Turbo_LDPC_decoder(LinearBlockCode.H_,
							LinearBlockCode.H,
							LinearBlockCode.G_Inner,
							Decoding_info,
							LinearBlockCode.Permutation_Seq_,
							LinearBlockCode.Permutation_Seq);
					}
					//cout << "Hi" << endl;
					LLRTemp = 0;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						LLRTemp += abs(Decoding_info.rx_signal_seq.at(i) * 2 / Decoding_info.var);
						//if (Decoding_info.rx_signal_seq.at(i) * 2 / Decoding_info.var > 300) cout << Decoding_info.rx_signal_seq.at(i) * 2 / Decoding_info.var << ",";
						//cout << Decoding_info.rx_signal_seq.at(i) * 2 / Decoding_info.var << ", ";
					}
					//cout << endl <<endl;

					LLRTemp /= Message_Seq.size();
					Decoding_info.Ave_LLR += LLRTemp;

					/*
					for (int i = 0; i < LinearBlockCode.Permutation_Seq.size(); i++) {
						cout << LinearBlockCode.Permutation_Seq.at(i) << " ";
					}
					cout << endl;
					*/
					//cout <<"A: "<<  LLRTemp << endl;
					//cout << Decoding_info.Ave_LLR / Transmitted_Block << endl;
					if (Diff_Candidate_Test) {
						//A_star_PC_out_CBC_OSC_Verified_Function(G, decoding_info, Channel_LLR_CodeWord_Temp, Hard_Result, temp, Average_LLR);
						//cout << Decoding_info.Candidate_Metric.size() << endl;
						//int BestOne = 0;
						//double BestOneMetric = Decoding_info.Candidate_Metric.at(0);
						//for (int index = 0; index < Decoding_info.Candidate_Metric.size(); ++index) {
							//cout << Decoding_info.Candidate_Metric.at(index) << ",";
							//if (Decoding_info.Candidate_Metric.at(index) <= BestOneMetric) {
								//BestOne = index;
								//BestOneMetric = Decoding_info.Candidate_Metric.at(index);
							//}
							//cout << 
						//}
						//cout << endl << Decoding_info.Candidate_Metric.size();
						//cout <<endl<< BestOne << "(" << Decoding_info.Candidate_Metric.size() << ")," << BestOneMetric << endl;
						//Decoding_info.estimated_codeword = Decoding_info.Iter_Decode_Candidates_Est._matrix.at(BestOne);
						Decoding_info.rx_signal_seq.resize(192);
						cout << endl << Decoding_info.Candidate_Metric.size();
						if (Decoding_info.Candidate_Metric.size() < 40) {
							for (int j = 0; j < Decoding_info.rx_signal_seq.size(); ++j) 
							{
								Decoding_info.rx_signal_seq.at(j) = 0; // Decoding_info.rx_signal_seq.at(j)*(-1)*(Decoding_info.Candidate_Metric.size() - 1);
							}
							for (int i = 0; i < Decoding_info.Candidate_Metric.size(); ++i) {
								for (int j = 0; j < Decoding_info.rx_signal_seq.size(); ++j) {
									Decoding_info.rx_signal_seq.at(j) +=
										Decoding_info.Iterative_Decoding_Candidates._matrix[i][j];
								}
							}
						}
						else{

						}
						Decoder(LinearBlockCode.G, Decoding_info);

						vector<__int8> Mes(Message_Length, 0);
						for (int i = 0; i < Message_Length; i++) {
							Mes.at(i) = Decoding_info.estimated_codeword.at(i);
						}
						Block_DeInterleaver(Mes);
						for (int i = 0; i < Message_Length; i++) {
							Decoding_info.estimated_codeword.at(i) = Mes.at(i);
							//cout << (int)Mes.at(i) << ",";
						}
						//cout << endl;
					}
					if (Combine_Two_CodeWords) {
						//cout << "In";
						Decoder(LinearBlockCode.G_Inner, Decoding_info);
					}
					else if(Astar_after_TurboDecoding) // Astar_after_TurboDecoding
					{
						//cout << "I am here" << endl;
						vector<__int8> MM = Message_Seq;
						Block_Interleaver(MM);
						for (int i = 0; i < Message_Seq.size(); ++i) {
							Decoding_info.code_seq.at(i) = MM.at(i);
						}
						Decoding_info.code_seq.erase(Decoding_info.code_seq.begin() + Message_Seq.size(),
							Decoding_info.code_seq.begin() + LinearBlockCode.H_.Col_number);
						Decoding_info.code_seq.resize(LinearBlockCode.H.Col_number);

						Decoding_info.rx_signal_seq.erase(Decoding_info.rx_signal_seq.begin() + LinearBlockCode.H.Col_number, Decoding_info.rx_signal_seq.end());
						//Decoding_info.rx_signal_seq.resize(LinearBlockCode.H.Col_number);
						//cout << Decoding_info.rx_signal_seq.at(114) / Decoding_info.rx_signal_seq.at(145) << endl;
						/*
						double Avellr = 0;
						for (int i = 0; i < 128; i++) {
							Avellr += abs(Decoding_info.rx_signal_seq.at(i));
						}
						Avellr /= 128;
						if (Avellr < 30) {
							for (int i = 110; i < 128; i++) {
								Decoding_info.rx_signal_seq.at(i) /= 4;
							}
						}*/

						Decoder(LinearBlockCode.G, Decoding_info);

						/*
						MATRIX<__int8> Sorted_G(LinearBlockCode.G);
						vector <size_t> Location_Index(LinearBlockCode.G.Col_number, 0);
						MATRIX<double> Metric_Table(2, Decoding_info.rx_signal_seq.size());
						Pre_Procedure(Decoding_info.rx_signal_seq, LinearBlockCode.G, Sorted_G, Location_Index, Metric_Table, Decoding_info);
						vector<__int8> MM(Message_Seq.size()),CC(Decoding_info.rx_signal_seq.size());
						for (int i = 0; i < Message_Seq.size(); ++i) {
							if (Decoding_info.Sorted_R.at(i) > 0) MM.at(i) = 0;
							else MM.at(i) = 1;
						}
						Systematic_Linear_Block_Code_Encoder(Sorted_G, MM, CC);
						Desort_Function(Location_Index, CC, Decoding_info.estimated_codeword);
						*/

						vector<__int8> Mes(Message_Length, 0);
						for (int i = 0; i < Message_Length; i++) {
							Mes.at(i) = Decoding_info.estimated_codeword.at(i);
						}
						Block_DeInterleaver(Mes);
						for (int i = 0; i < Message_Length; i++) {
							Decoding_info.estimated_codeword.at(i) = Mes.at(i);
						}
					}
				}
				else if (LinearBlockCode.Code_number == RandomCode_with_LDPC) {
					Decoding_info.Return_Est_1_or_LLR_0 = 1;
					SPA(LinearBlockCode.H, Decoding_info);
					A_star_PC_out_CBC_OSC_Verified(LinearBlockCode.G, Decoding_info);
				}
				else if (LinearBlockCode.Code_number == Polar_with_Astar_middle_index) {
				Hybrid_Polar_Astar_Decoder(
					LinearBlockCode.G,
					LinearBlockCode.G_Inner,
					Decoding_info);
				}
				else {
			     	//Decoder = A_star_PC_out_CBC_OSC_Verified;
				    //Operation_Parameter.Decoder_Version == A_Star_PC_out_CBC_OSC_Verified;
					Decoder(LinearBlockCode.G, Decoding_info);  // 一般 Astar decoder!
					//Operation_Parameter.Decoder_Version == Polar_Code_BP_Decoder;
				}
				/*
				cout <<endl<< "R(afterDecoding):" << endl;
				for (int i = 0; i < 10; ++i) {
					cout << Decoding_info.rx_signal_seq.at(i) << " ";
				}
				cout << endl << endl;
				*/

			//	Decoding_info.estimated_codeword.erase(Decoding_info.estimated_codeword.begin(), Decoding_info.estimated_codeword.begin() + 64);
				
				if (LinearBlockCode.Punctured_Number != 0) {
					if (Decoding_info.Return_Est_1_or_LLR_0) {
						Decoding_info.estimated_codeword.erase(Decoding_info.estimated_codeword.begin(), Decoding_info.estimated_codeword.begin() + LinearBlockCode.H.Row_number);
					}
					else {
						int Add = LinearBlockCode.Punctured_Number + LinearBlockCode.Col_Number - LinearBlockCode.Row_Number;
						Decoding_info.rx_signal_seq.erase(Decoding_info.rx_signal_seq.begin(), Decoding_info.rx_signal_seq.begin() + Add);
						//cout << Decoding_info.rx_signal_seq.size() << endl;
						Decoding_info.estimated_codeword.resize(LinearBlockCode.Col_Number);
						Decoder(LinearBlockCode.G, Decoding_info);
					}
				}
				
			/*
				cout << "-----------------Start------------------------" << endl;
				cout << "M1:" << endl;
				for (int i = 0; i < 30; ++i) {
					cout << (int)Message_Seq.at(i) << " ";
				}
				cout << endl;
				cout << "M2:" << endl;
				for (int i = 64; i < 94; ++i) {
					cout << (int)Message_Seq.at(i) << " ";
				}
				cout << endl;

				cout << "R1:" << endl;
				for (int i = 0; i < 30; ++i) {
					cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
				}
				cout << endl;
				cout << "R2:" << endl;
				for (int i = 64; i < 94; ++i) {
					cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
				}
				cout << endl << endl;
				cout << "-------------------End----------------------" << endl << endl;*/
				/*
				cout << "M: " << endl;
				for (int i = 0; i < Message_Length; ++i) {
					cout << (int)Message_Seq.at(i) << " ";
				}
				cout << endl;
				cout << "Est: " << endl;
				for (int i = 0; i < Message_Length; ++i) {
					cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
				}
				cout << endl << endl << endl;
				*/

				if (CRC_Check) Message_Seq.erase(Message_Seq.end() - CRC_length, Message_Seq.end());

				{

					Total_Avg_STE += Decoding_info.STE;
					Total_Avg_Com += Decoding_info.COM;
					Total_Avg_DM_STE += Decoding_info.DM_STE;
					Total_Avg_DM_Com += Decoding_info.DM_COM;
					Total_Avg_DN += Decoding_info.DeletedNode_counter;
					Total_Avg_CandidateCodeWord += Decoding_info.CandidateCodeWord;

					// 20190711 new method for calculating the number of message errors  
					// do the pairwise addition over GF(2) for estimated message and transmitted message 
					//cout << Message_Seq.size() << endl;
					transform(
						Decoding_info.estimated_codeword.begin(),
						Decoding_info.estimated_codeword.begin() + Message_Seq.size(),  // Message_Seq.size()  Message_Length
						Message_Seq.begin(),  // Message_Seq 
						Error_Seq.begin(),
						std::bit_xor<__int8>());

					// sum the number of bit 1 
					Temp_ErrorBit_Counter = accumulate(Error_Seq.begin(), Error_Seq.end(), 0);
					//cout << Transmitted_Block << ":" << Temp_ErrorBit_Counter << endl;

					ErrorBit_Counter += Temp_ErrorBit_Counter;
					if (Temp_ErrorBit_Counter != 0)
						++ErrorBlock_Counter;

					if (Temp_ErrorBit_Counter != 0) cout << "ERROR: "<< Temp_ErrorBit_Counter << endl;
					/*
					for (int i = 0; i < Decoding_info.rx_signal_seq.size(); ++i) {
						cout << setprecision(2) << Decoding_info.rx_signal_seq.at(i) << " ";
					}
					system("pause");
					*/
					//20190712 Computation of ML Lower Bound 
					/*
					if (ML_LowerBound(
						Decoding_info.rx_signal_seq,
						Tx_signal_seq,
						Decoding_info.estimated_codeword)) {

						++ML_LB_ErrorBlock_Counter;
						ML_LB_ErrorBit_Counter += Temp_ErrorBit_Counter;
					}*/
				}
				//##CANCELED##//
				/*
				if (Temp_ErrorBit_Counter != 0) {
					cout << "***ERROR: " << Temp_ErrorBit_Counter << endl;
					for (int i = 0; i < Message_Length; ++i) {
						if (Message_Seq.at(i) != Decoding_info.estimated_codeword.at(i)) cout << i << ", ";
					}
					cout << endl;
					//cout << "R0: " << r0 / Transmitted_Block << ", R1: " << r1 / Transmitted_Block / 4 << ", R2: " << r2 / Transmitted_Block / 6 << endl;

				}*/
				//cout << endl;

				//Operation_Parameter.BlockStep = 1;
				//cout << "C";
				//
				// When each time of the counter arrives the configured steps of the block, 
				// show the real-time decoding results.  
				if (//(Block_Counter >= Operation_Parameter.BlockStep) 
					(((size_t)Transmitted_Block % Operation_Parameter.BlockStep) == 0) ||
					(Transmitted_Block == Operation_Parameter.block_number)) {

					// Statistical decoding parameters
						{
							//Block_Counter = 0;
							// error rate
							BER = ErrorBit_Counter / (Message_Seq.size()*Transmitted_Block);  // Message_Seq
							BLER = ErrorBlock_Counter / Transmitted_Block;

							// ML lower bound
							ML_LB_BER = ML_LB_ErrorBit_Counter / (Message_Seq.size()*Transmitted_Block); // Message_Seq
							ML_LB_BLER = ML_LB_ErrorBlock_Counter / Transmitted_Block;

							// statistics of decoder
							Avg_STE = Total_Avg_STE / Transmitted_Block;
							Avg_Com = Total_Avg_Com / Transmitted_Block;
							Avg_CandidateCodeWord = Total_Avg_CandidateCodeWord / Transmitted_Block;

							// statistics of decoder used in deleting mechanism
							Avg_DM_STE = Total_Avg_DM_STE / Transmitted_Block;
							Avg_DM_Com = Total_Avg_DM_Com / Transmitted_Block;
							Avg_DN = Total_Avg_DN / Transmitted_Block;
						}

						// Yin
						/*
						if ((size_t)Transmitted_Block % 100000 == 0) {
							cout << "#" << (int)Transmitted_Block << ":" << endl;
							for (int i = 0; i < 20; ++i) {
								for (int j = 0; j < 20; ++j) {
									cout << setw(8) << Decoding_info.FlipRecorder[i][j];
								}
								cout << endl;
							}
							cout << endl;
						}



						/*
						if ((size_t)Transmitted_Block % 1000 == 0) {
							for (int i = 0; i < 101; ++i) {
								cout << setw(5) << Decoding_info.Accumulate_Break.at(i) << " ";
								if (i % 20 == 0)cout << endl;
							}
							cout << endl << endl;
						}
						*/
						// Show the real-time statistics. 


						printf(
							"\r Block #%.0f (%.2f %%)"
							", Eb #%.0f"
							", BER：%.2e"
							", EB #%.0f"
							", BLER：%.2e"
							//", ML_Eb # %.0f "
							//", ML_BER : %.2e"
							//", ML_EB # %.0f "
							//", ML_BLER : %.2e"
							", Avg.STE = %.2e"
							", Avg.COM = %.2e"
							//", Avg.CD = %.2e"
							//", Worst.CD = %d"
							//", MUS = %d"
							", Worst.STE = %.2e"
							", Worst.COM = %.2e"
							", Ave_LLR_Before = %.3f"
							", Aver_LLR = %.3f"
							", ML(BER): %.2e"
							", ML(BLER): %.2e"
							//", Ave_Error = %.2f"
							, Transmitted_Block
							, (Transmitted_Block * 100. / Operation_Parameter.block_number)
							, ErrorBit_Counter
							, BER
							, ErrorBlock_Counter
							, BLER
							//, ML_LB_ErrorBit_Counter
							//, ML_LB_BER
							//, ML_LB_ErrorBlock_Counter
							//, ML_LB_BLER
							, Avg_STE
							, Avg_Com
							//, Avg_CandidateCodeWord
							//, Decoding_info.Worst_Case_Candidate
							//, Decoding_info.MaxUsed_Stack
							, Decoding_info.Worst_Case_STE
							, Decoding_info.Worst_Case_COM
							, AveLLR_Before / Transmitted_Block
							, Decoding_info.Ave_LLR / Transmitted_Block
							, ML_LB_BER
							, ML_LB_BLER
							//, Inner_Error / Transmitted_Block / Codeword_Length
						);/*
						  ", ML(BER): %.2e"
							", ML(BLER): %.2e"
							, ML_LB_BER
							, ML_LB_BLER
						  */
						/*
						printf(
							"\r Block #%.0f (%.2f %%)"
							", Eb #%.0f"
							", BER：%.2e"
							", EB #%.0f"
							//", ML_Eb # %.0f "
							//", ML_BER : %.2e"
							//", ML_EB # %.0f "
							//", ML_BLER : %.2e"
							", Avg.STE = %.2e"
							", C0: %d"
							", C1: %d"
							", C2: %d"
							", C3: %d"
							", C4: %d"
							", C5: %d  "
							", D1: %d"
							", D2: %d"
							", D3: %d"
							", D4: %d"
							", D5: %d"
							", CBC_Candidate: %d"
							", Total_Candidate: %d"
							", CBC_STE: %d"
							, Transmitted_Block
							, (Transmitted_Block * 100. / Operation_Parameter.block_number)
							, ErrorBit_Counter
							, BER
							, ErrorBlock_Counter
							//, ML_LB_ErrorBit_Counter
							//, ML_LB_BER
							//, ML_LB_ErrorBlock_Counter
							//, ML_LB_BLER
							, Avg_STE
							, Decoding_info.Dz_0_Number
							, Decoding_info.Dz_1_Number
							, Decoding_info.Dz_2_Number
							, Decoding_info.Dz_3_Number
							, Decoding_info.Dz_4_Number
							, Decoding_info.Dz_5_Number
							, Decoding_info.D1
							, Decoding_info.D2
							, Decoding_info.D3
							, Decoding_info.D4
							, Decoding_info.D5
							, (Decoding_info.CBC_Candidate / (int)Transmitted_Block)
							, (Decoding_info.TotalCounter / (int)Transmitted_Block)
							, (Decoding_info.CBC_STE / (int)Transmitted_Block / Message_Length)
						);*/

						double Execution_Time((clock() - start_time_for_SNR) / 1000.);

						// Store the real-time statistics of decoder 
						{
							ClearFile(Realtime_Data);
							WriteFile(Realtime_Data, snr_dB);
							WriteFile(Realtime_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
							WriteFile(Realtime_Data, BER);
							WriteFile(Realtime_Data, BLER);
							WriteFile(Realtime_Data, Avg_STE);
							WriteFile(Realtime_Data, Avg_Com);
							WriteFile(Realtime_Data, Avg_CandidateCodeWord);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_STE);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_COM);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_Candidate);
							WriteFile(Realtime_Data, Decoding_info.MaxUsed_Stack);
							WriteFile(Realtime_Data, "-----------");
							//WriteFile(Realtime_Data, Avg_DM_STE);
							//WriteFile(Realtime_Data, Avg_DM_Com);
							//WriteFile(Realtime_Data, Avg_DN);
							WriteFile(Realtime_Data, ML_LB_BER);
							WriteFile(Realtime_Data, ML_LB_BLER);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, (Decoding_info.CBC_Candidate / Transmitted_Block));
							WriteFile(Realtime_Data, Decoding_info.TotalCounter / Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.CBC_STE / Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, (Decoding_info.DM_STE*Decoding_info.CBC_length) / Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, "\n\n");
						}
						/*
						{
							ClearFile(Realtime_Data);
							WriteFile(Realtime_Data, snr_dB);
							WriteFile(Realtime_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
							WriteFile(Realtime_Data, BER);
							WriteFile(Realtime_Data, BLER);
							WriteFile(Realtime_Data, Avg_STE);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.Dz_0_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_1_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_2_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_3_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_4_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_5_Number);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.D1);
							WriteFile(Realtime_Data, Decoding_info.D2);
							WriteFile(Realtime_Data, Decoding_info.D3);
							WriteFile(Realtime_Data, Decoding_info.D4);
							WriteFile(Realtime_Data, Decoding_info.D5);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.CBC_Candidate / (int)Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.TotalCounter / (int)Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.CBC_STE / (int)Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, "\n\n");
						}*/
						if (ErrorBlock_Counter >= Operation_Parameter.ErrorBlock_Thr) break;
					    //if (Transmitted_Block >= 1000000) break;
				}
			}
			/*
			cout << endl << endl << "Result:" << endl;
			for (int x = 0; x < 9; ++x) {
				for (int y = 0; y < 13; ++y) {
					cout << Decoding_info.Amplitude_Flipping_number._matrix[x][y] / Decoding_info.Amplitude_Flipping_number._matrix[x][13] << "  ";
				}
				cout << endl;
				if (x % 3 == 2) cout << endl;
			}
			cout << endl;
			cout << Decoding_info.Amplitude_Flipping_number._matrix[0][13] / Transmitted_Block << endl;
			cout << Decoding_info.Amplitude_Flipping_number._matrix[3][13] / Transmitted_Block << endl;
			cout << Decoding_info.Amplitude_Flipping_number._matrix[6][13] / Transmitted_Block << endl;
			*/

			if (EXIT_Chart_Record) {
				string C1_priori = "C1_priori.txt";
				string C1_extrinsic = "C1_extrinsic.txt";
				string C2_priori = "C2_priori.txt";
				string C2_extrinsic = "C2_extrinsic.txt";

				ClearFile(C1_priori);
				ClearFile(C1_extrinsic);
				ClearFile(C2_priori);
				ClearFile(C2_extrinsic);

				cout << endl << Decoding_info.LLR_1_priori.size() << "," << Decoding_info.LLR_1_extrinsic.size() << ","
					<< Decoding_info.LLR_2_priori.size() << "," << Decoding_info.LLR_2_extrinsic.size() << endl;

				for (int i = 0; i < Decoding_info.LLR_1_extrinsic.size(); ++i) {
					WriteFile(C1_priori, Decoding_info.LLR_1_priori.at(i));
					WriteFile(C1_extrinsic, Decoding_info.LLR_1_extrinsic.at(i));
					WriteFile(C2_priori, Decoding_info.LLR_2_priori.at(i));
					WriteFile(C2_extrinsic, Decoding_info.LLR_2_extrinsic.at(i));
					if (i % 1000 == 0) cout << i << endl;
				}
				cout << endl << "END" << endl;
				system("pause");
			}
			/*
			string LLR_0 = "LLR_0.txt";
			string LLR_1 = "LLR_1.txt";
			ClearFile(LLR_0);
			ClearFile(LLR_1);
			for (int i = 0; i < 2000; ++i) {
				WriteFile(LLR_0, LLR_Distribution_0.at(i));
				WriteFile(LLR_1, LLR_Distribution_1.at(i));
			}*/
		}// for loop : Transmitted Block

		
		//      ************************************                  general A*                   ********************************************* 
		else {
			Decoding_info.snr = snr_dB;
			Decoding_info.var = (Es / (2.*LinearBlockCode.Code_Rate)) * pow(10., (-(snr_dB / 10.)));
			Decoding_info.CBC_length = Decoding_info.Control_Level - Message_Length;
			//cout << snr_dB<<":  "<<Decoding_info.var << endl;
			//system("pause");
			vector<double> LLR_Distribution_0(2000, 0);
			vector<double> LLR_Distribution_1(2000, 0);

			Decoding_info.Return_Est_1_or_LLR_0 = 1;
			long double r0 = 0,r1 = 0, r2 = 0, r3 = 0;
			long double BERr0 = 0, BERr1 = 0, BERr2 = 0, BERr3 = 0;
			//vector<int> Error_Array(64, 0);
			for (
				Transmitted_Block = 1;
				Transmitted_Block <= Operation_Parameter.block_number;
				++Transmitted_Block) {
				Decoding_info.initial_para();
				if (RandomMessage_flag) {
					Generate_Message(Message_Length, Message_Seq);
										
						// For LLR Distribution Test
						/*
						for (int i = 0; i < Message_Seq.size(); ++i) {
							Message_Seq.at(i) = 0;
						}*/
					
					Systematic_Linear_Block_Code_Encoder(LinearBlockCode.G, Message_Seq, Codeword_Seq);
					
					Decoding_info.message_seq = Message_Seq;
					BPSK_Modulation(Codeword_Seq, Tx_signal_seq);
					Decoding_info.code_seq = Codeword_Seq;
					/*
					for (int i = 0; i < Message_Seq.size(); ++i) {
						Tx_signal_seq.at(i) *= 2;
					}
					for (int i = Message_Seq.size(); i < Codeword_Seq.size(); ++i) {
						Tx_signal_seq.at(i) *= 0;
					}*/
				}
			
				AWGN_Channel(
					AWGN_Mean,
					snr_dB,  // snr_dB
					LinearBlockCode.Code_Rate, // LinearBlockCode.Code_Rate System_CodeRate
					Tx_signal_seq,
					Decoding_info.rx_signal_seq,
					display_channel_info_flag);
				
				/*
				for (int i = 0; i < Message_Seq.size(); ++i) {
					Decoding_info.rx_signal_seq.at(i) *= 2 / Decoding_info.var;
				}
				for (int i = Message_Seq.size(); i < Codeword_Seq.size(); ++i) {
					Decoding_info.rx_signal_seq.at(i) *= 0 / Decoding_info.var;
				}*/

				//Decoder(LinearBlockCode.H, Decoding_info);
				Decoding_info.Code_Rate = LinearBlockCode.Code_Rate;
				Decoder(LinearBlockCode.G, Decoding_info);
				
				
				//Majority_Hard_decoding(LinearBlockCode.G, LinearBlockCode.G_, Decoding_info);

				//Majority_Soft_decoding(LinearBlockCode.G, Decoding_info.rx_signal_seq, Decoding_info.var);
				//Decoding_info.estimated_codeword.erase(Decoding_info.estimated_codeword.begin(), Decoding_info.estimated_codeword.begin() + 4);
				
				/*
				// RM(64,42) Decoder
				
				Majority_Soft_decoding_ver3(LinearBlockCode.G, LinearBlockCode.G_, Decoding_info.rx_signal_seq, Decoding_info.var);
				for (int i = 0; i < Message_Seq.size(); ++i) {
					if (Decoding_info.rx_signal_seq.at(i) < 0) Decoding_info.estimated_codeword.at(i) = 1;
					else Decoding_info.estimated_codeword.at(i) = 0;
				}*/
				

				//cout << Message_Seq.size() << endl;

			    // For LLR Distribution Test
				/*
				int Temp;
				for (int i = 0; i < Message_Seq.size(); ++i) {
					Temp = ROUND_2_INT(( Decoding_info.rx_signal_seq.at(i) + 100) * 10);
					if (Temp < 0) Temp = 0;
					else if (Temp >= 2000) Temp = 1999;
					LLR_Distribution_1.at(Temp)++;
				}
				cout << endl;
				*/
				
				/*
				r0 += (Decoding_info.rx_signal_seq.at(0));
				for (int i = 1; i < 7; ++i) {
					r1 += (Decoding_info.rx_signal_seq.at(i));
				}
				for (int i = 7; i < 22; ++i) {
					r2 += (Decoding_info.rx_signal_seq.at(i));
				}
				for (int i = 22; i < 42; ++i) {
					r3 += (Decoding_info.rx_signal_seq.at(i));
				}

				
				for (int i = 0; i < Message_Seq.size(); ++i) {
					if (i == 0 && Decoding_info.estimated_codeword.at(i) != Message_Seq.at(i)) ++BERr0;
					else if (i < 7 && Decoding_info.estimated_codeword.at(i) != Message_Seq.at(i)) ++BERr1;
					else if (i < 22 && Decoding_info.estimated_codeword.at(i) != Message_Seq.at(i)) ++BERr2;
					else if(i < 42 && Decoding_info.estimated_codeword.at(i) != Message_Seq.at(i)) ++BERr3;
				}*/
				
				
				/*
				cout << "block #" << Transmitted_Block << ":                ";
				for (int i = 0; i < Message_Seq.size(); ++i) {
					if (Decoding_info.estimated_codeword.at(i) != Message_Seq.at(i)) cout << i << ",";
				}
				cout << endl;
				*/
				/*
				cout << "C:" << endl;
				for (int i = 0; i < Message_Seq.size(); ++i) {
					cout << (int)Message_Seq.at(i) << " ";
				}
				cout << endl;
				cout << "R:" << endl;
				for (int i = 0; i < Decoding_info.rx_signal_seq.size(); ++i) {
					if(Decoding_info.rx_signal_seq.at(i)>0) cout << "0 ";
					else cout << "1 ";
				}
				cout << endl;
				cout << "Result:" << endl;
				for (int i = 0; i < Message_Seq.size(); ++i) {
					cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
				}
				cout << endl << endl;
				*/
				
				{
					Total_Avg_STE += Decoding_info.STE;
					Total_Avg_STE_1 = Decoding_info.STE_1;
					Total_Avg_STE_2 = Decoding_info.STE_2;
					Total_Avg_STE_3 = Decoding_info.STE_3;
					Total_Avg_Com += Decoding_info.COM;
					Total_Avg_Com_1 += Decoding_info.COM_1;
					Total_Avg_DM_STE += Decoding_info.DM_STE;
					Total_Avg_DM_Com += Decoding_info.DM_COM;
					Total_Avg_DN += Decoding_info.DeletedNode_counter;
					Total_Avg_CandidateCodeWord += Decoding_info.CandidateCodeWord;
					Total_Avg_Alpha += Decoding_info.New_OSC_Alpha;
					Total_Avg_Binary_STE += Decoding_info.Binary_STE;
					Total_Updated_Best_Node += Decoding_info.num_Best_Node_Update;
					Total_Deleted_Candidate += Decoding_info.num_Deleted_Candidate_in_Stack;
					
					

					// 20190711 new method for calculating the number of message errors  
					// do the pairwise addition over GF(2) for estimated message and transmitted message 
					
					transform(
						Decoding_info.estimated_codeword.begin(),
						Decoding_info.estimated_codeword.begin() + Message_Seq.size(),  // Message_Seq.size()  Message_Length
						Message_Seq.begin(),
						Error_Seq.begin(),
						std::bit_xor<__int8>());
					/*
					for (int i = 0; i < Message_Seq.size(); ++i) {
						cout << (int)Message_Seq.at(i) << " ";
					}
					cout << endl; 
					int Z0 = 0;
					for (int i = 0; i < Codeword_Seq.size(); ++i) {
						cout << (int)Codeword_Seq.at(i) << " ";
						if (Codeword_Seq.at(i) == 0) ++Z0;
					}
					cout << endl;
					cout << Z0 << endl;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
					}
					cout << endl << endl;*/

					// sum the number of bit 1 
					Temp_ErrorBit_Counter = accumulate(Error_Seq.begin(), Error_Seq.end(), 0);
					
					
					ErrorBit_Counter += Temp_ErrorBit_Counter;
					if (Temp_ErrorBit_Counter != 0)
						++ErrorBlock_Counter;
				
					//20190712 Computation of ML Lower Bound 
					
					if (ML_LowerBound(
						Decoding_info.rx_signal_seq,
						Tx_signal_seq,
						Decoding_info.estimated_codeword)) {

						++ML_LB_ErrorBlock_Counter;
						ML_LB_ErrorBit_Counter += Temp_ErrorBit_Counter;
					}
				}
				/*
				Decoder = A_star_PC_out_CBC_OSC_Adaptive_i_NEW;
				Decoder(LinearBlockCode.G, Decoding_info);

				Total_Avg_STE_2 = Decoding_info.STE_2;
				Total_Avg_Com_2 += Decoding_info.COM_2;

				transform(
					Decoding_info.estimated_codeword.begin(),
					Decoding_info.estimated_codeword.begin() + Message_Seq.size(),  // Message_Seq.size()  Message_Length
					Message_Seq.begin(),
					Error_Seq.begin(),
					std::bit_xor<__int8>());

				Temp_ErrorBit_Counter = accumulate(Error_Seq.begin(), Error_Seq.end(), 0);

				ErrorBit_Counter_1 += Temp_ErrorBit_Counter;
				*/
				
				//##CANCELED##//
				/*
				if (Temp_ErrorBit_Counter != 0) {
					cout << "ERROR: " << Temp_ErrorBit_Counter;
					cout << endl;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						cout << (int)Message_Seq.at(i) << " ";
					}
					cout << endl;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						cout << (int)Decoding_info.estimated_codeword.at(i) << " ";
					}
					cout << endl << endl;
				}*/

				//
				// When each time of the counter arrives the configured steps of the block, 
				// show the real-time decoding results.  

				if (//(Block_Counter >= Operation_Parameter.BlockStep) 
					(((size_t)Transmitted_Block % Operation_Parameter.BlockStep) == 0) ||
					(Transmitted_Block == Operation_Parameter.block_number)) {
					/*
					cout << endl;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						cout << setw(6) << i;
					}
					cout << endl;
					for (int i = 0; i < Message_Seq.size(); ++i) {
						printf("%.3f ", Error_Array.at(i) / Transmitted_Block);
					}
					cout << endl << endl;
					*/
					
					//cout << endl << r0 / Transmitted_Block << ", " << r1 / 6 / Transmitted_Block << ", " << r2 / 15 / Transmitted_Block << ", " << r3 / 20 / Transmitted_Block << endl;
					//printf("r0: %.4f, r1: %.4f, r2: %.4f, r3: %.4f\n", BERr0 / Transmitted_Block, BERr1 / 6 / Transmitted_Block, BERr2 / 15 / Transmitted_Block, BERr3 / 20 / Transmitted_Block);
					// Statistical decoding parameters
					
						{
							//Block_Counter = 0;
							// error rate
							BER = ErrorBit_Counter / (Message_Seq.size()*Transmitted_Block);
							BLER = ErrorBlock_Counter / Transmitted_Block;

							// ML lower bound
							ML_LB_BER = ML_LB_ErrorBit_Counter / (Message_Seq.size()*Transmitted_Block);
							ML_LB_BLER = ML_LB_ErrorBlock_Counter / Transmitted_Block;

							// statistics of decoder
							Avg_STE = Total_Avg_STE / Transmitted_Block;
							Avg_STE_1 = Total_Avg_STE_1 / Transmitted_Block;
							Avg_STE_2 = Total_Avg_STE_2 / Transmitted_Block;
							Avg_STE_3 = Total_Avg_STE_3 / Transmitted_Block;
							Avg_Com = Total_Avg_Com / Transmitted_Block;
							Avg_Com_1 = Total_Avg_Com_1 / Transmitted_Block;
							Avg_Com_2 = Total_Avg_Com_2 / Transmitted_Block;
							Avg_CandidateCodeWord = Total_Avg_CandidateCodeWord / Transmitted_Block;
							Avg_Alpha = Total_Avg_Alpha / Transmitted_Block;
							Avg_Binary_STE = Total_Avg_Binary_STE / Transmitted_Block;
							

							// statistics of decoder used in deleting mechanism
							Avg_DM_STE = Total_Avg_DM_STE / Transmitted_Block;
							Avg_DM_Com = Total_Avg_DM_Com / Transmitted_Block;
							Avg_DN = Total_Avg_DN / Transmitted_Block;

							//
							Avg_Deleted_Candidate = Total_Deleted_Candidate / Total_Updated_Best_Node;
						}
						
						// Yin
						/*
						if ((size_t)Transmitted_Block % 100000 == 0) {
							cout << "#" << (int)Transmitted_Block << ":" << endl;
							for (int i = 0; i < 20; ++i) {
								for (int j = 0; j < 20; ++j) {
									cout << setw(8) << Decoding_info.FlipRecorder[i][j];
								}
								cout << endl;
							}
							cout << endl;
						}



						/*
						if ((size_t)Transmitted_Block % 1000 == 0) {
							for (int i = 0; i < 101; ++i) {
								cout << setw(5) << Decoding_info.Accumulate_Break.at(i) << " ";
								if (i % 20 == 0)cout << endl;
							}
							cout << endl << endl;
						}
						*/
						// Show the real-time statistics. 


						printf(
							"\r Block #%.0f (%.2f %%)"
							", Eb #%.0f"
							", BER：%.2e"
							//", BER67：%.2e"
							", EB #%.0f"
							", BLER：%.2e"
							//", ML_Eb # %.0f "
							//", ML_BER : %.2e"
							//", ML_EB # %.0f "
							//", ML_BLER : %.2e"
							", Avg.STE = %.2e"
							//", Avg.STE66 = %.2e"
							//", Avg.STE67 = %.2e"
							", Avg.COM = %.2e"
							//", Avg.COM66 = %.2e"
							//", Avg.COM67 = %.2e"
							", Avg.Binary STE = %.2e"
							//", Avg.CD = %.2e"
							//", Worst.CD = %d"
							//", MUS = %d"
							", Worst.STE = %.2e"
							", Worst.COM = %.2e"
							", Worst.OSC Ratio = %.2e"						//PoHan
							", CBC.STE = %.2e"
							", LLR: %.3f"
							", ML(BER): %.2e"
							", ML(BLER): %.2e"
							, Transmitted_Block
							, (Transmitted_Block * 100. / Operation_Parameter.block_number)
							, ErrorBit_Counter
							, BER
							//, BER_1
							, ErrorBlock_Counter
							, BLER
							//, ML_LB_ErrorBit_Counter
							//, ML_LB_BER
							//, ML_LB_ErrorBlock_Counter
							//, ML_LB_BLER
							, Avg_STE
							//, Avg_STE_1
							//, Avg_STE_2
							, Avg_Com
							//, Avg_Com_1
							//, Avg_Com_2
							, Avg_Binary_STE
							//, Avg_CandidateCodeWord
							//, Decoding_info.Worst_Case_Candidate
							//, Decoding_info.MaxUsed_Stack
							, Decoding_info.Worst_Case_STE
							, Decoding_info.Worst_Case_COM
							, Decoding_info.Worst_OSC_Ratio							//PoHan
							, Decoding_info.CBC_STE / Transmitted_Block / Message_Length
							, Decoding_info.Ave_LLR / Transmitted_Block
							, ML_LB_BER
							, ML_LB_BLER
						);
						/*
						printf(
							"\r Block #%.0f (%.2f %%)"
							", Eb #%.0f"
							", BER：%.2e"
							", EB #%.0f"
							//", ML_Eb # %.0f "
							//", ML_BER : %.2e"
							//", ML_EB # %.0f "
							//", ML_BLER : %.2e"
							", Avg.STE = %.2e"
							", C0: %d"
							", C1: %d"
							", C2: %d"
							", C3: %d"
							", C4: %d"
							", C5: %d  "
							", D1: %d"
							", D2: %d"
							", D3: %d"
							", D4: %d"
							", D5: %d"
							", CBC_Candidate: %d"
							", Total_Candidate: %d"
							", CBC_STE: %d"
							, Transmitted_Block
							, (Transmitted_Block * 100. / Operation_Parameter.block_number)
							, ErrorBit_Counter
							, BER
							, ErrorBlock_Counter
							//, ML_LB_ErrorBit_Counter
							//, ML_LB_BER
							//, ML_LB_ErrorBlock_Counter
							//, ML_LB_BLER
							, Avg_STE
							, Decoding_info.Dz_0_Number
							, Decoding_info.Dz_1_Number
							, Decoding_info.Dz_2_Number
							, Decoding_info.Dz_3_Number
							, Decoding_info.Dz_4_Number
							, Decoding_info.Dz_5_Number
							, Decoding_info.D1
							, Decoding_info.D2
							, Decoding_info.D3
							, Decoding_info.D4
							, Decoding_info.D5
							, (Decoding_info.CBC_Candidate / (int)Transmitted_Block)
							, (Decoding_info.TotalCounter / (int)Transmitted_Block)
							, (Decoding_info.CBC_STE / (int)Transmitted_Block / Message_Length)
						);*/

						double Execution_Time((clock() - start_time_for_SNR) / 1000.);

						// Store the real-time statistics of decoder 
						{
							ClearFile(Realtime_Data);
							WriteFile(Realtime_Data, snr_dB);
							WriteFile(Realtime_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
							WriteFile(Realtime_Data, BER);
							WriteFile(Realtime_Data, BLER);
							WriteFile(Realtime_Data, Avg_STE);
							WriteFile(Realtime_Data, Avg_Com);
							WriteFile(Realtime_Data, Avg_CandidateCodeWord);
							WriteFile(Realtime_Data, Avg_Binary_STE);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_STE);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_COM);
							WriteFile(Realtime_Data, Decoding_info.Worst_Case_Candidate);
							WriteFile(Realtime_Data, Decoding_info.Worst_OSC_Ratio);
							WriteFile(Realtime_Data, Decoding_info.MaxUsed_Stack);
							WriteFile(Realtime_Data, "-----------");
							//WriteFile(Realtime_Data, Avg_DM_STE);
							//WriteFile(Realtime_Data, Avg_DM_Com);
							//WriteFile(Realtime_Data, Avg_DN);
							WriteFile(Realtime_Data, ML_LB_BER);
							WriteFile(Realtime_Data, ML_LB_BLER);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, (Decoding_info.CBC_Candidate / Transmitted_Block));
							WriteFile(Realtime_Data, Decoding_info.TotalCounter / Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.CBC_STE / Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, (Decoding_info.DM_STE*Decoding_info.CBC_length) / Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, "\n\n");

						}
						/*
						{
							ClearFile(Realtime_Data);
							WriteFile(Realtime_Data, snr_dB);
							WriteFile(Realtime_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
							WriteFile(Realtime_Data, BER);
							WriteFile(Realtime_Data, BLER);
							WriteFile(Realtime_Data, Avg_STE);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.Dz_0_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_1_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_2_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_3_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_4_Number);
							WriteFile(Realtime_Data, Decoding_info.Dz_5_Number);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.D1);
							WriteFile(Realtime_Data, Decoding_info.D2);
							WriteFile(Realtime_Data, Decoding_info.D3);
							WriteFile(Realtime_Data, Decoding_info.D4);
							WriteFile(Realtime_Data, Decoding_info.D5);
							WriteFile(Realtime_Data, "-----------");
							WriteFile(Realtime_Data, Decoding_info.CBC_Candidate / (int)Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.TotalCounter / (int)Transmitted_Block);
							WriteFile(Realtime_Data, Decoding_info.CBC_STE / (int)Transmitted_Block / Message_Length);
							WriteFile(Realtime_Data, "\n\n");
						}*/
						if (ErrorBlock_Counter >= Operation_Parameter.ErrorBlock_Thr) break;
						//if (Transmitted_Block > 1000) break;
				}
			}
			
			/*
			string LLR_0 = "LLR_0.txt";
			string LLR_1 = "LLR_1.txt";
			vector<double> LLR_TemP = LLR_Distribution_1;
			reverse(LLR_TemP.begin(), LLR_TemP.end());
			LLR_Distribution_0 = LLR_TemP;
			for (int i = 0; i < 2000; ++i) {
				WriteFile(LLR_0, LLR_Distribution_0.at(i));
				WriteFile(LLR_1, LLR_Distribution_1.at(i));
			}
			cout << "END";*/

		}
		/*
		string Error_Acc = "Error_Acc.txt";
		for (int i = 0; i < 192; ++i) {
			WriteFile(Error_Acc, Decoding_info.Error_Accumulation.at(i));
		}
		cout << "End of simulation!" << endl;
		*/

		// for loop : Transmitted Block

		//
		//
		//
		
		// Store the complete statistics of the decoder 
		
		{
			WriteFile(Complete_Data, snr_dB);
			WriteFile(Complete_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
			WriteFile(Complete_Data, BER);
			WriteFile(Complete_Data, BLER);
			WriteFile(Complete_Data, Avg_STE);
			WriteFile(Complete_Data, Avg_Com);
			WriteFile(Complete_Data, Avg_CandidateCodeWord);
			WriteFile(Complete_Data, Avg_Binary_STE);
			WriteFile(Complete_Data, Decoding_info.Worst_Case_STE);
			WriteFile(Complete_Data, Decoding_info.Worst_Case_COM);
			WriteFile(Complete_Data, Decoding_info.Worst_Case_Candidate);
			WriteFile(Complete_Data, (double)Decoding_info.Worst_OSC_Ratio); //PoHan
			WriteFile(Complete_Data, Decoding_info.MaxUsed_Stack);
			WriteFile(Complete_Data, "-----------");
			WriteFile(Complete_Data, Avg_DM_STE);
			WriteFile(Complete_Data, Avg_DM_Com);
			WriteFile(Complete_Data, Avg_DN);
			WriteFile(Complete_Data, "-----------");
			WriteFile(Complete_Data, ML_LB_BER);
			WriteFile(Complete_Data, ML_LB_BLER);
			WriteFile(Complete_Data, (Decoding_info.CBC_Candidate / Transmitted_Block));
			WriteFile(Complete_Data, Decoding_info.TotalCounter / Transmitted_Block);
			WriteFile(Complete_Data, Decoding_info.CBC_STE / Transmitted_Block / Message_Length);
			WriteFile(Complete_Data, (Decoding_info.DM_STE*Decoding_info.CBC_length) / Transmitted_Block / Message_Length);
			WriteFile(Complete_Data, "-----------");
			/*
			//Pin test
			WriteFile("Hard_test.txt", snr_dB);
			WriteFile("Hard_test.txt", "-----------");
			for (int i = 0; i < 128; i++) {
				WriteFile("Hard_test.txt", Decoding_info.err_count[i]);
			}
			WriteFile("Hard_test.txt", "-----------");
			*/
			/*
			WriteFile(RSS_Data, snr_dB);
			WriteFile(RSS_Data, (double)Avg_Alpha);
			WriteFile(RSS_Data, (double)Decoding_info.Worst_OSC_Ratio);
			WriteFile(RSS_Data, "-----------");
			
			WriteFile(RSS_Data, snr_dB);
			WriteFile(RSS_Data, Decoding_info.Dz_0_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_1_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_2_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_3_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_4_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_5_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_6_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_7_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_8_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_9_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_10_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_11_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_12_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_13_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_14_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_15_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_16_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_17_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_18_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_19_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_20_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_21_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_22_Number);
			WriteFile(RSS_Data, Decoding_info.Dz_23_Number);
			WriteFile(RSS_Data, "-----------");
			*/
		}

		/*
		{
			WriteFile(Complete_Data, snr_dB);
			WriteFile(Complete_Data, snr_dB + 10 * log10(LinearBlockCode.Code_Rate));
			WriteFile(Complete_Data, BER);
			WriteFile(Complete_Data, BLER);
			WriteFile(Complete_Data, Avg_STE);
			WriteFile(Complete_Data, "-----------");
			WriteFile(Complete_Data, Decoding_info.Dz_0_Number);
			WriteFile(Complete_Data, Decoding_info.Dz_1_Number);
			WriteFile(Complete_Data, Decoding_info.Dz_2_Number);
			WriteFile(Complete_Data, Decoding_info.Dz_3_Number);
			WriteFile(Complete_Data, Decoding_info.Dz_4_Number);
			WriteFile(Complete_Data, Decoding_info.Dz_5_Number);
			WriteFile(Complete_Data, "-----------");
			WriteFile(Complete_Data, Decoding_info.D1);
			WriteFile(Complete_Data, Decoding_info.D2);
			WriteFile(Complete_Data, Decoding_info.D3);
			WriteFile(Complete_Data, Decoding_info.D4);
			WriteFile(Complete_Data, Decoding_info.D5);
			WriteFile(Realtime_Data, "-----------");
			WriteFile(Complete_Data, Decoding_info.CBC_Candidate / (int)Transmitted_Block);
			WriteFile(Complete_Data, Decoding_info.TotalCounter / (int)Transmitted_Block);
			WriteFile(Complete_Data, Decoding_info.CBC_STE / (int)Transmitted_Block / Message_Length);
			WriteFile(Complete_Data, "\n\n");
		}*/
		

		// Finish time for this simulation  
		double end_time_for_SNR(clock());
		Show_Current_Time();

		// Show the average spending time per block at this SNR 
		cout
			<< "\n*** Time Information for " << snr_dB << " [dB] : \n"
			<< "    Spending Time : " << ((end_time_for_SNR - start_time_for_SNR) / 1000.) << " [s] \n"
			<< "    Decoding Rate : " << (Transmitted_Block / ((end_time_for_SNR - start_time_for_SNR) / 1000.)) << " [ block / s ] \n\n\n";
		{
			WriteFile(Complete_Data, ((end_time_for_SNR - start_time_for_SNR) / 1000.));
			WriteFile(Complete_Data, (Transmitted_Block / ((end_time_for_SNR - start_time_for_SNR) / 1000.)));
			WriteFile(Complete_Data, "\n\n");
		}
	} // for end : SNR 

	// Show the total spending time of entire simulation procedure.
	cout << " Finish Time :";
	Show_Current_Time();
	cout
		<< "\n\n Total Execution Time = "
		<< (clock() - start_time) / 1000.
		<< " [s] " << endl;

	system("pause");
} // end main()



bool ML_LowerBound(
	vector<double> &received_signal, 
	vector<double> &transmitted_signal,
	vector<__int8> &estimated_codeword) {
	double
		Distance_TransCodeword = 0,
		Distance_EstCodeword = 0;
	vector<double> Estimated_Signal(received_signal.size(), 0);
	
	BPSK_Modulation(estimated_codeword, Estimated_Signal);

	Distance_TransCodeword = vectorDistance(transmitted_signal, received_signal);
	Distance_EstCodeword = vectorDistance(Estimated_Signal, received_signal);
	
	if (Distance_TransCodeword > Distance_EstCodeword) {
		return TRUE; //error
	}
	else {
		return FALSE;
	}
}

void Sort_function(
	vector<__int8>		&Rx_signal,
	vector<size_t>		&permutation_seq,
	vector<__int8>		&sorted_rx_signal)
{
	size_t counter(sorted_rx_signal.size());
	for (size_t i(0); i < counter; ++i)
		sorted_rx_signal.at(i) = Rx_signal.at(permutation_seq.at(i));
}