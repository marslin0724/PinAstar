#include "PinAstar/DefineClass.h"

void OPER_PARA::ShowDecoder()
{
		std::cout

		<< " *** Basic:"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << Hard_Decision_Decoder << ": Hard Decision Decoder "
		<< "\n " << A_Star_Original << ": A* Algorithm Original "
		<< "\n " << A_Star_Han << ": A* Algorithm, Han's, No function "	
		<< "\n " << A_Star_Parity << ": A* Algorithm + Parity "
		<< "\n " << A_Star_Parity_f << ": A* Algorithm + Parity-f "
	

		<< "\n\n *** PC:"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC << ": A* Algorithm PC "
		<< "\n " << A_Star_PC_out << ": A* Algorithm PC-out "
		<< "\n " << A_Star_2_Stack << ": A* Algorithm 2-Stack "
		<< "\n " << A_Star_3_Stack << ": A* Algorithm 3-Stack "
		<< "\n " << A_Star_4_Stack << ": A* Algorithm 4-Stack "


		<< "\n\n *** PC + Parity Heuristic Metric:"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_Parity_PC << ": A* Algorithm PC + Parity "
		<< "\n " << A_Star_Parity_PC_out << ": A* Algorithm PC-out + Parity "
		<< "\n " << A_Star_2_Parity_Stack << ": A* Algorithm 2-Stack + Parity "
		<< "\n " << A_Star_3_Parity_Stack << ": A* Algorithm 3-Stack + Parity "
		<< "\n " << A_Star_4_Parity_Stack << ": A* Algorithm 4-Stack + Parity "


		<< "\n\n *** DM-I:"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_out_DM_I << ": A* Algorithm PC-out DM-I "
		<< "\n " << A_Star_2_Stack_DM_I << ": A* Algorithm 2-Stack DM-I "
		<< "\n " << A_Star_3_Stack_DM_I << ": A* Algorithm 3-Stack DM-I "
		<< "\n " << A_Star_4_Stack_DM_I << ": A* Algorithm 4-Stack DM-I "


		<< "\n\n *** DM-I + Parity Heuristic Metric:"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_out_DM_I_Parity << ": A* PC-out Algorithm DM-I + Parity "
		<< "\n " << A_Star_2_Stack_DM_I_Parity << ": A* Algorithm 2-Stack DM-I + Parity "
		<< "\n " << A_Star_3_Stack_DM_I_Parity << ": A* Algorithm 3-Stack DM-I + Parity "
		<< "\n " << A_Star_4_Stack_DM_I_Parity << ": A* Algorithm 4-Stack DM-I + Parity "


		<< "\n\n *** DM-II"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_DM_II << ": A* Algorithm PC-DM-II "
		<< "\n " << A_Star_PCout_DM_II << ": A* Algorithm PCout-DM-II "
		<< "\n " << A_Star_2_Stack_DM_II << ": A* Algorithm 2-Stack-DM-II "
		<< "\n " << A_Star_3_Stack_DM_II << ": A* Algorithm 3-Stack-DM-II "
		<< "\n " << A_Star_4_Stack_DM_II << ": A* Algorithm 4-Stack-DM-II "
		<< "\n " << A_Star_BMA << ": A* Algorithm BMA "


		<< "\n\n *** OSC"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_OSC << ": A* Algorithm PC-OSC"
		<< "\n " << A_Star_PCout_OSC << ": A* Algorithm PC-out-OSC"
		<< "\n " << A_Star_2_Stack_OSC << ": A* Algorithm 2-Stack-OSC"
		<< "\n " << A_Star_3_Stack_OSC << ": A* Algorithm 3-Stack-OSC"
		<< "\n " << A_Star_4_Stack_OSC << ": A* Algorithm 4-Stack-OSC"


		<< "\n\n *** CBC"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PCout_CBC << ": A* Algorithm PCout-CBC "
		<< "\n " << A_Star_2_Stack_CBC << ": A* Algorithm 2-Stack-CBC "
		<< "\n " << A_Star_3_Stack_CBC << ": A* Algorithm 3-Stack-CBC "

		<< "\n\n *** Combinations"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PCout_CBC_OSC << ": A* Algorithm PCout-CBC-OSC "
		<< "\n " << A_Star_2_Stack_CBC_OSC << ": A* Algorithm 2-Stack-CBC-OSC "
		<< "\n " << A_Star_3_Stack_CBC_OSC << ": A* Algorithm 3-Stack-CBC-OSC "

		<< "\n\n *** Yin's Decoding Algorithm"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_out_Limited_CBC_OSC << ": A* Algorithm PCout-Limited CBC-OSC"
		<< "\n " << A_Star_PC_out_CBC_OSC_Dynamic_I << ": A* Algorithm PCout CBC-OSC Dynamic-I"
		<< "\n " << A_Star_PC_out_CBC_OSC_MultiDecoder << ": A* Algorithm PCout CBC-OSC DoubleDecoder"
		<< "\n " << A_Star_PC_out_CBC_OSC_Verified << ": A* Algorithm PCout CBC-OSC(lastest version of No.37)"
		<< "\n " << A_Star_PC_out_CBC_OSC_Adaptive_i << ": A* Algorithm PCout CBC-OSC Adaptive-i"

		<< "\n\n *** PoHan's Decoding Algorithm"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_OSC_II << ": A* Algorithm PC-OSCII"
		<< "\n " << A_Star_PC_out_Dynamic_CBC_OSC << ": A* Algorithm PC-out-Dynamic-CBC-OSC"
		<< "\n " << A_Star_PC_out_CBC_OSC_II << ": A* Algorithm PC-out-CBC-OSCII"
		<< "\n " << A_Star_PC_out_CBC_OSC_III << ": A* Algorithm PC-out-CBC-OSCIII (Sufficient condition)"
		<< "\n " << A_Star_RSS_test << ": A* Algorithm RSS Test"
		<< "\n " << A_Star_PC_out_CBC_OSC_WorstMetric << ": A* Algorithm PC-out-CBC-OSC-WorstMetric"
		<< "\n " << A_Star_2_Multi_Stack << ":  A* Algorithm 2-Multi-Stack"
		<< "\n " << A_Star_2_Multi_Stack_OSC << ":  A* Algorithm 2-Multi-Stack-OSC"
		<< "\n " << A_Star_PC_out_CBC_OSC_Adaptive_i_New << ": A* Algorithm PCout CBC-OSC Adaptive-i(lastest version of No.66)"
		<< "\n " << A_Star_Fano << ":  A* Algorithm with Fano Metric"
		<< "\n " << A_Star_PC_out_CBC_OSC_Fano << ":  A* Algorithm PCout CBC-OSC with Fano Metric"
		<< "\n " << A_Star_PC_out_CBC_OSC_Adaptive_i_Fano << ": A* Algorithm PCout CBC-OSC Adaptive-i with Fano Metric"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC << ":  A* Algorithm 2-Base-PCout-CBC-OSC"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Latest << ":  A* Algorithm 2-Base-PCout-CBC-OSC(lastest version of No.101)"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Parallel << ":  A* Algorithm 2-Base-PCout-CBC-OSC-Parallel"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i Parallel"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i Fano Sufficient Condition"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i Parallel Fano Sufficient Condition"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i Fano"
		<< "\n " << A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano << ":  A* Algorithm 2-Base-PCout-CBC-OSC Adaptive-i Parallel Fano"

		<< "\n\n *** Pin's Decoding Algorithm"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_Segment_Orignal << ": A* Algorithm using Segment"
		<< "\n " << A_Star_Segment_Orignal_ver2 << ": A* Algorithm using Segment ver2"
		<< "\n " << Hard_test << ": Hard_desicion_test"
		<< "\n " << MinD_test << ": Minimum distance test"
		<< "\n " << A_Star_Section_PC << ": A_Star_Section_PC"
		<< "\n " << A_Star_Section_PC_out << ": A_Star_Section_PC_out"
		<< "\n " << SPA_A_Star << ": SPA decoder + A*"
		<< "\n " << SPA_A_Star_ver2 << ": SPA decoder + A* version2"
		<< "\n " << A_Star_N_Section << ": A* N Section decoder"
		<< "\n " << A_Star_5_Stack << ": A* 5 Stack"
		<< "\n " << A_Star_6_Stack << ": A* 6 Stack"
		<< "\n " << A_Star_7_Stack << ": A* 7 Stack"
		<< "\n " << SPA_AvgLLR        << ": SPA Average LLR"
		<< "\n " << A_Star_Rep        << ": A* repeation"
		<< "\n " << A_Star_Section_MultiStack << ": A* Section MultiStack"
		<< "\n " << SPA_AStar << ": SPA + A* using concatenated"
		<< "\n " << Brute_Force_SISO << ": Brute_Force Soft-in/Soft-out"
		<< "\n " << Brute_Force_AStar << ": Brute_Force + A* "
		<< "\n " << Brute_Force_AStar_ver2 << ": Brute_Force + A* version2"
		<< "\n " << BF_AStar_interleaver_ver2 << ": Brute Force + A* + interleaver version2"
		<< "\n " << BF_SPA_AStar << ":Brute Force + SPA + A* "
		<< "\n " << ABP_Decoder << ": ABP decoder"
		<< "\n " << BF_ABP_AStar << ": Brute Force + ABP(BGA) + A* decoder"
		<< "\n " << ABP_AStar	 << ": ABP + A* decoder"
		<< "\n " << BF_AStar_interleaver << ": Brute_Force + A* + interleaver"
		<< "\n " << AStar_Adaptive_i	 << ": AStar_Adaptive_i"
		<< "\n " << RS_Hard_Desicion	 << ": RS Hard Desicion"
		<< "\n " << MRIP_BF_AStar		<< ": MRIP BF AStar"
		<< "\n " << Subcode_MRIP_AStar	<< ": Using result of SISO+Astar as inital path in all Astar"
		<< "\n " << AStar_MultipleStack << ": AStar MultipleStack"
		<< "\n " << AStar_constraint_i_noComparison << ": AStar_constraint_i_noComparison"
		<< "\n " << Subcode_AStar_AStar				<< ": Subcode AStar for AStar using"
		<< "\n " << BF_SPA_BF_AStar					<< ":BF + SPA + BF + AStar"
		<< "\n " << ABP_MRIP_AStar		<< ":using ABP as the initial Path A* decoder"
		<< "\n " << Subcode24_12_MRIP_AStar << ":Subcode24_12_MRIP_AStar"
		<< "\n " << SISO_Path_AStar		<< ":SISO_Path_AStar-CBC-OSC"
		<< "\n " << SISO_Path_AStar_CBC2_OSC <<":SISO_Path_AStar_CBC2_OSC"
		<< "\n " << AStar_constraint_i_MinHeap<< ": AStar_constraint_i_MinHeap"
		<< "\n " << SISO_Path_AStar_MinHeap	<< ": SISO_Path_AStar_MinHeap"
		<< "\n " << SISO_Path_AStar_Adaptive<< ": SISO_Path_AStar_Adaptive(first decoder-i = 2)"
		<< "\n " << SISO_Path_AStar_MultiBase_Adaptive << ": SISO_Path_AStar_MultiBase_Adaptive"
		<< "\n " << AStar_MultipleStack_SC	<< ": AStar (i-1)-Stack_Sufficient Condition"
		<< "\n " << AStar_section_PC_MinHeap<< ": AStar_section_PC_MinHeap"
		<< "\n " << SISO_Path_AStar_Adaptive_MinHeap <<": SISO_Path_AStar_Adaptive_MinHeap"

		<<"\n\n *** Other Decoding Algorithm"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << Polar_Code_BP_Decoder << ": BP Decoding Algorithm"
		<< "\n " << LDPC_Sum_Product_Decoder << ": SPA Algorithm"
		<< "\n " << Majority_Soft_Decoding << ": Majority Soft Decoding Algorithm"
		<< "\n " << Majority_Hard_Decoding << ": Majority Hard Decoding Algorithm"
		<< "\n " << Majority_Soft_Decoding_ver2 << ": Majority Soft Decoding Algorithm ver2"
		<< "\n " << LDPC_VC_RBP_Algorithm << ": LDPC decoding method: Residual BP based on V2C"
		

		<< "\n\n *** Test"
		<< "\n ---------------------------------------------------------------"
		<< "\n " << A_Star_PC_out_CBC_OSC_Test      << ": A* Algorithm PCout-CBC-OSC-Test "
		<< "\n " << A_Star_PC_out_CBC_OSC_Test2     << ": A* Algorithm PCout-CBC-OSC-Test2 "
		

		<< "\n\n Select a decoding algorithm : ";
}

void OPER_PARA::DecodingParameterKeyIn() {

	/*if (
		(Decoder_Version == A_Star_Parity_f)
		) {
		cout << " Alpha = ";
		cin >> Alpha;
	}
	else if (
		(Decoder_Version == A_Star_PC) ||
		(Decoder_Version == A_Star_PC_out) ||
		(Decoder_Version == A_Star_Parity_PC) ||
		(Decoder_Version == A_Star_Parity_PC_out)
		) {
		cout << " i = ";
		cin >> Constraint_i;
	}
	else if (
		(Decoder_Version == A_Star_PC_out_DM_I) ||
		(Decoder_Version == A_Star_PC_out_DM_I_Parity)
		) {
		cout << " i = ";
		cin >> Constraint_i;

		cout << " j = ";
		cin >> Constraint_j;

		cout << " Control Level = ";
		cin >> Control_Level;
	}
	else if (
		(Decoder_Version == A_Star_PC_DM_II) ||
		(Decoder_Version == A_Star_PCout_DM_II)
		) {
		cout << " i = ";
		cin >> Constraint_i;

		cout << " j = ";
		cin >> Constraint_j;

		cout << " Check Level = ";
		cin >> Check_Level;

		cout << " Control Level = ";
		cin >> Control_Level;

		cout << " DM Stack Size = ";
		cin >> DM_StackSize;
	}
	else if (
		(Decoder_Version == A_Star_2_Stack_DM_II) ||
		(Decoder_Version == A_Star_3_Stack_DM_II) ||
		(Decoder_Version == A_Star_4_Stack_DM_II)
		) {
		cout << " j = ";
		cin >> Constraint_j;

		cout << " Check Level = ";
		cin >> Check_Level;

		cout << " Control Level = ";
		cin >> Control_Level;

		cout << " DM Stack Size = ";
		cin >> DM_StackSize;
	}
	else if (
		(Decoder_Version == A_Star_2_Stack_DM_I) ||
		(Decoder_Version == A_Star_3_Stack_DM_I) ||
		(Decoder_Version == A_Star_4_Stack_DM_I) ||
		(Decoder_Version == A_Star_2_Stack_DM_I_Parity) ||
		(Decoder_Version == A_Star_3_Stack_DM_I_Parity) ||
		(Decoder_Version == A_Star_4_Stack_DM_I_Parity)
		) {
		cout << " j = ";
		cin >> Constraint_j;

		cout << " Control Level = ";
		cin >> Control_Level;
	}
	else if (
		(Decoder_Version == A_Star_BMA)
		) {
		cout << " j = ";
		cin >> Constraint_j;

		cout << " Check Level = ";
		cin >> Check_Level;

		cout << " Control Level = ";
		cin >> Control_Level;

		cout << " DM Stack Size = ";
		cin >> DM_StackSize;
	}
	else if (
		(Decoder_Version == A_Star_PCout_CBC)
		) {
	}
	else if(
		(Decoder_Version == A_Star_2_Stack_CBC) ||
		(Decoder_Version == A_Star_3_Stack_CBC)
		){

	}*/

	//
	if (Decoder_Version == A_Star_Parity_f) {
		// Parity_f_Flag = TRUE;
		cout << " Alpha ( 0 < x < 1 , R ) = ";
		cin >> Alpha;
	}
	// RSS
	if (
		(Decoder_Version == A_Star_PC) ||
		(Decoder_Version == A_Star_PC_out) ||
		(Decoder_Version == A_Star_Parity_PC) ||
		(Decoder_Version == A_Star_Parity_PC_out) ||
		(Decoder_Version == A_Star_PC_out_DM_I) ||
		(Decoder_Version == A_Star_PC_out_DM_I_Parity) ||
		(Decoder_Version == A_Star_PC_DM_II) ||
		(Decoder_Version == A_Star_PCout_DM_II) ||
		(Decoder_Version == A_Star_PC_OSC) ||
		(Decoder_Version == A_Star_PCout_OSC) ||
		(Decoder_Version == A_Star_PCout_CBC) ||
		(Decoder_Version == A_Star_PCout_CBC_OSC) || 
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Verified) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_PC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_III) ||
		(Decoder_Version == A_Star_PC_out_Dynamic_CBC_OSC) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_WorstMetric) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_New) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Fano)	||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Latest) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano) ||
		(Decoder_Version == A_Star_Section_PC) ||
		(Decoder_Version == A_Star_Section_PC_out||
		 Decoder_Version == AStar_Adaptive_i ||
			Decoder_Version == AStar_MultipleStack ||
			Decoder_Version == AStar_constraint_i_noComparison ||
			Decoder_Version == Subcode_AStar_AStar ||
			Decoder_Version == ABP_MRIP_AStar ||
			Decoder_Version == Subcode24_12_MRIP_AStar ||
			Decoder_Version == SISO_Path_AStar ||
			Decoder_Version == SISO_Path_AStar_CBC2_OSC ||
			Decoder_Version == AStar_constraint_i_MinHeap ||
			Decoder_Version == SISO_Path_AStar_MinHeap ||
			Decoder_Version == SISO_Path_AStar_Adaptive ||
			Decoder_Version == SISO_Path_AStar_MultiBase_Adaptive ||
			Decoder_Version == AStar_MultipleStack_SC||
			Decoder_Version == AStar_section_PC_MinHeap ||
			Decoder_Version == SISO_Path_AStar_Adaptive_MinHeap)
		) {
		cout << " i = ";
		cin >> Constraint_i;
	}
	
	// DM-I
	if (
		(Decoder_Version == A_Star_PC_out_DM_I) ||
		(Decoder_Version == A_Star_2_Stack_DM_I) ||
		(Decoder_Version == A_Star_3_Stack_DM_I) ||
		(Decoder_Version == A_Star_4_Stack_DM_I) ||
		(Decoder_Version == A_Star_PC_out_DM_I_Parity) ||
		(Decoder_Version == A_Star_2_Stack_DM_I_Parity) ||
		(Decoder_Version == A_Star_3_Stack_DM_I_Parity) ||
		(Decoder_Version == A_Star_4_Stack_DM_I_Parity) ||
		(Decoder_Version == A_Star_PC_DM_II) ||
		(Decoder_Version == A_Star_PCout_DM_II) ||
		(Decoder_Version == A_Star_2_Stack_DM_II) ||
		(Decoder_Version == A_Star_3_Stack_DM_II) ||
		(Decoder_Version == A_Star_4_Stack_DM_II)||
		Decoder_Version == SISO_Path_AStar ||
		Decoder_Version == SISO_Path_AStar_CBC2_OSC ||
		Decoder_Version == SISO_Path_AStar_MinHeap ||
		Decoder_Version == SISO_Path_AStar_Adaptive ||
		Decoder_Version == SISO_Path_AStar_MultiBase_Adaptive ||
		Decoder_Version == SISO_Path_AStar_Adaptive_MinHeap
		) {
		cout << " j = ";
		cin >> Constraint_j;

		cout << " Control Level ( k < x < n , N ) = ";
		cin >> Control_Level;
	}
	// DM-II
	if (
		(Decoder_Version == A_Star_PC_DM_II) ||
		(Decoder_Version == A_Star_PCout_DM_II) ||
		(Decoder_Version == A_Star_2_Stack_DM_II) ||
		(Decoder_Version == A_Star_3_Stack_DM_II) ||
		(Decoder_Version == A_Star_4_Stack_DM_II) ||
		(Decoder_Version == A_Star_BMA) 
		) {
		cout << " Check Level ( 0 < x < k , N ) = ";
		cin >> Check_Level;

		cout << " DM Stack Size = ";
		cin >> DM_StackSize;
	}
	// OSC
	if (
		(Decoder_Version == A_Star_PC_OSC) ||
		(Decoder_Version == A_Star_PCout_OSC) ||
		(Decoder_Version == A_Star_2_Stack_OSC) ||
		(Decoder_Version == A_Star_3_Stack_OSC) ||
		(Decoder_Version == A_Star_4_Stack_OSC) ||
		(Decoder_Version == A_Star_PCout_CBC_OSC) ||
		(Decoder_Version == A_Star_2_Stack_CBC_OSC) ||
		(Decoder_Version == A_Star_3_Stack_CBC_OSC) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Verified)  ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_PC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_Dynamic_CBC_OSC) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_III) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_WorstMetric) ||
		(Decoder_Version == A_Star_2_Multi_Stack_OSC) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_New) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Fano)	||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Latest) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano||
		 Decoder_Version == AStar_Adaptive_i ||
		 Decoder_Version == Brute_Force_AStar||
		 Decoder_Version == SISO_Path_AStar ||
		 Decoder_Version == SISO_Path_AStar_CBC2_OSC ||
		 Decoder_Version == SISO_Path_AStar_MinHeap ||
		 Decoder_Version == SISO_Path_AStar_Adaptive ||
		 Decoder_Version == SISO_Path_AStar_MultiBase_Adaptive ||
		 Decoder_Version == AStar_constraint_i_MinHeap ||
		 Decoder_Version == AStar_constraint_i_noComparison || 
		 Decoder_Version == SISO_Path_AStar_Adaptive_MinHeap)
		) {
		cout << " OSC-Alpha ( 0 < x < 1 , R ) = ";
		cin >> OSC_Alpha;
	}
	// CBC
	if (
		(Decoder_Version == A_Star_PCout_CBC) ||
		(Decoder_Version == A_Star_2_Stack_CBC) ||
		(Decoder_Version == A_Star_3_Stack_CBC) ||
		(Decoder_Version == A_Star_PCout_CBC_OSC) ||
		(Decoder_Version == A_Star_2_Stack_CBC_OSC) ||
		(Decoder_Version == A_Star_3_Stack_CBC_OSC) || 
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Verified) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_WorstMetric) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_New) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Fano)	||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Latest) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_Dynamic_CBC_OSC) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_III ||
		 Decoder_Version == Brute_Force_AStar)
		) {

		cout << " j = ";
		cin >> Constraint_j;

		cout << " Control Level ( k < x < n , N ) = ";
		cin >> Control_Level;

		cout << " CBC-Flipping-Bit (1 or 2 bits) = ";
		cin >> CBC_FlippingBit;
	}
	//OSCII
	if (
		(Decoder_Version == A_Star_PC_OSC_II) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_II) 
		) {
		cout << " New-OSC-Alpha ( 0 < x < 1 , R ) = ";
		cin >> New_OSC_Alpha;
	}
	//Worst Metric
	if (
		(Decoder_Version == A_Star_PC_out_CBC_OSC_WorstMetric)
		) {
		cout << " Worst Metric Ratio = ";
		cin >> Worst_metric_ratio;
	}
	//Fano Metric
	if (
		(Decoder_Version == A_Star_Fano) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Fano) ||
		(Decoder_Version == A_Star_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano) ||
		(Decoder_Version == A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano)
		) {
		cout << " Fano Metric Parameter = ";
		cin >> Fano_Metric_Parameter;
	}

	//Multi-Stack
	if (
		(Decoder_Version == A_Star_2_Multi_Stack) ||
		(Decoder_Version == A_Star_2_Multi_Stack_OSC)
		) {
		cout << " Maximum multi-stack size : ";
		cin >> Multi_Stack_Size;
	}
	//section
	if (Decoder_Version == A_Star_N_Section) {
		cout << "N_Section : ";
		cin >> N_section;
	}
	if ((Decoder_Version == A_Star_Section_PC) ||
		(Decoder_Version == A_Star_Section_PC_out) ||
		(Decoder_Version == A_Star_Section_MultiStack ||
		 Decoder_Version == AStar_section_PC_MinHeap)
		) {
		cout << " Section1_i : ";
		cin >> section1_i;
		cout << " Section2_i : ";
		cin >> section2_i;
	}
	if (Decoder_Version == SPA_A_Star ||
		Decoder_Version == SPA_A_Star_ver2||
		Decoder_Version == ABP_Decoder||
		Decoder_Version == BF_ABP_AStar ||
		Decoder_Version == ABP_AStar ||
		Decoder_Version == BF_SPA_BF_AStar ||
		Decoder_Version == ABP_MRIP_AStar ||
		Decoder_Version == Subcode24_12_MRIP_AStar) {
		cout << "SPA_iteration :";
		cin >> SPA_I;
	}

	
	//<< "\n " << Polar_Code_BP_Decoder << ": BP Decoding Algorithm"
	//<< "\n " << LDPC_Sum_Product_Decoder << ": SPA Algorithm"
	if (Decoder_Version == Polar_Code_BP_Decoder) cout << "BP-Decoding: Iteration = " << BP_Decode_Iter;
	if (Decoder_Version == LDPC_Sum_Product_Decoder) cout << "SPA-Decoding: Iteration = " << SPA_Iter;
}

void OPER_PARA::DecodingParameterShow() {
	system("cls"); // clear screen
	switch (Decoder_Version)
	{
	case A_Star_Original:
		cout << " \n A* I (Original)";
		break;
	case A_Star_Parity:
		cout << " \n A* Parity";
		break;
	case A_Star_Parity_f:
		cout << " \n A* Parity-f, Alpha = " << Alpha;
		break;

		// RSS
	case A_Star_PC:
		cout
			<< " \n A* PC-"
			<< Constraint_i << " ";
		break;
	case A_Star_PC_out:
		cout
			<< " \n A* PC-out-"
			<< Constraint_i;
		break;
	case A_Star_2_Stack:
		cout << " \n A* 2-Stack";
		break;
	case A_Star_3_Stack:
		cout << " \n A* 3-Stack";
		break;
	case A_Star_4_Stack:
		cout << " \n A* 4-Stack";
		break;

		// RSS + Parity Heuristic Function
	case A_Star_Parity_PC:
		cout
			<< " \n A* Parity + PC-"
			<< Constraint_i << " ";
		break;
	case A_Star_Parity_PC_out:
		cout
			<< " \n A* Parity + PC-out-"
			<< Constraint_i << " ";
		break;
	case A_Star_2_Parity_Stack:
		cout << " \n A* Parity + 2-Stack";
		break;
	case A_Star_3_Parity_Stack:
		cout << " \n A* Parity + 3-Stack";
		break;
	case A_Star_4_Parity_Stack:
		cout << " \n A* Parity + 4-Stack";
		break;

		// DM-I 
	case A_Star_PC_out_DM_I:
		cout
			<< " \n A* PC-out-" << Constraint_i
			<< "-" << Constraint_j << " DM-I "
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Stack_DM_I:
		cout
			<< " \n A* 2-Stack-" << Constraint_j << " DM-I "
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_3_Stack_DM_I:
		cout
			<< " \n A* 3-Stack-" << Constraint_j << " DM-I "
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_4_Stack_DM_I:
		cout
			<< " \n A* 4-Stack-" << Constraint_j << " DM-I "
			<< ", Control Level = " << Control_Level;
		break;

		// DM-I + Parity Heuristic Function 
	case A_Star_PC_out_DM_I_Parity:
		cout
			<< " \n A* PC-out-" << Constraint_i
			<< "-" << Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Stack_DM_I_Parity:
		cout
			<< " \n A* 2-Stack-" << Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_3_Stack_DM_I_Parity:
		cout
			<< " \n A* 3-Stack-" << Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_4_Stack_DM_I_Parity:
		cout
			<< " \n A* 4-Stack-" << Constraint_j
			<< " DM-I + Parity"
			<< ", Control Level = " << Control_Level;
		break;

		// DM-II
	case A_Star_PC_DM_II:
		cout
			<< " \n  A* PC-DM-II-"
			<< Constraint_i << "-" << Constraint_j << " "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
	case A_Star_PCout_DM_II:
		cout
			<< " \n A* PC-out-DM-II-"
			<< Constraint_i << "-" << Constraint_j << " "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
	case A_Star_2_Stack_DM_II:
		cout
			<< " \n A* 2-Stack-"
			<< Constraint_j << "-DM-II "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
	case A_Star_3_Stack_DM_II:
		cout
			<< " \n A* 3-Stack-"
			<< Constraint_j << "-DM-II "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
	case A_Star_4_Stack_DM_II:
		cout
			<< " \n A* 4-Stack-"
			<< Constraint_j << "-DM-II "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
	case A_Star_BMA:
		cout
			<< " \n A* BMA-"
			<< Constraint_j << " "
			<< "\n Check Level = " << Check_Level
			<< "\n Control Level = " << Control_Level
			<< "\n DM Stack Size = " << DM_StackSize;
		break;
		
		// OSC
	case A_Star_PC_OSC:
		cout
			<< "\n A* PC-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha;
		break;
	case A_Star_PCout_OSC:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha;
		break;
	case A_Star_2_Stack_OSC:
		cout
			<< "\n A* 2-Stack"
			<< "\n OSC-" << OSC_Alpha;
		break;
	case A_Star_3_Stack_OSC:
		cout
			<< "\n A* 3-Stack"
			<< "\n OSC-" << OSC_Alpha;
		break;
	case A_Star_4_Stack_OSC:
		cout
			<< "\n A* 4-Stack"
			<< "\n OSC-" << OSC_Alpha;
		break;

		// CBC
	case A_Star_PCout_CBC:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Stack_CBC:
		cout
			<< "\n A* 2-Stack"
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_3_Stack_CBC:
		cout
			<< "\n A* 3-Stack"
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;

		// CBC + OSC
	case A_Star_PCout_CBC_OSC:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Stack_CBC_OSC:
		cout
			<< "\n A* 2-Stack" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_3_Stack_CBC_OSC:
		cout
			<< "\n A* 3-Stack" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
		// Chou-Yin
	case A_Star_PC_out_CBC_OSC_MultiDecoder:
		cout
			<< "\n A* Pcout-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< "\n (DoubleDecoder)";
		break;
	case A_Star_PC_out_CBC_OSC_Verified:
		cout
			<< "\n A* PC-out" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< "\n (latest version)";
		break;
	case A_Star_PC_out_CBC_OSC_Adaptive_i:
		cout
			<< "\n A* PC-out" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< "\n Adaptive-i";
		break;
		// OSCII
	case A_Star_PC_OSC_II:
		cout 
			<< "\n A* PC-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n New OSC-" << New_OSC_Alpha;
		break;
	case A_Star_PC_out_Dynamic_CBC_OSC:
		cout
			<< "\n A* PC-out-Dynamic" << Constraint_i_ratio
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit;
		break;
	case A_Star_PC_out_CBC_OSC_II:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n OSC-" << OSC_Alpha
			<< "\n New OSC-" << New_OSC_Alpha;
		break;
	case A_Star_PC_out_CBC_OSC_III:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
		// Worst Metric Ratio
	case A_Star_PC_out_CBC_OSC_WorstMetric:
		cout
			<< "\n A* PC-" << Constraint_i
			<< "\n CBC-" << CBC_FlippingBit
			<< "\n OSC-" << OSC_Alpha
			<< "\n Worst Metric Ratio-" << Worst_metric_ratio;
		break;
		// Multi-Stack
	case A_Star_2_Multi_Stack:
		cout
			<< " \n A* 2-Multi-Stack"
			<< " \n Max size - " << Multi_Stack_Size;
		break;
		// Multi-Base
	case A_Star_2_Base_PC_out_CBC_OSC:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Latest:
		cout 
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Parallel:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base";
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base"
			<< "\n Sufficient Condition";
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base"
			<< "\n Fano Parameter = " << Fano_Metric_Parameter;
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base"
			<< "\n Parallel";
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base"
			<< "\n 3 phase"
			<< "\n Parallel";
		break;
	case A_Star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano:
		cout
			<< "\n A* PC-out-" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Adaptive "
			<< "\n 2 Base"
			<< "\n Parallel"
			<< "\n Fano Parameter = " << Fano_Metric_Parameter;;
		break;
		// Fano-metric
	case A_Star_Fano:
		cout 
			<< " \n A* Original with Fano Metric";
		break;
	case A_Star_PC_out_CBC_OSC_Fano:
		cout
			<< "\n A* PC-out" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Fano Parameter = " << Fano_Metric_Parameter;
		break;
	case A_Star_PC_out_CBC_OSC_Adaptive_i_Fano:
		cout
			<< "\n A* PC-out" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< ", Control Level = " << Control_Level
			<< "\n Fano Parameter = " << Fano_Metric_Parameter
			<< "\n Adaptive-i";
		break;
		// Adaptive-i
	case A_Star_PC_out_CBC_OSC_Adaptive_i_New:
		cout
			<< "\n A* PC-out" << Constraint_i
			<< "\n OSC-" << OSC_Alpha
			<< "\n CBC-" << CBC_FlippingBit
			<< "\n Adaptive-i";
		break;
		// Other Decoding Algorithm
	case Polar_Code_BP_Decoder:
		cout
			<< "\n BP decoding Algorithm: Iteration = " << BP_Decode_Iter;
		break;
	case  LDPC_Sum_Product_Decoder:
		cout
			<< "\n Sum Product Decoding Algorithm: Iteration = " << SPA_Iter;
		break;
	case A_Star_Segment_Orignal:
		cout
			<< "\n Section Original decoding Algorithm:  " ;
		break;
	case A_Star_Section_PC:
		cout
			<< "\n A* Section PC-" << section1_i << "-" << section2_i
			<< "\n section1_i = " << section1_i
			<< "\n section2_i = " << section2_i;
		break;
	case A_Star_Section_PC_out:
		cout
			<< "\n A* Section PC-out-" << section1_i << "-" << section2_i
			<< "\n section1_i = " << section1_i
			<< "\n section2_i = " << section2_i;
		break;
	case SPA_A_Star:
		cout << "\n SPA_A*:";
		break;
	case SPA_A_Star_ver2:
		cout << "\n SPA_A* version2" << ",SPA_Iteration = " << SPA_I;
		break;
	case A_Star_N_Section:
		cout << "\n A* " << N_section << "-Section ";
		break;
	case A_Star_5_Stack:
		cout << "\n A* 5 Stack ";
		break;
	case A_Star_6_Stack:
		cout << "\n A* 6 Stack ";
		break;
	case A_Star_7_Stack:
		cout << "\n A* 7 Stack ";
		break;
	case A_Star_Rep:
		cout << "\n A* repeatition";
		break;
	case A_Star_Section_MultiStack:
		cout 
			<< "\n A* Section MultiStack" << section1_i << "-" << section2_i
			<< "\n section1_i = " << section1_i
			<< "\n section2_i = " << section2_i;
		break;
	case SPA_AStar:
		cout << "\n SPA + A* using concatenated code";
		break;
	case Brute_Force_SISO:
		cout << "\n Brute Force SISO";
		break;
	case Brute_Force_AStar:
		cout << "\n Brute Force + A*";
		break;
	case Brute_Force_AStar_ver2:
		cout << "\n Brute Force + A* version2";
		break;
	case BF_AStar_interleaver_ver2:
		cout << "\n BF_AStar_interleaver_ver2";
		break;
	case ABP_Decoder:
		cout << "\n ABP_Decoder";
		break;
	case BF_ABP_AStar:
		cout << "\n SISO + ABP + AStar";
		break;
	case ABP_AStar:
		cout << "\n ABP + A*";
		break;
	case BF_AStar_interleaver:
		cout << "BF_AStar_interleaver";
		break;
	case Subcode_MRIP_AStar:
		cout << "Subcode_MRIP_AStar";
		break;
	case AStar_MultipleStack:
		cout << "AStar_MultipleStack";
		break;
	case AStar_constraint_i_noComparison:
		cout << "AStar_constraint_i_noComparison";
		break;
	case Subcode_AStar_AStar:
		cout << " Subcode AStar for AStar using";
		break;
	case BF_SPA_BF_AStar:
		cout << "BF+SPA+BF+AStar";
		break;
	case ABP_MRIP_AStar:
		cout << "ABP_MRIP_AStar";
		break;
	case Subcode24_12_MRIP_AStar:
		cout << "Subcode24_12_MRIP_AStar";
		break;
	case SISO_Path_AStar:
		cout << "SISO_Path_AStar,i =" << Constraint_i
			<< ",OSC-" << OSC_Alpha
			<< ",Control Level = " << Control_Level;
		break;
	case SISO_Path_AStar_CBC2_OSC:
		cout << "SISO_Path_AStar,i' =" << Constraint_i+2
			<< ",OSC-" << OSC_Alpha
			<< ",Control Level = " << Control_Level;
		break;
	case SISO_Path_AStar_MinHeap:
		cout << "SISO_Path_AStar,i =" << Constraint_i
			<< ",OSC-" << OSC_Alpha
			<< ",Control Level = " << Control_Level;
		break;
	case SISO_Path_AStar_Adaptive:
		cout << "SISO_Path_AStar_Adaptive,i1=2,i =" << Constraint_i
			<< ",OSC-" << OSC_Alpha
			<< ",Control Level = " << Control_Level;
		break;
	case SISO_Path_AStar_MultiBase_Adaptive:
		cout << "SISO_Path_AStar_MultiBase_Adaptive,i1=2,i =" << Constraint_i
			<< ",OSC-" << OSC_Alpha
			<< ",Control Level = " << Control_Level;
		break;
	case AStar_MultipleStack_SC:
		cout << "AStar_MultipleStack_SC,i = " << Constraint_i << "," << Constraint_i - 1 << "-Stack ";
		break;
	case AStar_section_PC_MinHeap:
		cout << "AStar_section_PC_MinHeap,i = " << Constraint_i << ",i1 = "<<section1_i<<",i2 = "<<section2_i;
		break;
	case SISO_Path_AStar_Adaptive_MinHeap:
		cout << "SISO_Path_AStar_Adaptive_MinHeap,i = " << Constraint_i << ", OSC-" << OSC_Alpha << ",";
	}
}

OPER_PARA::OPER_PARA(CODE Code)
{
	using namespace std;
	bool ReKey_Flag = TRUE;
	while (ReKey_Flag)
	{
		system("cls");
		cout
			<< "\n ----------------------------------------------------- \n"
			<< " A* Decoding Simulation for " << Code.Title
			<< endl;

		// Select Decoder
		ShowDecoder();
		cin >> Decoder_Version;
		
		// Select the way for which input the decoding data
		cout << " Read decoding parameters from 'SimulatedParameters.txt' ? (No - 0 / Yes - 1) = "; 
		cin >> Parameter_ReadTxt_Flag;

		if (!Parameter_ReadTxt_Flag) {
			DecodingParameterKeyIn();
			cout << " Block Number = ";					cin >> block_number;
			cout << " Start SNR [dB] = ";				cin >> SNR_dB_start;
			cout
				<< " ?? Es/N0 [dB] = "
				<< SNR_dB_start + 10. * log10(Code.Code_Rate)
				<< " [dB]\n" << endl;
			cout 
				<< " End SNR [dB] = ";					cin >> SNR_dB_end;
			cout
				<< " ?? Es/N0 [dB] = "
				<< SNR_dB_end + 10. * log10(Code.Code_Rate)
				<< " [dB]\n" << endl;
			cout
				<< " Step SNR) [dB] = ";				cin >> SNR_dB_step;
			cout 
				<< " Block Step = ";					cin >> BlockStep;
			cout 
				<< " Error Block Threshold = ";			cin >> ErrorBlock_Thr;
			cout
				<< " Stack Size = ";					cin >> StackSize;
			cout
				<< "\n -----------------------------------------------------\n"
				<< " Do you want to change any decoding parameters ? (No - 0 / Yes - 1) = ";

			cin >> ReKey_Flag;
		}
		// Read decoding parameters from .txt
		else { 
			string filename("SimulatedParameters.txt");
			fstream fp;
			fp.open(filename, ios::out | ios::in);
			if (!fp) cout << "Fail to open file: " << filename << endl;

			vector<double> ReadParameter(15, 0);
			for (size_t i(0); i < 15; ++i) {
				fp >> ReadParameter.at(i);
			}
			Eb_N0_1_or_Es_N0_2      = ReadParameter[0];
			block_number			= (__int64) ReadParameter[1];
			SNR_dB_start			= ReadParameter[2];
			SNR_dB_end				= ReadParameter[3];
			SNR_dB_step				= ReadParameter[4];
			StackSize				= (size_t) ReadParameter[5];
			BlockStep				= (size_t) ReadParameter[6];
			ErrorBlock_Thr			= (size_t) ReadParameter[7];
			Constraint_i			= (size_t) ReadParameter[8];
			Constraint_j			= (size_t) ReadParameter[9];
			Control_Level			= (size_t) ReadParameter[10];
			OSC_Alpha				= ReadParameter[11];
			CBC_FlippingBit			= (size_t) ReadParameter[12];

			cout 
				<< "\n  1. Block Number = "; Present_Large_Number(block_number);
			cout
				<< "\n  2. Start - SNR [dB] = "			<< SNR_dB_start
				<< "\n  3. End - SNR [dB] = "			<< SNR_dB_end
				<< "\n  4. Step - SNR [dB] = "			<< SNR_dB_step
				<< "\n  5. Stack Size = ";	Present_Large_Number(StackSize);
			cout 
				<< "\n  6. Block Step = "				<< BlockStep
				<< "\n  7. Error Block Threshold = "	<< ErrorBlock_Thr
				<< "\n  8. Constraint i = "				<< Constraint_i
				<< "\n  9. Constraint j = "				<< Constraint_j
				<< "\n 10. Control Level = "			<< Control_Level
				<< "\n 11. OSC-Alpha = "				<< OSC_Alpha
				<< "\n 12. CBC-x = "					<< CBC_FlippingBit
				<< "\n -----------------------------------------------------\n"
				<< " Do you want to change any parameters ? (No - 0 / Yes - 1) = "
				<< endl;
			cin >> ReKey_Flag;
		}
	}

	// Show the operating parameters of decoder
	DecodingParameterShow();
	
	// Show the parameters of this simulation
	cout
		<< " Simulation of Soft-Decision Decoding Algorithm for " << Code.Title << "\n"
		<< "\n -----------------------------------------------------\n";
	cout 
		<< "\n The Number of Transmitted Blocks = ";
	Present_Large_Number(block_number);

	cout
		<< "\n The Number of Transmitted Symbols = ";
	Present_Large_Number((__int64)Code.Row_Number *  (__int64)block_number);

	cout
		<< "\n Start SNR [dB] = "	<< SNR_dB_start << " [dB] "
		<< " ?? Es/N0 [dB] = "		<< SNR_dB_start + 10 * log10(Code.Code_Rate) << " [dB] \n"
		<< " End SNR [dB] = "		<< SNR_dB_end	<< " [dB] "
		<< " ?? Es/N0 [dB] = "		<< SNR_dB_end + 10 * log10(Code.Code_Rate) << " [dB] \n"
		<< " Step SNR [dB] = "		<< SNR_dB_step	<< " [dB] \n"
		<< " Block Step = "			<< BlockStep	<< "\n"
		<< " Stack Size = ";
	Present_Large_Number(StackSize);
	cout
		<< "\n"
		<< " Error Block Threshold =  " << ErrorBlock_Thr;
	//Present_Large_Number(ErrorBlock_Thr);
	if (Eb_N0_1_or_Es_N0_2 == 1) cout << " \n Eb/No Test";
	else if (Eb_N0_1_or_Es_N0_2 == 2) cout << " \n Es/No Test";
	cout
		<< "\n\n -----------------------------------------------------\n"
		<< endl;
}


// For example, 1000 -> 1,000, 1234567 -> 1,234,567
void Present_Large_Number(__int64 input_number) {
	using namespace std;
	vector<size_t> number;
	while (input_number >= 10) {
		number.push_back(input_number % 10);
		input_number /= 10;
	}
	number.push_back(input_number);

	int initial_position = 0;
	if ((number.size() % 3) == 1) {
		initial_position = 1;
	}
	else if ((number.size() % 3) == 2) {
		initial_position = 2;
	}
	else if ((number.size() % 3) == 0) {
		initial_position = 3;
	}

	int i(number.size() - 1);
	for (int j(0); j < initial_position; ++j) {
		cout << number.at(i - j);
	}

	for (i -= initial_position; i >= 0; --i) {
		if (((i + 1) % 3) == 0) cout << ",";
		cout << number.at(i);
	}
}


void DECODING_INFO::CRC_GMatrix_Genenerator(int Mes_Length, vector<__int8> CRC, MATRIX<__int8> &G) {
	//int CRC_length = CRC.size() - 1;
	int Mes_CRC_Length = Mes_Length + CRC_length;
	G.Building_Empty_Matrix(Mes_Length, Mes_CRC_Length);
	for (int row = 0; row < Mes_Length; ++row) {
		G._matrix[row][row] = 1;
		vector <__int8> polymomial(Mes_CRC_Length);
		polymomial.at(row) = 1;
		int index = row;
		while (index < Mes_Length) {
			int temp;
			bool status = FALSE;
			for (int add = 0; add < CRC_length + 1; add++) {
				polymomial.at(index + add) ^= CRC.at(add);
				if (status == FALSE) {
					if (polymomial.at(index + add) == 1) {
						temp = index + add;
						status = TRUE;
					}
				}
			}
			index = temp;
		}
		for (int col = Mes_Length; col < Mes_CRC_Length; col++) {
			G._matrix[row][col] = polymomial.at(col);
		}
	}
}

vector<__int8> DECODING_INFO::Mes_To_MesCrc(vector<__int8> Message, MATRIX<__int8> G) {
	vector<__int8> CodeWord(G.Col_number);
	//cout << CodeWord.size() << Message.size() << endl;
	for (int col = 0; col < G.Col_number; col++) {
		for (int row = 0; row < G.Row_number; row++) {
			CodeWord.at(col) ^= (Message.at(row) & G._matrix[row][col]);
		}
	}
	return CodeWord;
}

MATRIX <__int8> DECODING_INFO::G_times_G(MATRIX <__int8> G1, MATRIX <__int8> G2) {
	MATRIX <__int8> G;
	if (G1.Col_number != G2.Row_number) cout << "Multiplication Failed! Matrice size do not fit!" << endl;
	else {
		G.Building_Empty_Matrix(G1.Row_number, G2.Col_number);
		for (int i = 0; i < G1.Row_number; ++i) {
			for (int j = 0; j < G2.Col_number; ++j) {
				for (int k = 0; k < G1.Col_number; ++k) {
					G._matrix[i][j] ^= (G1._matrix[i][k] & G2._matrix[k][j]);
				}
			}
		}
	}
	return G;
}