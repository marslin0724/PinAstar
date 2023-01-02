#include "AStarDecode.h"
#include <iostream>
#include <vector>
void A_star_RSS_Test(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		message_state_number(pow(2, message_length));;

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		Sorted_codeword(codeword_length, 0);

	
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	Best_Goal.metric = FLT_MAX;

	vector<NODE_PATH> Stack(1, Pointer);
	int Distance = 0, Minimum_Distance = INT_MAX;
		
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	double OSC_metric_thr(0), Metric_difference(0);
	int D_z(0);
	for (size_t i(0); i < message_state_number; ++i) {		//�p��minimum distance
		//���w message sequence
		for (size_t j(0); j < message_length; ++j)
			message_seq.at(message_length - j - 1) = (i >> j) & 1;

		//��l�� codeword sequence
		codeword_seq.clear();
		Systematic_Linear_Block_Code_Encoder(G, message_seq, codeword_seq);
		Distance = 0;
		for (size_t k(0); k < codeword_length; k++) {
			if (codeword_seq.at(k) != 0)Distance++;
		}
		if ((Distance < Minimum_Distance) && Distance != 0) Minimum_Distance = Distance;
		std::cout << std::endl;
		for (size_t i(0); i < message_length; ++i) std::cout << " " << message_seq.at(i);
		std::cout << "Minimum Distance : " << Minimum_Distance;
	};
	/*
	for (size_t i(0); i < (size_t)codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		if ((Hard_RX.at(i) != Sorted_codeword.at(i)) && (i < message_length)) {
			Metric_difference += abs(decoding_info.rx_signal_seq.at(i));
			D_z++;
		}
			// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	*/
	decoding_info.New_OSC_Alpha = Metric_difference / OSC_metric_thr;
	if (decoding_info.Worst_OSC_Ratio < decoding_info.New_OSC_Alpha) decoding_info.Worst_OSC_Ratio = decoding_info.New_OSC_Alpha;
	/*
	if (D_z == 0)decoding_info.Dz_0_Number++;
	else if (D_z == 1)decoding_info.Dz_1_Number++;
	else if (D_z == 2)decoding_info.Dz_2_Number++;
	else if (D_z == 3)decoding_info.Dz_3_Number++;
	else if (D_z == 4)decoding_info.Dz_4_Number++;
	else if (D_z == 5)decoding_info.Dz_5_Number++;
	else if (D_z == 6)decoding_info.Dz_6_Number++;
	else if (D_z == 7)decoding_info.Dz_7_Number++;
	else if (D_z == 8)decoding_info.Dz_8_Number++;
	else if (D_z == 9)decoding_info.Dz_9_Number++;
	else if (D_z == 10)decoding_info.Dz_10_Number++;
	else if (D_z == 11)decoding_info.Dz_11_Number++;
	else if (D_z == 12)decoding_info.Dz_12_Number++;
	else if (D_z == 13)decoding_info.Dz_13_Number++;
	else if (D_z == 14)decoding_info.Dz_14_Number++;
	else if (D_z == 15)decoding_info.Dz_15_Number++;
	else if (D_z == 16)decoding_info.Dz_16_Number++;
	else if (D_z == 17)decoding_info.Dz_17_Number++;
	else if (D_z == 18)decoding_info.Dz_18_Number++;
	else if (D_z == 19)decoding_info.Dz_19_Number++;
	else if (D_z == 20)decoding_info.Dz_20_Number++;
	else if (D_z == 21)decoding_info.Dz_21_Number++;
	else if (D_z == 22)decoding_info.Dz_22_Number++;
	else if (D_z == 23)decoding_info.Dz_23_Number++;
	*/
	/*
	do {
	Pointer = Stack.at(0);

	if (Pointer.level == (message_length - 1))
		Stack.erase(Stack.begin());

	for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
		Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
		++decoding_info.STE;

		if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)) {
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
			for (size_t j(message_length); j < codeword_length; ++j)
				Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
			//
			Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
		}
		else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)) {
			if (Child_Node.metric != Pointer.metric)
				Place_Node(Stack, Child_Node, decoding_info);
			else {
				Stack.at(0) = Child_Node;
				++decoding_info.COM;
			}
		}
	}
	} while (!Stack.empty());
	*/


	//
	//Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	//Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	decoding_info.estimated_codeword = decoding_info.code_seq;
	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	
}

void A_star_PC_out_Dynamic_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)  // �ץ��������p�S���M��DM-I�����D
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq ��sorting rx���G
	// Location_index �����Ƨ�
	//system("pause");

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	double Constraint_i(0); //PoHan
	bool Update_Flag;
	size_t Update_Num(0), Non_Update_Num(0);
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence�Һ�X��codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // ��XOSC threshold
	

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//decoding_info.DoubleDecoder = TRUE;
	// �}�l Tree Search
	//vector<__int8> BestGoalTemp;
	//decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());

	do {

		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// �o�̪�pointer���O�u��pointer, �u�O��pointer�h����Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // ��pointer��level�F��k-1����, ���U�Ӫ����child node���|�Ok, �]���b�o�ӨB�J��stack�̤W�����ȵ��d��(��pointer�s��Ƥ���N��pop��������)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Update_Flag = FALSE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.metric < decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //�o�̥u�����control�Ӥw�ҥH�S����systemetic encoder
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, decoding_info);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
					if (Best_Goal.metric >= OSC_metric_thr) {
						decoding_info.num_Best_Node_Update += decoding_info.temp_num_Deleted_Candidate_in_Stack;
					}
				}

				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.metric >= decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq����Fi��message������extension��level k, �é�����level n��seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							OSC_metric_thr,
							1,
							decoding_info);
					/*
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
					*/
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							OSC_metric_thr,
							decoding_info);
					/*
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
					*/
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
				}
				//cout << "p";
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Place_Node(Stack, Child_Node, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;

	} while (!Stack.empty());
	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);

	//cout << decoding_info.code_seq.size() << "," << decoding_info.Sorted_R.size() << endl;
	//system("pause");


	//Error Test
	/*
	for (int i = 0; i < 192; ++i) {
		//cout << decoding_info.Sorted_R.at(i) << ",";
		if (decoding_info.Sorted_R.at(i) < 0 && decoding_info.code_seq.at(i) == 0) decoding_info.Error_Accumulation.at(i)++;
		else if (decoding_info.Sorted_R.at(i) > 0 && decoding_info.code_seq.at(i) == 1) decoding_info.Error_Accumulation.at(i)++;
	}*/
	//cout << endl;

	/*
	int Total_Error = 0;
	int MRIP_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size() / 2; ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++MRIP_Error;
			//cout << "index: " << i << ", ";
		}
	}
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++Total_Error;
			//cout << "index: " << i << ", ";
		}
	}
	//cout << MRIP_Error << endl;
	*/
	/*
	int Total_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != codeword_seq.at(i)) {
			++Total_Error;
			cout << "index: (" << i << ", " << decoding_info.Sorted_R.at(i) << ") / ";
		}
	}*/

	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}