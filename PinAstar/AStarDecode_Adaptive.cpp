#include "AStarDecode.h"
#include <iostream>
#include <vector>

void A_star_PC_out_CBC_OSC_Adaptive_i_NEW(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {//兩層

	Adaptive_I_Decode_Info Adaptive_Info;
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
	vector<NODE_PATH>
		Stack_CBC1(1, Pointer),
		Stack_CBC2(1, Pointer);

	Adaptive_Info.Hard_RX = Hard_RX;
	Adaptive_Info.MRIP_codeword = MRIP_codeword;
	Adaptive_Info.Sorted_G = Sorted_G;
	Adaptive_Info.Metric_Table = Metric_Table;
	Adaptive_Info.Best_Goal = Best_Goal;
	Adaptive_Info.Best_Goal.metric = DBL_MAX;

	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Adaptive_Info.Sorted_G, Location_Index, Adaptive_Info.Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	Adaptive_Info.OSC_metric_thr = 0;

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Adaptive_Info.Metric_Table._matrix[0][i] != 0) Adaptive_Info.Hard_RX.at(i) = 1;

		// for OSC threshold
		Adaptive_Info.OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq.assign(Adaptive_Info.Hard_RX.begin(), Adaptive_Info.Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, message_seq, Adaptive_Info.MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	Adaptive_Info.OSC_metric_thr = decoding_info.OSC_Alpha*Adaptive_Info.OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		Adaptive_Info.OSC_metric_thr += Adaptive_Info.Metric_Table._matrix[0][i] + Adaptive_Info.Metric_Table._matrix[1][i];
	}
	decoding_info.DoubleDecoder = TRUE;

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')
	
	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
	decoding_info.phase = 1;
	A_star_PC_out_CBC_OSC_Block_First(G, decoding_info, Adaptive_Info, Stack_CBC1);  //   1st
	//cout << "(1):" << decoding_info.Counter << endl;
	//cout << "(1):" << decoding_info.Counter << endl;
	if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
		//cout << "(X2)" << endl;
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
		decoding_info.phase = 2;
		//A_star_PC_out_CBC_OSC_Block_New(G, decoding_info, Adaptive_Info, Stack_CBC1, Stack_CBC2); //   2nd
		A_star_PC_out_CBC_OSC_Block_Last(G, decoding_info, Adaptive_Info, Stack_CBC1);  //    2rd 2132
		//cout << "(2):" << decoding_info.Counter << endl;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder3_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder3_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder2_i;
			//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
			decoding_info.phase = 3;
			A_star_PC_out_CBC_OSC_Block_Last(G, decoding_info, Adaptive_Info, Stack_CBC2);  //    3rd
		}
	}
	decoding_info.CBC_FlippingBit = Temp_CBC;
	decoding_info.Constraint_i = Temp_i;

	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, Adaptive_Info.Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	decoding_info.STE = decoding_info.STE / (double)message_length;
	/*
	decoding_info.STE_1 = decoding_info.STE_1 / (double)message_length;
	decoding_info.STE_2 = decoding_info.STE_2 / (double)message_length;
	decoding_info.STE_3 = decoding_info.STE_3 / (double)message_length;
	*/
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

void A_star_PC_out_CBC_OSC_Block_First(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	bool First_Flag = TRUE;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	
	//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
	// 開始 Tree Search
	do {
		//cout << Stack.size() << endl;
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
			//cout << "C:" << Stack.size() << endl;
			// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
			
			// End
			Extend_Node_Procedure(Pointer, Child_Node, Adaptive_info.Metric_Table, new_bit);
			if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now
			//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "(A1)";
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
				// End
				++decoding_info.Counter;
				codeword_seq = Adaptive_info.MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
						Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				if (First_Flag == TRUE) {
					Stack_CBC1.at(0) = Child_Node;
					First_Flag = FALSE;
				}
				else {
					//Place_Node(Stack_CBC1, Child_Node, decoding_info);
					Stack_CBC1.insert(Stack_CBC1.begin() + position_number, Child_Node);
				}
				position_number++;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = Adaptive_info.MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}
					//cout << Child_Node.metric << endl;
					if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
				}
				else Child_Node = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "(A3)";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}
		if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
			//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
			break;
		}
	} while (!Stack.empty());
	//cout << decoding_info.Counter << endl;
	if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	//cout << endl << endl;
}

void A_star_PC_out_CBC_OSC_Block_New(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_Current, vector<NODE_PATH> &Stack_Next) {
	{
		size_t
			message_length(G.Row_number),
			codeword_length(G.Col_number),
			position_number(0),
			error_counter(0);

		vector<__int8>
			codeword_seq(codeword_length, 0),
			message_seq(message_length, 0);

		NODE_PATH
			Pointer(message_length),
			Child_Node(message_length),
			temp_Node(message_length);

		bool First_Flag = TRUE;
		bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
		int Level_k_previous, Difference;


		if (decoding_info.Cancelled_Candidate_i != 0) {
			Operater_Deletion = TRUE;
			Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
			Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
		}

		//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
		// 開始 Tree Search
		do {
			//cout << Stack.size() << endl;
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer = Stack_Current.at(0);
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack_Current.erase(Stack_Current.begin());
			}
			//++decoding_info.Counter;
			if ((Pointer.level < message_length) && (Pointer.metric < Adaptive_info.Best_Goal.metric) && (Pointer.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Pointer.level); j < message_length; ++j) {
					Pointer.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = Adaptive_info.MRIP_codeword;
				for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Pointer,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Pointer,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Pointer,
							Adaptive_info.Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				if (Operater_Deletion == FALSE || Pointer.D_z > decoding_info.Cancelled_Candidate_i) {
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Pointer.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Pointer.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}
					//cout << Pointer.metric << endl;
					if (temp_Node.metric < Pointer.metric)Pointer = temp_Node;
				}
				else Pointer = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer, Adaptive_info.Best_Goal, Stack_Current);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
					//cout << "C:" << Stack.size() << endl;
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer.level == Level_k_previous) && (Pointer.D_z <= Difference)) {
						if (Pointer.level != (message_length - 1)) Stack_Current.erase(Stack_Current.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure(Pointer, Child_Node, Adaptive_info.Metric_Table, new_bit);
					if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
						++Child_Node.D_z;
						Child_Node.Diff_Index.push_back(Pointer.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now
					//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
					// Reach level k
					if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
						//cout << "(A1)";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = Adaptive_info.MRIP_codeword;

						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
							}
						}
						error_counter = Child_Node.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
								Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
								if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
							}
						}

						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_Current);
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
						//cout << "(A2)";
						//cout << "B";
						++decoding_info.Counter;
						if (First_Flag == TRUE) {
							Stack_Next.at(0) = Child_Node;
							First_Flag = FALSE;
						}
						else {
							//Place_Node(Stack_Next, Child_Node, decoding_info);
							Stack_Next.insert(Stack_Next.begin() + position_number, Child_Node);
						}
						position_number++;
						for (__int16 j(Child_Node.level); j < message_length; ++j) {
							Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = Adaptive_info.MRIP_codeword;
						for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Adaptive_info.Sorted_G,
									Adaptive_info.Metric_Table,
									codeword_seq,
									Adaptive_info.Hard_RX,
									Child_Node,
									Adaptive_info.Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Adaptive_info.Sorted_G,
									Adaptive_info.Metric_Table,
									codeword_seq,
									Adaptive_info.Hard_RX,
									Child_Node,
									Adaptive_info.Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Adaptive_info.Sorted_G,
									Adaptive_info.Metric_Table,
									codeword_seq,
									Adaptive_info.Hard_RX,
									Child_Node,
									Adaptive_info.Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
									Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
									if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
						}
						else Child_Node = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_Current);
						//cout << Adaptive_info.Best_Goal.metric << endl;
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
						//cout << "(A3)";
						if (Child_Node.metric != Pointer.metric)
							Place_Node(Stack_Current, Child_Node, decoding_info);
						else {
							Stack_Current.at(0) = Child_Node;
							++decoding_info.COM;
						}
					}
				}
			}
			if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
				//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
				break;
			}
		} while (!Stack_Current.empty());
		//cout << decoding_info.Counter << endl;
		if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
		//cout << endl << endl;
	}
}

void A_star_PC_out_CBC_OSC_Block_Last(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);



	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	if (decoding_info.Cancelled_Candidate_i != 0) {
		Operater_Deletion = TRUE;
		Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
		Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
	}
	//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
	// 開始 Tree Search
	do {
		//cout << Stack.size() << endl;
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack_CBC1.at(0);
		//Stack_CBC1.erase(Stack_CBC1.begin());
		
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack_CBC1.erase(Stack_CBC1.begin());
		}
		
		//++decoding_info.Counter;
		if ((Pointer.level < message_length) && (Pointer.metric < Adaptive_info.Best_Goal.metric) && (Pointer.D_z == decoding_info.Constraint_i)) {
			//cout << "(A2)";
			//cout << "B";
			++decoding_info.Counter;
			for (__int16 j(Pointer.level); j < message_length; ++j) {
				Pointer.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
			}
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = Adaptive_info.MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (__int16 j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}
			//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

			decoding_info.STE += (codeword_length - Pointer.level);
			//++decoding_info.CandidateCodeWord;

			// CBC
			if (decoding_info.CBC_FlippingBit == 1) {
				temp_Node
					= Control_Band_Check_1bit(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 2) {
				temp_Node
					= Control_Band_Check_2bits(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 3) {
				temp_Node
					= Control_Band_Check_3bits(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info,
						1);
			}
			else temp_Node.metric = DBL_MAX;
			if (Operater_Deletion == FALSE || Pointer.D_z > decoding_info.Cancelled_Candidate_i) {
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
						Pointer.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Pointer.metric > Adaptive_info.Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer.metric)Pointer = temp_Node;
			}
			else Pointer = temp_Node;
			//cout << Adaptive_info.Best_Goal.metric << endl;
			Update_Best_Goal_Procedure(Pointer, Adaptive_info.Best_Goal, Stack_CBC1);
			//cout << Adaptive_info.Best_Goal.metric << endl;
		}
		else {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
				//cout << "C:" << Stack.size() << endl;
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if ((Operater_Deletion == TRUE) && (Pointer.level == Level_k_previous) && (Pointer.D_z <= Difference)) {
					if (Pointer.level != (message_length - 1)) Stack_CBC1.erase(Stack_CBC1.begin());
					//cout << "!";
					break;
				}
				// End
				Extend_Node_Procedure(Pointer, Child_Node, Adaptive_info.Metric_Table, new_bit);
				if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}
				++decoding_info.STE;
				++decoding_info.Binary_STE;
				// Child_Node: The node we are examining now
				//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
				// Reach level k
				if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					//cout << "(A1)";
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
					// End
					++decoding_info.Counter;
					codeword_seq = Adaptive_info.MRIP_codeword;

					// DM-I: Reach Control level to check hamming distance
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					error_counter = Child_Node.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
					}
					decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}

					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_CBC1);
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
					//cout << "(A2)";
					//cout << "B";
					++decoding_info.Counter;
					for (__int16 j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

					codeword_seq = Adaptive_info.MRIP_codeword;
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node.level);
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info,
								1);
					}
					else temp_Node.metric = DBL_MAX;
					if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
								Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
								if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
							}
						}
						//cout << Child_Node.metric << endl;
						if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
					}
					else Child_Node = temp_Node;
					//cout << Adaptive_info.Best_Goal.metric << endl;
					Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_CBC1);
					//cout << Adaptive_info.Best_Goal.metric << endl;
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
					//cout << "(A3)";
					if (Child_Node.metric != Pointer.metric)
						Place_Node(Stack_CBC1, Child_Node, decoding_info);
					else {
						//Place_Node(Stack_CBC1, Child_Node, decoding_info);
						Stack_CBC1.at(0) = Child_Node;
						++decoding_info.COM;
						
					}
				}
			}
			if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
				//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
				break;
			}
		}
	} while (!Stack_CBC1.empty());
	//cout << decoding_info.Counter << endl;
	if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	//cout << endl << endl;
}

void A_star_PC_out_CBC_OSC_Adaptive_i_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {//兩層

	Adaptive_I_Decode_Info Adaptive_Info;
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
	MATRIX<double> Fano_Metric_Table(2, message_length);


	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH>
		Stack_CBC1(1, Pointer),
		Stack_CBC2(1, Pointer);

	Adaptive_Info.Hard_RX = Hard_RX;
	Adaptive_Info.MRIP_codeword = MRIP_codeword;
	Adaptive_Info.Sorted_G = Sorted_G;
	Adaptive_Info.Metric_Table = Metric_Table;
	Adaptive_Info.Fano_Metric_Table = Fano_Metric_Table;
	Adaptive_Info.Best_Goal = Best_Goal;
	Adaptive_Info.Best_Goal.metric = DBL_MAX;

	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Adaptive_Info.Sorted_G, Location_Index, Adaptive_Info.Metric_Table, Adaptive_Info.Fano_Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	Adaptive_Info.OSC_metric_thr = 0;

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Adaptive_Info.Metric_Table._matrix[0][i] != 0) Adaptive_Info.Hard_RX.at(i) = 1;

		// for OSC threshold
		Adaptive_Info.OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq.assign(Adaptive_Info.Hard_RX.begin(), Adaptive_Info.Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, message_seq, Adaptive_Info.MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	Adaptive_Info.OSC_metric_thr = decoding_info.OSC_Alpha*Adaptive_Info.OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		Adaptive_Info.OSC_metric_thr += Adaptive_Info.Metric_Table._matrix[0][i] + Adaptive_Info.Metric_Table._matrix[1][i];
	}
	decoding_info.DoubleDecoder = TRUE;

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
	decoding_info.phase = 1;
	A_star_PC_out_CBC_OSC_Block_First(G, decoding_info, Adaptive_Info, Stack_CBC1);  //   1st
	//cout << "(1):" << decoding_info.Counter << endl;
	//cout << "(1):" << decoding_info.Counter << endl;
	if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
		//cout << "(X2)" << endl;
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
		decoding_info.phase = 2;
		//A_star_PC_out_CBC_OSC_Block_New(G, decoding_info, Adaptive_Info, Stack_CBC1, Stack_CBC2); //   2nd
		A_star_PC_out_CBC_OSC_Block_Last(G, decoding_info, Adaptive_Info, Stack_CBC1);  //    2rd 2132
		//cout << "(2):" << decoding_info.Counter << endl;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder3_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder3_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder2_i;
			//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
			decoding_info.phase = 3;
			A_star_PC_out_CBC_OSC_Block_Last(G, decoding_info, Adaptive_Info, Stack_CBC2);  //    3rd
		}
	}
	decoding_info.CBC_FlippingBit = Temp_CBC;
	decoding_info.Constraint_i = Temp_i;

	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, Adaptive_Info.Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	decoding_info.STE = decoding_info.STE / (double)message_length;
	/*
	decoding_info.STE_1 = decoding_info.STE_1 / (double)message_length;
	decoding_info.STE_2 = decoding_info.STE_2 / (double)message_length;
	decoding_info.STE_3 = decoding_info.STE_3 / (double)message_length;
	*/
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

void A_star_PC_out_CBC_OSC_Block_Fano_First(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	bool First_Flag = TRUE;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除

	//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
	// 開始 Tree Search
	do {
		//cout << Stack.size() << endl;
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
			//cout << "C:" << Stack.size() << endl;
			// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***

			// End
			Extend_Node_Procedure_Fano(Pointer, Child_Node, Adaptive_info.Metric_Table, Adaptive_info.Fano_Metric_Table, new_bit);
			if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now
			//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "(A1)";
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
				// End
				++decoding_info.Counter;
				codeword_seq = Adaptive_info.MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
						Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				if (First_Flag == TRUE) {
					Stack_CBC1.at(0) = Child_Node;
					First_Flag = FALSE;
				}
				else {
					//Place_Node(Stack_CBC1, Child_Node, decoding_info);
					Stack_CBC1.insert(Stack_CBC1.begin() + position_number, Child_Node);
				}
				position_number++;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = Adaptive_info.MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}
					//cout << Child_Node.metric << endl;
					if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
				}
				else Child_Node = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "(A3)";
				Place_Node_Fano(Stack, Child_Node, decoding_info);
				/*
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
				*/
			}
		}
		if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
			//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
			break;
		}
	} while (!Stack.empty());
	//cout << decoding_info.Counter << endl;
	if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	//cout << endl << endl;
}

void A_star_PC_out_CBC_OSC_Block_Fano_Last(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info, vector<NODE_PATH> &Stack_CBC1) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);



	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	if (decoding_info.Cancelled_Candidate_i != 0) {
		Operater_Deletion = TRUE;
		Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
		Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
	}
	//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
	// 開始 Tree Search
	do {
		//cout << Stack.size() << endl;
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack_CBC1.at(0);
		Stack_CBC1.erase(Stack_CBC1.begin());
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack_CBC1.erase(Stack_CBC1.begin());
		}
		*/
		//++decoding_info.Counter;
		if ((Pointer.level < message_length) && (Pointer.metric < Adaptive_info.Best_Goal.metric) && (Pointer.D_z == decoding_info.Constraint_i)) {
			//cout << "(A2)";
			//cout << "B";
			++decoding_info.Counter;
			for (__int16 j(Pointer.level); j < message_length; ++j) {
				Pointer.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
			}
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = Adaptive_info.MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (__int16 j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}
			//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

			decoding_info.STE += (codeword_length - Pointer.level);
			//++decoding_info.CandidateCodeWord;

			// CBC
			if (decoding_info.CBC_FlippingBit == 1) {
				temp_Node
					= Control_Band_Check_1bit(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 2) {
				temp_Node
					= Control_Band_Check_2bits(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 3) {
				temp_Node
					= Control_Band_Check_3bits(
						Adaptive_info.Sorted_G,
						Adaptive_info.Metric_Table,
						codeword_seq,
						Adaptive_info.Hard_RX,
						Pointer,
						Adaptive_info.Best_Goal,
						decoding_info,
						1);
			}
			else temp_Node.metric = DBL_MAX;
			if (Operater_Deletion == FALSE || Pointer.D_z > decoding_info.Cancelled_Candidate_i) {
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
						Pointer.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Pointer.metric > Adaptive_info.Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer.metric)Pointer = temp_Node;
			}
			else Pointer = temp_Node;
			//cout << Adaptive_info.Best_Goal.metric << endl;
			Update_Best_Goal_Procedure(Pointer, Adaptive_info.Best_Goal, Stack_CBC1);
			//cout << Adaptive_info.Best_Goal.metric << endl;
		}
		else {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
				//cout << "C:" << Stack.size() << endl;
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if ((Operater_Deletion == TRUE) && (Pointer.level == Level_k_previous) && (Pointer.D_z <= Difference)) {
					if (Pointer.level != (message_length - 1)) Stack_CBC1.erase(Stack_CBC1.begin());
					//cout << "!";
					break;
				}
				// End
				Extend_Node_Procedure_Fano(Pointer, Child_Node, Adaptive_info.Metric_Table, Adaptive_info.Fano_Metric_Table, new_bit);
				if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}
				++decoding_info.STE;
				++decoding_info.Binary_STE;
				// Child_Node: The node we are examining now
				//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
				// Reach level k
				if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					//cout << "(A1)";
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
					// End
					++decoding_info.Counter;
					codeword_seq = Adaptive_info.MRIP_codeword;

					// DM-I: Reach Control level to check hamming distance
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					error_counter = Child_Node.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
					}
					decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}

					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_CBC1);
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
					//cout << "(A2)";
					//cout << "B";
					++decoding_info.Counter;
					for (__int16 j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

					codeword_seq = Adaptive_info.MRIP_codeword;
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node.level);
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Adaptive_info.Sorted_G,
								Adaptive_info.Metric_Table,
								codeword_seq,
								Adaptive_info.Hard_RX,
								Child_Node,
								Adaptive_info.Best_Goal,
								decoding_info,
								1);
					}
					else temp_Node.metric = DBL_MAX;
					if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
								Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
								if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
							}
						}
						//cout << Child_Node.metric << endl;
						if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
					}
					else Child_Node = temp_Node;
					//cout << Adaptive_info.Best_Goal.metric << endl;
					Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack_CBC1);
					//cout << Adaptive_info.Best_Goal.metric << endl;
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
					//cout << "(A3)";
					if (Child_Node.metric != Pointer.metric)
						Place_Node_Fano(Stack_CBC1, Child_Node, decoding_info);
					else {
						Place_Node_Fano(Stack_CBC1, Child_Node, decoding_info);
						/*
						Stack_CBC1.at(0) = Child_Node;
						++decoding_info.COM;
						*/
					}
				}
			}
			if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
				//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
				break;
			}
		}
	} while (!Stack_CBC1.empty());
	//cout << decoding_info.Counter << endl;
	if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	//cout << endl << endl;
}

void A_star_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		distance(0);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Best_Goal_2(message_length), // 用來接收下一層 A* 回傳的結果
		Child_Node(message_length);

	deque<NODE_PATH> Stack1,Stack2;
	deque<NODE_PATH>* Stack, *Stack_next,*tmp_Stack;
	Stack = &Stack1, Stack_next = &Stack2;
	Stack->push_back(Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	double metric_thr(0);
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) {
			Hard_RX.at(i) = 1;
		}
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack->front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	do {
		Pointer = Stack->front();
		Stack->pop_front();
		if (Pointer.metric >= Best_Goal.metric) continue;
		//flipping node and push them into queue

		if (Pointer.D_z < decoding_info.Constraint_i && Pointer.level < message_length) {
			for (int i = Pointer.level; i < message_length; i++) {
				NODE_PATH tmp = Pointer;
				tmp.message_bits.at(i) ^= 1;
				++tmp.D_z;
				tmp.Diff_Index.push_back(i);
				tmp.level = i + 1;
				if (Metric_Table._matrix[tmp.message_bits.at(i)][i]) {
					tmp.metric += Metric_Table._matrix[tmp.message_bits.at(i)][i];
				}
				if (Pointer.metric < Best_Goal.metric)
					Stack_next->push_back(tmp);
			}
			decoding_info.STE += 2 * (message_length - Pointer.level);
		}
		//Update best goal
		codeword_seq = MRIP_codeword;
		for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
			codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
			for (size_t j(message_length); j < codeword_length; ++j) {
				codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
			}
		}

		for (size_t j(message_length); j < codeword_length; ++j) {
			Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < metric_thr) break;
		}
		
		if (Stack->empty()) {
			if (distance >= 3 && Best_Goal.metric < 1.2 * metric_thr) break;
			distance++;
			tmp_Stack = Stack;
			Stack = Stack_next;
			Stack_next = tmp_Stack;
			sort(Stack->begin(),Stack->end(),
				[&](const NODE_PATH& a, const NODE_PATH& b) {
					decoding_info.COM++;
					return a.metric < b.metric;
				});
		}
	} while (!Stack->empty());

	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;
}

