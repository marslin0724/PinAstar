#include "AStarDecode.h"

void A_star_PC_OSC_II (MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double First_nonzero_metric = 0;
	double metric_thr(0), Total_nonzero_metric(0), nonzero_metric_num(0), Avg_nonzero_metric(0);
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) {
			Hard_RX.at(i) = 1;
		}
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;
	bool Break_Flag(FALSE);

	do
	{
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
			++decoding_info.STE;

			if ((Child_Node.level == message_length) &&
				(Child_Node.metric < Best_Goal.metric) &&
				(Child_Node.D_z <= decoding_info.Constraint_i)) {

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack)) {
					if (Best_Goal.metric < metric_thr) {
						Break_Flag = TRUE;
						break;
					}
					if (First_nonzero_metric == 0) {		//PoHan
						First_nonzero_metric = Best_Goal.metric;
					}
					if (Best_Goal.metric <= (decoding_info.New_OSC_Alpha - (decoding_info.SNR * 0.05))* First_nonzero_metric) {
						Break_Flag = TRUE;
					}
				}
				
				if (Child_Node.metric != 0) {
					Total_nonzero_metric += Child_Node.metric;
					nonzero_metric_num++;
					Avg_nonzero_metric = Total_nonzero_metric / nonzero_metric_num;
				}
				
			}
			else if (
				(Child_Node.level < message_length) &&
				(Child_Node.metric < Best_Goal.metric) &&
				(Child_Node.D_z <= decoding_info.Constraint_i)) {

				Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		if (Break_Flag) break;
	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;
}

void A_star_PC_out_OSC_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double First_nonzero_metric = 0;
	double metric_thr(0);
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) {
			Hard_RX.at(i) = 1;
		}
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;
	bool Break_Flag(FALSE);

	do {
		Pointer = Stack.at(0);

		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z++;
			decoding_info.STE++;
			//
			if ((Child_Node.level == message_length)
				&& (Child_Node.metric < Best_Goal.metric)
				&& (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
				//
				if (Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack)) {
					if (Best_Goal.metric < metric_thr) {
						Break_Flag = TRUE;
						break;
					}
					if (First_nonzero_metric == 0) {		//PoHan
						First_nonzero_metric = Best_Goal.metric;
					}
					if (Best_Goal.metric <= (decoding_info.New_OSC_Alpha - (decoding_info.SNR * 0.05))* First_nonzero_metric) {
						Break_Flag = TRUE;
					}
				}
				

				/*double G_value = OSC_early_termination(Metric_Table, 22, codeword_seq, Hard_RX);
				if (Best_Goal.metric <= G_value) Break_Flag = TRUE;*/
			}
			else if ((Child_Node.level < message_length)
				&& (Child_Node.metric < Best_Goal.metric)
				&& (Child_Node.D_z == decoding_info.Constraint_i))
			{
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
					decoding_info.STE++;
				}
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
				//
				if (Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack)) {
					if (Best_Goal.metric < metric_thr) {
						Break_Flag = TRUE;
						break;
					}
					if (First_nonzero_metric == 0) {		//PoHan
						First_nonzero_metric = Best_Goal.metric;
					}
					if (Best_Goal.metric <= (decoding_info.New_OSC_Alpha - (decoding_info.SNR * 0.05))* First_nonzero_metric) {
						Break_Flag = TRUE;
					}
				}
				
				//
				

				/*double G_value = OSC_early_termination(Metric_Table, 22, codeword_seq, Hard_RX);
				if (Best_Goal.metric <= G_value) Break_Flag = TRUE;*/
			}
			else if ((Child_Node.level < message_length)
				&& (Child_Node.metric < Best_Goal.metric)
				&& (Child_Node.D_z < decoding_info.Constraint_i)) {
				//
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					decoding_info.COM++;
				}
			}
		}
		if (Break_Flag) break;
	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	//cout << endl << decoding_info.COM << endl;
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

}

void A_star_PC_out_CBC_OSC_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info)  // 修正部分狀況沒有套用DM-I的問題
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
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序
	//system("pause");

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	bool Update_Flag;
	double First_nonzero_metric = 0;
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);
	
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//decoding_info.DoubleDecoder = TRUE;
	// 開始 Tree Search
	//vector<__int8> BestGoalTemp;
	//decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());
	bool Break_Flag(FALSE);


	do {

		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
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
				if (Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, decoding_info)) {
					if (Best_Goal.metric < OSC_metric_thr) {
						Break_Flag = TRUE;
						break;
					}
					if (First_nonzero_metric == 0) {		//PoHan
						First_nonzero_metric = Best_Goal.metric;
					}
					if (Best_Goal.metric <= (decoding_info.New_OSC_Alpha - (decoding_info.SNR * 0.05))* First_nonzero_metric) {
						Break_Flag = TRUE;
						break;
					}
				}
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
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
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

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
							decoding_info);
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
							decoding_info);
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
				if (Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, decoding_info)) {
					if (Best_Goal.metric < OSC_metric_thr) {
						Break_Flag = TRUE;
						break;
					}
					if (First_nonzero_metric == 0) {		//PoHan
						First_nonzero_metric = Best_Goal.metric;
					}
					if (Best_Goal.metric <= (decoding_info.New_OSC_Alpha - (decoding_info.SNR * 0.05))* First_nonzero_metric) {
						Break_Flag = TRUE;
						break;
					}
				}
				//cout << "p";
				
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
		if (Break_Flag) break;

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

void A_star_PC_out_CBC_OSC_III(MATRIX<__int8> &G, DECODING_INFO &decoding_info) // Sufficient condition
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0), // rx_signal hard decision 結果
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<__int8> Sorted_H;
	MATRIX<double> Metric_Table(2, codeword_length);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length),
		temp_Node(message_length); //CBC 回傳的best node
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	int Min_Dz = INT_MAX;
	bool Update_Flag = FALSE;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	
	//cout << endl;
	// for OSD-i
	double OSC_metric_thr(0);
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	bool Break_Flag(FALSE);

	decoding_info.Counter = 0;
	do {
		Pointer = Stack.front();
		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			Update_Flag = FALSE;
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			else {
				Child_Node.Same_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) &&
				(Child_Node.metric < Best_Goal.metric) &&
				(Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}

				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
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


				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						Child_Node.D_z++;
						if (Child_Node.metric > Best_Goal.metric) break;
					}
					else {
						Child_Node.Same_Index.push_back(j);
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;

				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if ((Update_Flag == TRUE) && (Min_Dz >= Best_Goal.D_z)) {
					Min_Dz = Best_Goal.D_z;
					q_i = (dmin*1.2) - Best_Goal.D_z;
					if (q_i < 0) q_i = 0;
					for (size_t j(0); j < q_i; ++j) {
						Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
					}
					if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
					if (Best_Goal.metric < Sufficient_Condition_thr) Break_Flag = TRUE;
					Temp_Sufficient_Condition_thr = 0;
				}
			}
			else if ((Child_Node.level < message_length) &&
				(Child_Node.metric < Best_Goal.metric) &&
				(Child_Node.D_z == decoding_info.Constraint_i)) {
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

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
							decoding_info);
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
							decoding_info);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						Child_Node.D_z++;
						if (Child_Node.metric > Best_Goal.metric) break;
					}
					else {
						Child_Node.Same_Index.push_back(j);
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;

				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if ((Update_Flag == TRUE) && (Min_Dz >= Best_Goal.D_z)) {
					Min_Dz = Best_Goal.D_z;
					q_i = (dmin*1.2) - Best_Goal.D_z;
					if (q_i < 0) q_i = 0;
					for (size_t j(0); j < q_i; ++j) {
						Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
					}
					if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
					if (Best_Goal.metric < Sufficient_Condition_thr) Break_Flag = TRUE;
					Temp_Sufficient_Condition_thr = 0;
				}
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
			
			if (Best_Goal.metric < OSC_metric_thr) {
				Break_Flag = TRUE;
				break;
			}
			
		}
		if (Break_Flag == TRUE)break;
	} while (!Stack.empty());

	

	//
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

