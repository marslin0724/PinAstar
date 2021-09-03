#include "PinAstar/AStarDecode.h"

void A_star_PC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
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

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Best_Goal.metric < metric_thr) Break_Flag = TRUE;
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

void A_star_PC_out_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//
				if (Best_Goal.metric < metric_thr) Break_Flag = TRUE;

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
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//
				if (Best_Goal.metric < metric_thr) Break_Flag = TRUE;
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

void A_star_1_stack_OSC(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number);
	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Node);

	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
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

			if ((Child_Node.level == message_length) 
				&& (Child_Node.metric < Pre_Best_Goal.metric) 
				&& (Child_Node.D_z <= pc_i)) {
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
				//
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Best_Goal.metric < metric_thr) Break_Flag = TRUE;
			}
			else if ((Child_Node.level < message_length) 
				&& (Child_Node.metric < Pre_Best_Goal.metric) 
				&& (Child_Node.D_z < pc_i)) {

				if (Pointer.metric != Child_Node.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					decoding_info.COM++;
				}
			}
			else if ((Child_Node.level < message_length) 
				&& (Child_Node.metric < Pre_Best_Goal.metric) 
				&& (Child_Node.D_z == pc_i)) {
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
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Best_Goal.metric < metric_thr) Break_Flag = TRUE;
			}
			//
		}
		if (Break_Flag) break;
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}
	
	double best_G_value(0);
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
			for (size_t j(message_length); j < codeword_length; ++j)
				Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr)	break;
			
			// The method of early termination of OSD according to textbook.
			// but the effect of using this method is very not obvious

			//Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
			//double G_value = OSC_early_termination(Metric_Table, 22, codeword_seq, Hard_RX);
			/*double new_G_value = OSC_early_termination2(Metric_Table, 22, codeword_seq, Hard_RX, Pointer, message_length);
			if (new_G_value > best_G_value) best_G_value = new_G_value;
			if (Best_Goal.metric <= best_G_value) break;
			else {
				//size_t n(0);
				//while (Stack.at(n).metric < metric_thr) n++;
				//if (n < Stack.size())Stack.erase(Stack.begin() + n, Stack.end());

				unsigned int n(0), Stack_Size(Stack.size());
				while (n < Stack_Size) {
					if (Stack.at(n).metric >= best_G_value) {
						Stack.erase(Stack.begin() + n);
						Stack_Size--;
					}
					else n++;
				}
			}*/
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {
			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)) {

					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_OSC(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 3, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			// 
			if (Best_Goal.metric < metric_thr)	break;
			
			//Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
			//double G_value = OSC_early_termination(Metric_Table, 22, codeword_seq, Hard_RX);
			/*double new_G_value = OSC_early_termination2(Metric_Table, 22, codeword_seq, Hard_RX, Best_Goal_2, message_length);
			if (new_G_value > best_G_value) best_G_value = new_G_value;
			if (Best_Goal.metric <= best_G_value) break;
			else {
				//size_t n(0);
				//while (Stack.at(n).metric < metric_thr) n++;
				//if (n < Stack.size())Stack.erase(Stack.begin() + n, Stack.end());
				unsigned int n(0), Stack_Size(Stack.size());
				while (n < Stack_Size) {
					if (Stack.at(n).metric >= best_G_value) {
						Stack.erase(Stack.begin() + n);
						Stack_Size--;
					}
					else n++;
				}
			}*/
		}
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
}

void A_star_3_stack_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

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
		Best_Goal_2(message_length), //用來接收 下一層 A* 所跑的最佳結果
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}
	
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	double best_G_value(0);

	double metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			Stack.erase(Stack.begin());
			//
			Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
			for (size_t j(message_length); j < codeword_length; ++j)
				Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Best_Goal.metric < metric_thr)	break;
			
			/*double new_G_value = OSC_early_termination2(Metric_Table, 22, codeword_seq, Hard_RX, Pointer, message_length);
			if (new_G_value > best_G_value) best_G_value = new_G_value;
			if (Best_Goal.metric <= new_G_value) break;*/
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level))
					Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_2_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr)	break;
			
			/*double new_G_value = OSC_early_termination2(Metric_Table, 22, codeword_seq, Hard_RX, Pointer, message_length);
			if (new_G_value > best_G_value) best_G_value = new_G_value;
			if (Best_Goal.metric <= new_G_value) break;*/
		}
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
}

void A_star_4_stack_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector<size_t>
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
		Best_Goal_2(message_length), // 用來接收下一層 A* 所回傳結果
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
			for (size_t j(message_length); j < codeword_length; ++j)
				Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);

				if (new_bit != Hard_RX.at(Pointer.level))
					Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_3_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, 2, decoding_info);
			//
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}
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
}

void A_star_2_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	bool Jump_flag = FALSE;

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}

	// OSC
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
			//
			if (Best_Goal.metric < metric_thr) break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_DM_I(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 3, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}

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
}

void A_star_3_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	bool Jump_flag = FALSE;

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}

	// OSC
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
			if (Best_Goal.metric < metric_thr) break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;
				
				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_2_stack_DM_I(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}

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
}

void A_star_4_stack_DM_I_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Best_Goal_2(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	bool Jump_flag = FALSE;

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}

	// OSC
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
			if (Best_Goal.metric < metric_thr) break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_3_stack_DM_I(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}

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
}


void A_star_2_multiple_stack_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);
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
		Multi_Stack_Pointer(message_length),
		Best_Goal(message_length),
		Best_Goal_2(message_length), // 傳入下一個stack時用來告知目前最小的metric大小
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	vector<NODE_PATH> Multi_Stack;
	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double metric_thr(0);
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			for (size_t j(message_length); j < codeword_length; ++j) {
				Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
				if (Pointer.metric > Best_Goal.metric) break;
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Best_Goal.metric < metric_thr)	break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}
				if ((Child_Node.metric < Best_Goal.metric)) {
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Stack.erase(Stack.begin());
			//
			//Best_Goal_2 = Best_Goal;
			//A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info); //3是因為第一層的pc-1加上第二層的pc-out-2
			Place_Node_multistack(Multi_Stack, Pointer, decoding_info);
			if (Multi_Stack.size() == decoding_info.Multi_Stack_Size) {
				Best_Goal_2 = Best_Goal;
				Multi_Stack_Pointer = Multi_Stack.at(0);
				Multi_Stack.erase(Multi_Stack.begin());
				A_star_1_stack_OSC(Sorted_G, Metric_Table, Hard_RX, Multi_Stack_Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack, Multi_Stack);
				if (Best_Goal.metric < metric_thr)	break;
			} // 整個stack的大小為STACK Size +1 
			//Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

		if (Stack.empty() && !Multi_Stack.empty()) {
			do {
				Best_Goal_2 = Best_Goal;
				Multi_Stack_Pointer = Multi_Stack.at(0);
				Multi_Stack.erase(Multi_Stack.begin());
				if (Multi_Stack_Pointer.metric < Best_Goal_2.metric) {
					A_star_1_stack_OSC(Sorted_G, Metric_Table, Hard_RX, Multi_Stack_Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info);
					if ((Best_Goal_2.auxiliary()) < Best_Goal.metric) {
						Best_Goal = Best_Goal_2; // Update the best candidate c odeword.
					}
					if (Best_Goal.metric < metric_thr)	break;
				}
			} while (!Multi_Stack.empty());
			if (Best_Goal.metric < metric_thr)	break;
		}

	} while ((!Stack.empty()) || (!Multi_Stack.empty()));
	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

// methods of early termination of OSD
// Book : Lin Shu , Chapter 10, Reliability-Based Soft-Decision Decoding Algorithms for Linear Block Codes
double OSC_early_termination(MATRIX<double> &Metric_Table, size_t d_min, vector<__int8> &codeword, vector<__int8> &z_seq) {
	
	size_t q1(d_min);
	size_t Codeword_Length((size_t)z_seq.size());
	vector<size_t> D0_Index(Codeword_Length, -1);

	size_t counter(0);
	for (size_t i(0); i < Codeword_Length; i++) {
		if (codeword.at(i) != z_seq.at(i)) q1--;
		else D0_Index.at(counter++) = i;
	}

	if (q1 < 0) q1 = 0;

	double G_value(0);
	for (size_t i(counter - 1); i >= (counter - q1); i--) {
		G_value += (Metric_Table._matrix[0][D0_Index.at(i)] + Metric_Table._matrix[1][D0_Index.at(i)]);
	}

	return G_value;
}

double OSC_early_termination2(
	MATRIX<double>		&Metric_Table, 
	size_t				d_min,
	vector<__int8>		&codeword,
	vector<__int8>		&z_seq,
	NODE_PATH			&candidate, 
	size_t				message_length) {
	
	size_t delta(0), k(message_length), n((size_t)z_seq.size());
	double G_value(0);
	vector<size_t> D0_Index(n, -1);
	
	for (size_t j = 1; j <= (candidate.D_z + 1); j++)
		G_value += (Metric_Table._matrix[0][k - j] + Metric_Table._matrix[1][k - j]);
	
	// determine set D_0 and the size of D_1
	size_t counter(0), D1_size(0);
	for (size_t i(0); i < n; i++) {
		if (codeword.at(i) != z_seq.at(i)) D1_size++;
		else D0_Index.at(counter++) = i;
	}

	// determine delta
	delta = d_min - D1_size - (candidate.D_z + 1);
	if (delta < 0) delta = 0;

	for (size_t i(counter - 1); i >= (counter - delta); i--) {
		G_value += (Metric_Table._matrix[0][D0_Index.at(i)] + Metric_Table._matrix[1][D0_Index.at(i)]);
	}

	return G_value;
}

// combined version 
// DM-I + Parity-Heuristic + OSC
// lower STE but its speed is almost unchanged.
void A_star_3_stack_DM_I_Parity_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX.at(i) = 0;
		else Hard_RX.at(i) = 1;
	}

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	// OSC
	double metric_thr(0);
	for (size_t i(0); i < codeword_length; i++) {
		metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			decoding_info.CandidateCodeWord++;
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;

				size_t y(0);
				__int8 x(0);
				while ((Q_table._matrix[Pointer.level][y] != 0)) {
					x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

					if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
						if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

						if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
							decoding_info.DeletedNode_counter++;
							break;
						}
					}
					y++;
					decoding_info.STE++;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) <= decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Best_Goal_2 = Best_Goal;
			//
			A_star_2_stack_DM_I_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, (size_t)2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
			if (Best_Goal.metric < metric_thr) break;
		}

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
}