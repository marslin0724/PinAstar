#include "PinAstar/AStarDecode.h"

//

void A_star_1_stack(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)//A_star_laststack
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
	

	do {
		Pointer = Stack.at(0);

		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;

			if ((Child_Node.level == message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z <= pc_i)) {
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; 
					}
				}
				
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric >= Best_Goal.metric) 
							break;
					}
				}
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z < pc_i)) {
				if (Pointer.metric != Child_Node.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z == pc_i)) {
				for (size_t j(Child_Node.level); j < message_length; ++j) 
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這邊的G是Sorted_G
					}
				}

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;
				//
				
	
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric >= Best_Goal.metric) break;
					}
				}
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);
	
	do
	{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}

				if ((Child_Node.metric < Best_Goal.metric)){
					if (Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_3_stack(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

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
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_2_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
			//
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_4_stack(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			//
			//Systematic_Linear_Block_Code_Encoder(G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
			}

			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;

			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}
				if ((Child_Node.metric < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			Best_Goal_2 = Best_Goal;
			A_star_3_stack(G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_PC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
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
			}
			else if (
				(Child_Node.level < message_length) && 
				(Child_Node.metric < Best_Goal.metric) && 
				(Child_Node.D_z <= decoding_info.Constraint_i)) {

				Place_Node(Stack, Child_Node, decoding_info);	

			}
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

void A_star_PC_out(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) 
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do{
		Pointer = Stack.at(0);

		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}
		
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.Binary_STE;
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;
				//
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}

				for (size_t j(message_length); j < codeword_length; ++j) {
					Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (Child_Node.metric > Best_Goal.metric) break;
				}
				//
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i))
			{
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) 
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				
				//
				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				//
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				
				for (size_t j(message_length); j < codeword_length; ++j) {
					Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (Child_Node.metric > Best_Goal.metric) break;
				}
				//
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}
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
}

void A_star_2_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length), // 傳入下一個stack時用來告知目前最小的metric大小
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)){
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
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {
			
			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
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
			Best_Goal_2 = Best_Goal;
			A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t) (1 + 2), decoding_info); //3是因為第一層的pc-1加上第二層的pc-out-2
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

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
}

void A_star_3_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length), // 用來接收下一層 A* 回傳的結果
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
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
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

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
			Best_Goal_2 = Best_Goal;
			A_star_2_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
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

void A_star_4_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length), // 用來接收下一層 A* 回傳的結果
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
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
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
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
						++decoding_info.COM;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Stack.erase(Stack.begin());
			Best_Goal_2 = Best_Goal;
			A_star_3_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
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

void A_star_Parity_PC(MATRIX<__int8> &G, DECODING_INFO &decoding_info){
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
	
	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); 

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;

	do{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			decoding_info.STE++;
			
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z++;

			size_t y = 0;
			while (Q_table._matrix[Pointer.level][y] != 0) {
				__int8 x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
				y++;
				decoding_info.STE++;
			}

			if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				decoding_info.CandidateCodeWord++;

				Child_Node.metric = Child_Node.auxiliary();
				Child_Node.heuristic = 0;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) 
			{
				if (Child_Node.metric == Pointer.metric) { 
					Stack.insert(Stack.begin(), Child_Node);
					decoding_info.COM++;
				}
				else Place_Node(Stack, Child_Node, decoding_info);	
			}
		}
	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
}

void A_star_Parity_PC_out(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table);

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;

	do{
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			decoding_info.STE++;

			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z++;

			size_t y = 0;
			while (Q_table._matrix[Pointer.level][y] != 0) {

				__int8 x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
				y++;
				decoding_info.STE++;
			}

			if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				decoding_info.CandidateCodeWord++;

				Child_Node.metric = Child_Node.auxiliary();
				Child_Node.heuristic = 0;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				
				for (unsigned int j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits[j] = Hard_RX[j];
					decoding_info.STE++;
				}

				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;

				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (unsigned int j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
				Child_Node.heuristic = 0;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);

			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				Place_Node(Stack, Child_Node, decoding_info);
			}
		}
	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
}

void A_star_1_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		stack_size(G.Row_number + 1),
		message_length(G.Row_number),
		codeword_length(G.Col_number);
	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);
	do{
		Pointer = Stack.front();
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit){
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
			++decoding_info.STE;

			size_t y (0);
			while (Q_table._matrix[Pointer.level][y] != 0){
				__int8 x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
				++y;
				++decoding_info.STE;
			}

			if ((Child_Node.level == message_length) 
				&& (Child_Node.auxiliary() < Best_Goal.metric)
				&& (Child_Node.D_z <= pc_i)) {

				++decoding_info.CandidateCodeWord;
				Child_Node.metric = Child_Node.auxiliary();
				Child_Node.heuristic = 0;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) 
				&& (Child_Node.auxiliary() < Best_Goal.metric) 
				&& (Child_Node.D_z < pc_i)) {

				Place_Node(Stack, Child_Node, decoding_info);
			}
			else if ((Child_Node.level < message_length) 
				&& (Child_Node.auxiliary() < Best_Goal.metric) 
				&& (Child_Node.D_z == pc_i)) {
				
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits[j] = Hard_RX[j];
					decoding_info.STE++;
				}

				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;

				Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		stack_size(G.Row_number + 1),
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;

			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;

			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0;
				while (Q_table._matrix[Pointer.level][y] != 0) {
					__int8 x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_Parity(G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_3_stack_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		stack_size(G.Row_number + 1),
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);

				if (new_bit != Hard_RX[Pointer.level])
					Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0;
				while (Q_table._matrix[Pointer.level][y] != 0) {
					__int8 x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Best_Goal_2 = Best_Goal;
			A_star_2_stack_Parity(G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_Parity_2_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
	
	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來
	
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	do {
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {
			decoding_info.STE += 2;
			//
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level))	Child_Node.D_z = Pointer.D_z + 1;
				
				// compute heuristic value
				unsigned int y(0);
				while (Q_table._matrix[Pointer.level][y] != 0) {
					unsigned int x(Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node));
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
					++y;
					++decoding_info.STE;
				}
				//
				if ((Child_Node.auxiliary() < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 3, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
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

void A_star_Parity_3_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(message_length, 0);

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
	for (size_t i(0); i < message_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	do {
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;

			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;

			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level])
					Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0;
				while (Q_table._matrix[Pointer.level][y] != 0) {
					__int8 x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Best_Goal_2 = Best_Goal;
			A_star_2_stack_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
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

void A_star_Parity_4_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level])
					Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0;
				while (Q_table._matrix[Pointer.level][y] != 0) {
					__int8 x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Best_Goal_2 = Best_Goal;
			A_star_3_stack_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
	} while ((!Stack.empty()));

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

void A_star_1_stack_PC(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack[0];
		if (Pointer.level == (message_length - 1))
			Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z++;

			decoding_info.STE++;

			if ((Child_Node.level == message_length) 
				&& (Child_Node.metric < Best_Goal.metric) 
				&& (Child_Node.D_z <= pc_i)) {

				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;

				Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) 
				&& (Child_Node.metric < Best_Goal.metric) 
				&& (Child_Node.D_z <= pc_i)) {

				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					decoding_info.COM++;
				}
			}
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_multiple_stack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;

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
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
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
				A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Multi_Stack_Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack, Multi_Stack);
			} // 整個stack的大小為STACK Size +1 
			//Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
		
		if (Stack.empty() && !Multi_Stack.empty()) {
			do {
				Best_Goal_2 = Best_Goal;
				Multi_Stack_Pointer = Multi_Stack.at(0);
				Multi_Stack.erase(Multi_Stack.begin());
				if (Multi_Stack_Pointer.metric < Best_Goal_2.metric) {
					A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Multi_Stack_Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info);
					if ((Best_Goal_2.auxiliary()) < Best_Goal.metric) {
						Best_Goal = Best_Goal_2; // Update the best candidate c odeword.
					}
				}
			} while (!Multi_Stack.empty());
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

void A_star_2_multiple_stack(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info, vector<NODE_PATH> &Multi_Stack)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number);

	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Multi_Stack_Pointer(message_length),
		Best_Goal(Pre_Best_Goal),
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
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
						++decoding_info.COM;
					}
				}
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			//
			/*
			Best_Goal_2 = Best_Goal;
			A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			*/
			Place_Node_multistack(Multi_Stack, Pointer, decoding_info);
			if (Multi_Stack.size() == decoding_info.Multi_Stack_Size) {
				Best_Goal_2 = Best_Goal;
				Multi_Stack_Pointer = Multi_Stack.at(0);
				Multi_Stack.erase(Multi_Stack.begin());
				if (Multi_Stack_Pointer.D_z == 1) {
					A_star_2_multiple_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Multi_Stack_Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info, Multi_Stack);
				}
				else if (Multi_Stack_Pointer.D_z == 2) {
					A_star_1_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Multi_Stack_Pointer, Best_Goal_2, (size_t)(pc_i + 2), decoding_info);
				}
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack, Multi_Stack);
			}
			//
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}
