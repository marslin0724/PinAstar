#include "AStarDecode.h"

void A_star_section_PC(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / 2);

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

	Best_Goal.metric = FLT_MAX;
	/*
	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);
	*/
	vector<NODE_PATH> S_Stack1, S_Stack2;	//C_Stack為存放達到segment_length的node
	vector<NODE_COMB> C_Stack1, C_Stack2;
	S_Stack1.reserve(decoding_info.StackSize + 1);
	S_Stack2.reserve(decoding_info.StackSize + 1);
	C_Stack1.reserve(decoding_info.StackSize + 1);
	C_Stack2.reserve(decoding_info.StackSize + 1);
	S_Stack1.push_back(Pointer);
	Pointer.level = segment_length;
	S_Stack2.push_back(Pointer);
	//選擇Stack的參考
	size_t stack_flag = 2;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	
	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	while ((!S_Stack1.empty()) || (!S_Stack2.empty())) {
		if ((!S_Stack1.empty()) && (!S_Stack2.empty())) {
			if (S_Stack1.at(0).metric < S_Stack2.at(0).metric) {
				Pointer = S_Stack1.at(0);
				stack_flag = 1;
			}
			else {
				Pointer = S_Stack2.at(0);
				stack_flag = 2;
			}
			++decoding_info.COM;
		}
		else if (!S_Stack1.empty()) {
			Pointer = S_Stack1.at(0);
			stack_flag = 1;
		}
		else {
			Pointer = S_Stack2.at(0);
			stack_flag = 2;
		}
		if (stack_flag == 1) {
			if (Pointer.level == (segment_length - 1))
				S_Stack1.erase(S_Stack1.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section1_i)) {
					Combine_Segment(Child_Node, C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section1_i)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack1, Child_Node, decoding_info);
					else {
						S_Stack1.at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
			}
		}
		else {
			if (Pointer.level == (message_length - 1))
				S_Stack2.erase(S_Stack2.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section2_i)) {
					Combine_Segment(Child_Node, C_Stack2, C_Stack1, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section2_i)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack2, Child_Node, decoding_info);
					else {
						S_Stack2.at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
			}
		}
	}

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_section_PC_out(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / 2);

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

	Best_Goal.metric = FLT_MAX;
	/*
	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);
	*/
	vector<NODE_PATH> S_Stack1, S_Stack2;	//C_Stack為存放達到segment_length的node
	vector<NODE_COMB> C_Stack1, C_Stack2;
	S_Stack1.reserve(decoding_info.StackSize + 1);
	S_Stack2.reserve(decoding_info.StackSize + 1);
	C_Stack1.reserve(decoding_info.StackSize + 1);
	C_Stack2.reserve(decoding_info.StackSize + 1);
	S_Stack1.push_back(Pointer);
	Pointer.level = segment_length;
	S_Stack2.push_back(Pointer);
	//選擇Stack的參考
	size_t stack_flag = 2;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	while ((!S_Stack1.empty()) || (!S_Stack2.empty())) {
		if ((!S_Stack1.empty()) && (!S_Stack2.empty())) {
			if (S_Stack1.at(0).metric < S_Stack2.at(0).metric) {
				Pointer = S_Stack1.at(0);
				stack_flag = 1;
			}
			else {
				Pointer = S_Stack2.at(0);
				stack_flag = 2;
			}
			++decoding_info.COM;
		}
		else if (!S_Stack1.empty()) {
			Pointer = S_Stack1.at(0);
			stack_flag = 1;
		}
		else {
			Pointer = S_Stack2.at(0);
			stack_flag = 2;
		}
		if (stack_flag == 1) {
			if (Pointer.level == (segment_length - 1))
				S_Stack1.erase(S_Stack1.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section1_i)) {
					Combine_Segment(Child_Node, C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z < decoding_info.section1_i)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack1, Child_Node, decoding_info);
					else {
						S_Stack1.at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
				else if ((Child_Node.level < segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z == decoding_info.section1_i)) {
					decoding_info.STE += (segment_length - Child_Node.level);
					while (Child_Node.level < segment_length) {
						Child_Node.message_bits.at(Child_Node.level) = Hard_RX.at(Child_Node.level);
						Child_Node.level++;
					}
					Combine_Segment(Child_Node, C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
				}
			}
		}
		else {
			if (Pointer.level == (message_length - 1))
				S_Stack2.erase(S_Stack2.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section2_i)) {
					Combine_Segment(Child_Node, C_Stack2, C_Stack1, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)
						&& (Child_Node.D_z < decoding_info.section2_i)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack2, Child_Node, decoding_info);
					else {
						S_Stack2.at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z == decoding_info.section2_i)) {
					decoding_info.STE += (message_length - Child_Node.level);
					while (Child_Node.level < message_length) {
						Child_Node.message_bits.at(Child_Node.level) = Hard_RX.at(Child_Node.level);
						Child_Node.level++;
					}
					Combine_Segment(Child_Node, C_Stack2, C_Stack1, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
				}
			}
		}
	}

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_Section_Multistack(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / 2);

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

	Best_Goal.metric = FLT_MAX;
	/*
	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);
	*/
	vector<NODE_PATH> S_Stack1, S_Stack2;	//C_Stack為存放達到segment_length的node
	vector<NODE_COMB> C_Stack1, C_Stack2;
	S_Stack1.reserve(decoding_info.StackSize + 1);
	S_Stack2.reserve(decoding_info.StackSize + 1);
	C_Stack1.reserve(decoding_info.StackSize + 1);
	C_Stack2.reserve(decoding_info.StackSize + 1);
	S_Stack1.push_back(Pointer);
	Pointer.level = segment_length;
	S_Stack2.push_back(Pointer);
	//選擇Stack的參考
	size_t stack_flag = 2;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;

	while ((!S_Stack1.empty()) || (!S_Stack2.empty())) {
		if ((!S_Stack1.empty()) && (!S_Stack2.empty())) {
			if (S_Stack1.at(0).metric < S_Stack2.at(0).metric) {
				Pointer = S_Stack1.at(0);
				stack_flag = 1;
			}
			else {
				Pointer = S_Stack2.at(0);
				stack_flag = 2;
			}
			++decoding_info.COM;
		}
		else if (!S_Stack1.empty()) {
			Pointer = S_Stack1.at(0);
			stack_flag = 1;
		}
		else {
			Pointer = S_Stack2.at(0);
			stack_flag = 2;
		}
		if (stack_flag == 1) {
			if ((Pointer.level == segment_length) && (Pointer.D_z <= 1)) {
				//
				S_Stack1.erase(S_Stack1.begin());
				//			
				Combine_Segment(Pointer, C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
					G, Sorted_G, Metric_Table);
			}
			else if ((Pointer.level < segment_length) && (Pointer.D_z < 1)) {

				decoding_info.STE += 2;
				for (__int8 new_bit(0); new_bit < 2; new_bit++) {
					Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
					if (new_bit != Hard_RX.at(Pointer.level)) {
						++Child_Node.D_z;
					}

					if ((Child_Node.metric < Best_Goal.metric)) {
						if (Pointer.metric != Child_Node.metric)
							Place_Node(S_Stack1, Child_Node, decoding_info);
						else {
							S_Stack1.at(0) = Child_Node;
							++decoding_info.COM;
						}
					}
				}
			}
			else if ((Pointer.level < segment_length) && (Pointer.D_z == 1)) {
				S_Stack1.erase(S_Stack1.begin());
				Best_Goal_2 = Best_Goal;
				if (1 < decoding_info.section1_i - 2)
					A_star_Section1_Multistack(G,Sorted_G, Metric_Table, Hard_RX,
					Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info, C_Stack1, C_Stack2);
				else
					A_star_Section1_Multistack(G, Sorted_G, Metric_Table, Hard_RX,
					Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info, C_Stack1, C_Stack2);
				
				if (Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, S_Stack1)) {
					Update_Stack(Best_Goal, S_Stack2);
				}
			}
		}
		else {
			if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
				//
				S_Stack2.erase(S_Stack2.begin());
				//			
				Combine_Segment(Pointer, C_Stack2, C_Stack1, S_Stack1, S_Stack2, Best_Goal, message_length, decoding_info,
					G, Sorted_G, Metric_Table);
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

				decoding_info.STE += 2;
				for (__int8 new_bit(0); new_bit < 2; new_bit++) {
					Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
					if (new_bit != Hard_RX.at(Pointer.level)) {
						++Child_Node.D_z;
					}

					if ((Child_Node.metric < Best_Goal.metric)) {
						if (Pointer.metric != Child_Node.metric)
							Place_Node(S_Stack2, Child_Node, decoding_info);
						else {
							S_Stack2.at(0) = Child_Node;
							++decoding_info.COM;
						}
					}
				}
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
				S_Stack2.erase(S_Stack2.begin());
				Best_Goal_2 = Best_Goal;
				if (1 < decoding_info.section2_i - 2)
					A_star_Section2_Multistack(G,Sorted_G, Metric_Table, Hard_RX,
						Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info,C_Stack1, C_Stack2);
				else
					A_star_Section2_Multistack(G, Sorted_G, Metric_Table, Hard_RX,
						Pointer, Best_Goal_2, (size_t)(1 + 2), decoding_info, C_Stack1, C_Stack2);
				
				if (Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, S_Stack2)) {
					Update_Stack(Best_Goal, S_Stack1);
				}
			}
		}
		
	}

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_Section1_Multistack(MATRIX<__int8> &G, MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX,
	NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i,
	DECODING_INFO &decoding_info, vector<NODE_COMB> &C_Stack1, vector<NODE_COMB> &C_Stack2) {
	
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		section_length(message_length / 2);

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

		if ((Pointer.level == section_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//			
			Combine_Segment(Pointer, C_Stack1, C_Stack2, Stack, Stack, Best_Goal, section_length, decoding_info,
				G, Sorted_G, Metric_Table);
		}
		else if ((Pointer.level < section_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
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
		else if ((Pointer.level < section_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			//PC-1
			if (pc_i < decoding_info.section1_i) {
				Best_Goal_2 = Best_Goal;
				if (pc_i < decoding_info.section1_i - 2)
					A_star_Section1_Multistack(G, Sorted_G, Metric_Table, Hard_RX, Pointer,
						Best_Goal_2, pc_i + 1, decoding_info, C_Stack1, C_Stack2);
				else A_star_Section1_Multistack(G, Sorted_G, Metric_Table, Hard_RX, Pointer,
					Best_Goal_2, pc_i + 2, decoding_info, C_Stack1, C_Stack2);

			}//PC-out-2
			else {
				decoding_info.STE += (section_length - Pointer.level);
				while (Pointer.level < section_length) {
					Pointer.message_bits.at(Pointer.level) = Hard_RX.at(Pointer.level);
					Pointer.level++;
				}
				Combine_Segment(Pointer, C_Stack1, C_Stack2, Stack, Stack, Best_Goal, section_length, decoding_info,
					G, Sorted_G, Metric_Table);
			}
			//
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_Section2_Multistack(MATRIX<__int8> &G,MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX,
	 NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i,
	DECODING_INFO &decoding_info, vector<NODE_COMB> &C_Stack1, vector<NODE_COMB> &C_Stack2) {
	
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
			Combine_Segment(Pointer, C_Stack2, C_Stack1, Stack, Stack, Best_Goal, message_length, decoding_info,
				G, Sorted_G, Metric_Table);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
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
			//PC-1
			if (pc_i < decoding_info.section2_i) {
				Best_Goal_2 = Best_Goal;
				if(pc_i < decoding_info.section2_i-2) 
					A_star_Section2_Multistack(G, Sorted_G, Metric_Table, Hard_RX, Pointer,
					Best_Goal_2, pc_i + 1, decoding_info, C_Stack1, C_Stack2);
				else A_star_Section2_Multistack(G, Sorted_G, Metric_Table, Hard_RX, Pointer,
						Best_Goal_2, pc_i + 2, decoding_info, C_Stack1, C_Stack2);

			}//PC-out-2
			else {
				decoding_info.STE += (message_length - Pointer.level);
				while (Pointer.level < message_length) {
					Pointer.message_bits.at(Pointer.level) = Hard_RX.at(Pointer.level);
					Pointer.level++;
				}
				Combine_Segment(Pointer, C_Stack2, C_Stack1, Stack, Stack, Best_Goal, message_length, decoding_info,
					G, Sorted_G, Metric_Table);
			}
			//
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

//use MinHeap
void A_star_section_PC_MinHeap(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / 2);

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

	Best_Goal.metric = FLT_MAX;
	/*
	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);
	*/
	auto comp = [&](const NODE_PATH& a, const NODE_PATH& b) {
		decoding_info.COM++;
		return a.metric > b.metric;
	};
	typedef priority_queue<NODE_PATH, vector<NODE_PATH>, decltype(comp)> Q;
	Q S_Stack1(comp), S_Stack2(comp);
	//C_Stack為存放達到segment_length的node
	vector<NODE_COMB> C_Stack1, C_Stack2;
	S_Stack1.push(Pointer);
	Pointer.level = segment_length;
	S_Stack2.push(Pointer);
	//選擇Stack的參考
	size_t stack_flag = 2;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	while ((!S_Stack1.empty()) || (!S_Stack2.empty())) {
		if ((!S_Stack1.empty()) && (!S_Stack2.empty())) {
			if (S_Stack1.top().metric < S_Stack2.top().metric) {
				Pointer = S_Stack1.top();
				stack_flag = 1;
			}
			else {
				Pointer = S_Stack2.top();
				stack_flag = 2;
			}
			++decoding_info.COM;
		}
		else if (!S_Stack1.empty()) {
			Pointer = S_Stack1.top();
			stack_flag = 1;
		}
		else {
			Pointer = S_Stack2.top();
			stack_flag = 2;
		}
		if (Pointer.metric >= Best_Goal.metric) break;
		if (stack_flag == 1) {
			S_Stack1.pop();

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section1_i)) {
					Combine_Segment_MinHeap(Child_Node, C_Stack1, C_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < segment_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section1_i)) {
					if (Child_Node.metric != Pointer.metric)
						S_Stack1.push(Child_Node);
					else {
						int tmp = decoding_info.COM;
						S_Stack1.push(Child_Node);
						decoding_info.COM = tmp;
					}
				}
			}
		}
		else {
			S_Stack2.pop();

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section2_i)) {
					Combine_Segment_MinHeap(Child_Node, C_Stack2, C_Stack1, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);

				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)
					&& (Child_Node.D_z <= decoding_info.section2_i)) {
					if (Child_Node.metric != Pointer.metric)
						S_Stack2.push(Child_Node);
					else {
						int tmp = decoding_info.COM;
						S_Stack2.push(Child_Node);
						decoding_info.COM = tmp;
					}
				}
			}
		}
	}

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}


inline void Combine_Segment_MinHeap(NODE_PATH &Update_node, vector<NODE_COMB>& C_Stack_Update, vector<NODE_COMB> &C_Stack,
	NODE_PATH &Best_Goal, size_t segment_length,
	DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		start, end;
	double metric = 0;
	NODE_COMB
		Combine_node(codeword_length);
	vector<__int8>
		codeword_seq(codeword_length);
	double metric_bound = Best_Goal.metric - Update_node.metric;
	Combine_node.codeword_bits.assign(Update_node.message_bits.begin(), Update_node.message_bits.end());
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Combine_node.codeword_bits, Combine_node.codeword_bits);
	decoding_info.STE += (codeword_length - message_length);
	Combine_node.metric = Update_node.metric;
	Place_C_Stack(C_Stack_Update, Combine_node, decoding_info);
	for (int i = 0; i < C_Stack.size(); i++) {
		if (C_Stack.at(i).metric > metric_bound) {
			break;
		}
		for (int j = message_length; j < codeword_length; j++) {
			codeword_seq.at(j) = Combine_node.codeword_bits.at(j) ^ C_Stack.at(i).codeword_bits.at(j);
		}
		metric = Update_node.metric + C_Stack.at(i).metric;
		/*
		for (int i = 0; i < 4; i++) {
			cout << (int)Combine_node.message_bits.at(i) << ",";
		}
		cout << "metric = " << (double)Combine_node.metric << endl;
		*/
		//decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//

		for (size_t j(message_length); j < codeword_length; ++j)
			metric += Metric_Table._matrix[codeword_seq.at(j)][j];
		//
		if (metric < Best_Goal.metric) {
			for (int j = 0; j < message_length; j++) {
				Best_Goal.message_bits.at(j) = Combine_node.codeword_bits.at(j) ^ C_Stack.at(i).codeword_bits.at(j);
			}
			Best_Goal.metric = metric;
			Update_Stack(Best_Goal, C_Stack);
			Update_Stack(Best_Goal, C_Stack_Update);
			//更新bound
			metric_bound = Best_Goal.metric - Update_node.metric;

		}


	}
}