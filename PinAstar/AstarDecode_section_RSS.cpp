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