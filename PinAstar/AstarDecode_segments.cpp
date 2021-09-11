#include "AStarDecode.h"

void A_star_Segment(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{
	//test
	/*
	decoding_info.rx_signal_seq = { 0.6, 1.8, -1.1, -1.3, 0.9, 1.6, -1.0, -0.5 };
	G._matrix = { 
		{1,0,0,0,1,0,1,1},
		{0,1,0,0,1,1,1,0},
		{0,0,1,0,1,1,0,1},
		{0,0,0,1,0,1,1,1} };
		*/
	// 
	//暫定segments 為2
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / 2);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

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
	vector<NODE_PATH> S_Stack1,S_Stack2,C_Stack1,C_Stack2;	//C_Stack為存放達到segment_length的node
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
	//先讓C_Stack1與C_Stack2有node
	Pre_Procedure_Segment(decoding_info, Metric_Table,S_Stack1,C_Stack1, segment_length,0);
	Pre_Procedure_Segment(decoding_info, Metric_Table,S_Stack2,C_Stack2, segment_length, segment_length);
	//結合兩個node成一個candidate codeword
	Combine_Segment(C_Stack1.at(0), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
		G, Sorted_G, Metric_Table);
	if(C_Stack1.size() >= 2)
		Combine_Segment(C_Stack1.at(1), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
	/*
	for (int i = 0; i < 4;i++) {
		cout << (int)Best_Goal.message_bits.at(i) << ",";
	}
	cout << endl;*/
	while ((!S_Stack1.empty()) || (!S_Stack2.empty())){
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
			if ( Pointer.level == (segment_length - 1))
				S_Stack1.erase(S_Stack1.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				++decoding_info.STE;

				if ((Child_Node.level == segment_length) && (Child_Node.metric < Best_Goal.metric)) {
					Place_C_Stack(C_Stack1, Child_Node, decoding_info);
					Combine_Segment(Child_Node, C_Stack1, C_Stack2,S_Stack1,S_Stack2 ,Best_Goal, segment_length, decoding_info,
							G, Sorted_G, Metric_Table);
					
				}
				else if ((Child_Node.level < segment_length) && (Child_Node.metric < Best_Goal.metric)) {
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
				++decoding_info.STE;
				
				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)) {
					Place_C_Stack(C_Stack2, Child_Node, decoding_info);
					Combine_Segment(Child_Node, C_Stack2, C_Stack1, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
							G, Sorted_G, Metric_Table);
					
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)) {
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

inline void Pre_Procedure_Segment(	DECODING_INFO& decoding_info,
							MATRIX<double>& Metric_Table,
							vector<NODE_PATH>& S_Stack,vector<NODE_PATH>& C_Stack,
							size_t segment_length , size_t offset) 
{
	NODE_PATH
		Pointer(segment_length),
		Child_Node(segment_length);
	do {
		Pointer = S_Stack.at(0);

		if (Pointer.level == (offset + segment_length - 1))
			S_Stack.erase(S_Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			++decoding_info.STE;

			if ((Child_Node.level == segment_length + offset) ) {
				if (C_Stack.size() == 0) C_Stack.push_back(Child_Node);
				else Place_C_Stack(C_Stack, Child_Node, decoding_info);
			}
			else if ((Child_Node.level < segment_length + offset)) {
				if (Child_Node.metric != Pointer.metric) {
					Place_Node(S_Stack, Child_Node, decoding_info);
				}
				else {
					S_Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}
	} while (C_Stack.empty());
}

inline void Combine_Segment(	NODE_PATH &Update_node, vector<NODE_PATH>& C_Stack_Update,vector<NODE_PATH> &C_Stack, vector<NODE_PATH>& S_Stack1, vector<NODE_PATH>& S_Stack2, 
						NODE_PATH &Best_Goal,size_t segment_length,
						DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		start, end;
	NODE_PATH
		Combine_node(message_length);
	vector<__int8>
		codeword_seq(codeword_length);
	double metric_bound = Best_Goal.metric - Update_node.metric;
	Combine_node.message_bits.assign(Update_node.message_bits.begin(), Update_node.message_bits.end());
	if (Update_node.level == message_length) {
		start = 0, end = segment_length;
	}
	else {
		start = segment_length, end = message_length;
	}
	for (int i = 0; i < C_Stack.size(); i++) {
		if (C_Stack.at(i).metric > metric_bound) {
			break;
		}
		for (int j = start; j < end; j++) {
			Combine_node.message_bits.at(j) = C_Stack.at(i).message_bits.at(j);
		}
		Combine_node.metric = Update_node.metric + C_Stack.at(i).metric;
		/*
		for (int i = 0; i < 4; i++) {
			cout << (int)Combine_node.message_bits.at(i) << ",";
		}
		cout << "metric = " << (double)Combine_node.metric << endl;
		*/
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//
		Systematic_Linear_Block_Code_Encoder(Sorted_G, Combine_node.message_bits, codeword_seq);
		for (size_t j(message_length); j < codeword_length; ++j)
			Combine_node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
		//
		if (Update_Best_Goal_Procedure(Combine_node, Best_Goal, C_Stack)) {
			Best_Goal.metric += 1;
			Update_Best_Goal_Procedure(Combine_node, Best_Goal, C_Stack_Update);
			Best_Goal.metric += 1;
			Update_Best_Goal_Procedure(Combine_node, Best_Goal, S_Stack1);
			Best_Goal.metric += 1;
			Update_Best_Goal_Procedure(Combine_node, Best_Goal, S_Stack2);
			//更新bound
			metric_bound = Best_Goal.metric - Update_node.metric;
		}
		

	}
}

inline void Place_C_Stack(vector<NODE_PATH>& Stack, NODE_PATH& child_node, DECODING_INFO& decoding_info)//把node按照順序放入stack
{
	int
		position_number(Stack.size()),
		Stack_Size(Stack.size());
	
	for (int j(Stack_Size - 1); j >= 0; --j) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).metric < child_node.metric) {	//arrange stack by original metric
			position_number = j + 1;
			break;
		}
	}

	// the number of compared nodes
	decoding_info.COM = decoding_info.COM + ((double)Stack_Size -  position_number + 1);
	if (position_number == Stack_Size) {
		Stack.push_back(child_node);
		--decoding_info.COM;
	}
	else Stack.insert(Stack.begin() + position_number, child_node);

	// If the current capacity of stack is larger than the configured stack size
	// , then delete the last node which stores in the stack.
	if (Stack_Size > decoding_info.StackSize) Stack.pop_back(); // 整個stack的大小為STACK Size +1 

	if (decoding_info.MaxUsed_Stack < Stack_Size)
		decoding_info.MaxUsed_Stack = (unsigned int)Stack_Size;
}