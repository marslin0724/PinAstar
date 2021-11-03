#include "AStarDecode.h"

void A_star_repeatition(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0),
		rep_idx(G.Row_number,0);
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

	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);

	Pre_Procedure_rep(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table,rep_idx);

	do {
		Pointer = Stack.at(0);

		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure_rep(Pointer, Child_Node, Metric_Table, new_bit, rep_idx);
			++decoding_info.STE;

			if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric)) {
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				//
				Child_Node.heuristic = 0;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric)) {
				Place_Node_rep(Stack, Child_Node, decoding_info);
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
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void Extend_Node_Procedure_rep(
	NODE_PATH		&Pointer,
	NODE_PATH		&child_node,
	MATRIX<double>	&metric_table,
	__int8			&new_bit,
	vector<size_t>& rep_idx)
{
	child_node = Pointer;
	child_node.message_bits.at(Pointer.level) = new_bit;
	++child_node.level;
	child_node.metric += metric_table._matrix[new_bit][child_node.level - 1];
	child_node.heuristic +=	metric_table._matrix[new_bit][rep_idx.at(child_node.level - 1)];
}

void Place_Node_rep(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info)//把node按照順序放入stack
{
	size_t
		position_number(Stack.size()),
		Stack_Size(Stack.size());

	for (size_t j(0); j < Stack_Size; ++j) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).auxiliary() > child_node.metric) {	//arrange stack by original metric
			position_number = j;
			break;
		}
	}

	// the number of compared nodes
	decoding_info.COM = decoding_info.COM + ((double)position_number + 1);
	if (position_number == Stack_Size) --decoding_info.COM;

	Stack.insert(Stack.begin() + position_number, child_node);

	// If the current capacity of stack is larger than the configured stack size
	// , then delete the last node which stores in the stack.
	if (Stack_Size > decoding_info.StackSize) Stack.pop_back(); // 整個stack的大小為STACK Size +1 

	if (decoding_info.MaxUsed_Stack < Stack_Size)
		decoding_info.MaxUsed_Stack = (unsigned int)Stack_Size;
}


void Pre_Procedure_rep(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	vector<size_t>		&rep_idx)
{
	vector<double> sorting_rx_signal_seq(rx_signal);
	vector<double> pair_idx(rx_signal.size(),0);
	vector<double> sorting_pair_idx(rx_signal.size(), 0);
	for (int i = 0; i < pair_idx.size()/2; i++) {
		pair_idx.at(i) = pair_idx.at(i + pair_idx.size() / 2) = i;
	}
	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation(rx_signal, G, Sorted_G, Permutation);

	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	//
	Sort_Function(pair_idx, Permutation, sorting_pair_idx);
	for (int i = 0; i < rep_idx.size(); i++) {
		rep_idx.at(i) = sorting_pair_idx.at(i);
	}
	//reuse pair_idx 
	for (int i = rep_idx.size(); i < sorting_pair_idx.size(); i++) {
		pair_idx.at(sorting_pair_idx.at(i)) = i;
	}
	//put pair index
	for (int i = 0; i < rep_idx.size(); i++) {
		rep_idx.at(i) = pair_idx.at(rep_idx.at(i));
	}
	//
	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	//Sorted_R(sorting_rx_signal_seq);
}