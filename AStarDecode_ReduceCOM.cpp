#include "PinAstar/AStarDecode.h"

// 2018/9/1 
// some methods for trying to reduce the avg.COM of A*-I decoder
// A_star_I_New_Comparsion -> huge saving of COM, but the execution time is almost unchange.
// A_star_I_QuickSort -> higher COM, it is used to present the difference if we apply another method
// to find the stack position of child nodes.

void Place_Node_Qucik_Sorting(vector<NODE_PATH> &stack, NODE_PATH child_node, DECODING_INFO &decoding_info)
{
	stack.insert(stack.begin(), child_node);
	QuickSort(stack, 0, ((unsigned int)stack.size() - 1), decoding_info.COM);

	if (stack.size() > decoding_info.StackSize) stack.pop_back();
	if (decoding_info.MaxUsed_Stack < stack.size())
		decoding_info.MaxUsed_Stack = (unsigned int)stack.size();
}

void Place_Node_New(vector<NODE_PATH> &Stack, NODE_PATH &Child_Node, DECODING_INFO &decoding_info, size_t &pre_position)
{
	size_t position_number = 0;
	if (pre_position > 0) {

		if (Stack[pre_position].metric > Child_Node.metric) // ¡ô{

			for (size_t j(pre_position); j > 0; j--) {
				decoding_info.COM++;
			
				if (Stack[j].metric < Child_Node.metric) {
					position_number = j;
					pre_position = position_number;
					break;
				}
				position_number = 1;
			}
		else {
			for (size_t j(pre_position); j < Stack.size(); j++){
				decoding_info.COM++;
			
				if (Stack[j].metric > Child_Node.metric) {
					position_number = j;
					pre_position = position_number;
					break;
				}
				position_number = (size_t)Stack.size();
			}
		}
		Stack.insert(Stack.begin() + position_number, Child_Node);
		if (Stack.size() > decoding_info.StackSize) Stack.pop_back();
	}
	else {
		Place_Node(Stack, Child_Node, decoding_info);
	}

	if (decoding_info.MaxUsed_Stack < Stack.size())
		decoding_info.MaxUsed_Stack = (size_t)Stack.size();
}


// Functions for Qucik Sorting 
/************************************/
void swap(NODE_PATH &a, NODE_PATH &b) {
	NODE_PATH temp = a;
	a = b;
	b = temp;
}
size_t Partition(vector<NODE_PATH> &arr, size_t front, size_t end, double &ComNN_temp) {
	NODE_PATH pivot = arr.at(end);
	size_t i = front - 1;

	for (size_t j(front); j < end; j++) {
		ComNN_temp++;
		if (arr.at(j).metric < pivot.metric) {
			i++;
			swap(arr.at(i), arr.at(j));
		}
	}
	i++;
	swap(arr.at(i), arr.at(end));
	return i;
}
void QuickSort(vector<NODE_PATH> &arr, size_t front, size_t end, double &ComNN_temp) {
	ComNN_temp++;
	if (front < end) {
		size_t pivot = Partition(arr, front, end, ComNN_temp);
		QuickSort(arr, front, pivot - 1, ComNN_temp);
		QuickSort(arr, pivot + 1, end, ComNN_temp);
	}
}

void A_star_I_New_Comparsion(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

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

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	size_t pre_position = -1; // bug
	do{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);

			decoding_info.STE++;

			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)) {
				
				decoding_info.STE += (codeword_length - message_length);
				pre_position = -1;

				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)) {
				
				if (Child_Node.metric==Pointer.metric){
					Stack.insert(Stack.begin(), Child_Node);
					decoding_info.COM++;
				}
				else{
					Place_Node_New(Stack, Child_Node, decoding_info, pre_position);
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

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_I_QuickSort(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

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

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq,G,Sorted_G,Location_Index,Metric_Table);

	do{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit){
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			decoding_info.STE++;
			
			if ((Child_Node.level == message_length) && (Child_Node.metric <= Best_Goal.metric)){
				
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq); 
				for (size_t j(message_length); j < codeword_length; ++j)
					Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

				decoding_info.STE += (codeword_length - message_length);
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric <= Best_Goal.metric)) { 

				if (Child_Node.metric == Pointer.metric) {

					Stack.insert(Stack.begin(), Child_Node);
					decoding_info.COM++;
				}
				else Place_Node_Qucik_Sorting(Stack, Child_Node, decoding_info);	
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

