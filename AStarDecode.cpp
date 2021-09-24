#include "PinAstar/AStarDecode.h"



void Generate_Message(size_t seq_length, vector<__int8> &message_seq)
{
	// Randomly generate message bits
	/*std::default_random_engine seed;
	std::bernoulli_distribution bernoulli(0.5);
	
	size_t counter(message_seq.size());
	for (size_t i(0); i < counter; ++i)
		message_seq.at(i) = bernoulli(seed);*/


	std::random_device rd;
	std::mt19937 generator(rd());
	std::bernoulli_distribution bernoulli(0.5);

	auto gen = [&bernoulli, &generator]() {
		return bernoulli(generator);
	};
	generate(begin(message_seq), end(message_seq), gen);


} // end Generate_Message()

void BPSK_Modulation(vector<__int8> &input_seq , vector<double> &tx_signal_seq)
{
	// bit 0 -> symbol +1
	// bit 1 -> symbol -1 
	size_t counter(tx_signal_seq.size());
	for (size_t i(0); i < counter; ++i)
		tx_signal_seq.at(i) = 1. - 2. * input_seq.at(i);

} // end BPSK()

// Basic functions for A* decoding algorithm
// step 1 & 2
void Determine_Permutation(
	vector<double>			&Rx_signal, 
	MATRIX<__int8>			&G, 
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&permutation_seq)
{
	size_t 
		Register(0),
		col_number(0);

	// Initialize permutation sequence
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;
	
	// Compute permutation by using "quick sorting algorithm" 
	sort(
		permutation_seq.begin(), 
		permutation_seq.end(), 
		[&](const size_t& a, const size_t& b) 
		{return (abs(Rx_signal.at(a)) > abs(Rx_signal.at(b))); });
	
	// Col. exchange based on permutation sequence
	Sort_Matrix_Col(G, permutation_seq, sorted_G);
		
	// Obtain systematic form
	GJ_Elimination(sorted_G);
	
	// Check whether Sorted_G is systemtic form or not. 
	// If not, adjust unsigned interleaver sequence and Sorted_G
	// Check whether the first half of matrix is identity matrix or not
	for (size_t i(0); i < sorted_G.Row_number; ++i){	
		if (sorted_G._matrix[i][i] != 1){
			col_number = 0;
			
			// 同一列，由下一行開始，碰到第一個此列為 "1" 的行就交換 unsigned interleaver index
			for (size_t j(i + 1); j < sorted_G.Col_number; ++j){
				if (sorted_G._matrix[i][j] == 1){
					col_number = j;
					break;
				}
			}

			//Exchange the columns of permutation sequence
			Register = permutation_seq.at(i);
			permutation_seq.at(i) = permutation_seq.at(col_number);
			permutation_seq.at(col_number) = Register;
		
			// Exchange the columns of generator matrix
			Register = 0;
			for (size_t k(0); k < sorted_G.Row_number; ++k){
				Register = sorted_G._matrix[k][i];
				sorted_G._matrix[k][i] = sorted_G._matrix[k][col_number];
				sorted_G._matrix[k][col_number] = Register;
			}
		}
	}
} // end Determining_unsigned interleaver()*/

void Determine_Permutation_MultiBase(
	vector<double>			&Rx_signal,
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&permutation_seq,
	int						Total_Basis,
	int						Current_Basis)
{
	size_t
		Register(0),
		col_number(0);

	// Initialize permutation sequence
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	// Compute permutation by using "quick sorting algorithm" 
	sort(
		permutation_seq.begin(),
		permutation_seq.end(),
		[&](const size_t& a, const size_t& b)
	{return (abs(Rx_signal.at(a)) > abs(Rx_signal.at(b))); });

	// Col. exchange based on permutation sequence
	Sort_Matrix_Col(G, permutation_seq, sorted_G);

	// Obtain systematic form
	GJ_Elimination(sorted_G);

	/*
	for (size_t i(0); i < sorted_G.Row_number; i++) {
		for (size_t j(0); j < sorted_G.Col_number; j++) {
			cout << (int)sorted_G._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	*/

	// Check whether Sorted_G is systemtic form or not. 
	// If not, adjust unsigned interleaver sequence and Sorted_G
	// Check whether the first half of matrix is identity matrix or not
	for (size_t i(0); i < sorted_G.Row_number; ++i) {
		if (sorted_G._matrix[i][i] != 1) {
			col_number = 0;

			// 同一列，由下一行開始，碰到第一個此列為 "1" 的行就交換 unsigned interleaver index
			for (size_t j(i + 1); j < sorted_G.Col_number; ++j) {
				if (sorted_G._matrix[i][j] == 1) {
					col_number = j;
					break;
				}
			}

			//Exchange the columns of permutation sequence
			Register = permutation_seq.at(i);
			permutation_seq.at(i) = permutation_seq.at(col_number);
			permutation_seq.at(col_number) = Register;

			// Exchange the columns of generator matrix
			Register = 0;
			for (size_t k(0); k < sorted_G.Row_number; ++k) {
				Register = sorted_G._matrix[k][i];
				sorted_G._matrix[k][i] = sorted_G._matrix[k][col_number];
				sorted_G._matrix[k][col_number] = Register;
			}
		}
	}
		
	vector <int> interleaver(sorted_G.Col_number, 0);
	int j;
	/*
	// Paper
	for (size_t i(0); i < (sorted_G.Row_number - Multiple_Basis_Bits); i++) {
		interleaver.at(i) = i;
	}
	for (size_t i(sorted_G.Row_number); i < sorted_G.Col_number; i++) {
		interleaver.at(i - Multiple_Basis_Bits) = i;
	}
	for (size_t i(0); i < Multiple_Basis_Bits; i++) {
		interleaver.at(sorted_G.Col_number - (Multiple_Basis_Bits - i)) = sorted_G.Row_number - (Multiple_Basis_Bits - i);
	}
	*/
	
	// MC
	for (size_t i(0); i < (sorted_G.Row_number - Multiple_Basis_Bits); i++) {
		interleaver.at(i) = i;
	}
	for (size_t i(0); i < Multiple_Basis_Bits; i++) {
		interleaver.at(i + sorted_G.Row_number - Multiple_Basis_Bits) = i + sorted_G.Row_number;
	}
	for (size_t i(0); i < Multiple_Basis_Bits; i++) {
		interleaver.at(i + sorted_G.Row_number) = i + sorted_G.Row_number - Multiple_Basis_Bits;
	}
	for (size_t i(sorted_G.Row_number + Multiple_Basis_Bits); i < sorted_G.Col_number; i++) {
		interleaver.at(i) = i;
	}
	
	/*
	// SC
	for (size_t i(0); i < Multiple_Basis_Bits; i++) {
		interleaver.at(i) = sorted_G.Row_number - Multiple_Basis_Bits + i;
	}
	for (size_t i(Multiple_Basis_Bits); i < sorted_G.Row_number - Multiple_Basis_Bits; i++) {
		interleaver.at(i) = i;
	}
	for (size_t i(0); i < Multiple_Basis_Bits; i++) {
		interleaver.at(sorted_G.Row_number - Multiple_Basis_Bits + i) = i;
	}
	for (size_t i(sorted_G.Row_number); i < sorted_G.Col_number; i++) {
		interleaver.at(i) = i;
	}
	*/
	vector<size_t> temp_permutation_seq(permutation_seq);
	MATRIX<__int8> temp_Sorted_G(sorted_G);
	for (size_t i(0); i < sorted_G.Col_number; ++i) {
		permutation_seq.at(i) = temp_permutation_seq.at(interleaver.at(i));
		for (size_t j(0); j < sorted_G.Row_number; ++j) {
			sorted_G._matrix[j][i] = temp_Sorted_G._matrix[j][interleaver.at(i)];
		}
	}
	
	// Obtain systematic form
	GJ_Elimination(sorted_G);
	
	// Check whether the first half of matrix is identity matrix or not
	for (size_t i(0); i < sorted_G.Row_number; ++i) {
		if (sorted_G._matrix[i][i] != 1) {
			col_number = 0;

			// 同一列，由下一行開始，碰到第一個此列為 "1" 的行就交換 unsigned interleaver index
			for (size_t j(i + 1); j < sorted_G.Col_number; ++j) {
				if (sorted_G._matrix[i][j] == 1) {
					col_number = j;
					break;
				}
			}

			//Exchange the columns of permutation sequence
			Register = permutation_seq.at(i);
			permutation_seq.at(i) = permutation_seq.at(col_number);
			permutation_seq.at(col_number) = Register;

			// Exchange the columns of generator matrix
			Register = 0;
			for (size_t k(0); k < sorted_G.Row_number; ++k) {
				Register = sorted_G._matrix[k][i];
				sorted_G._matrix[k][i] = sorted_G._matrix[k][col_number];
				sorted_G._matrix[k][col_number] = Register;
			}
		}
	}
	
} // end Determining_unsigned interleaver()*/


void ParityCheck_Permutation(
	vector<double>			&ThirdParityBitLocation,
	MATRIX<__int16>			&ParityCheckbits,
	MATRIX<__int16>			&sorted_ParityCheckbits,
	vector<size_t>			&permutation_seq)
{
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;
	// Compute permutation by using "quick sorting algorithm" 
	sort(
		permutation_seq.begin(),
		permutation_seq.end(),
		[&](const size_t& a, const size_t& b)
	{return (ThirdParityBitLocation.at(a) < ThirdParityBitLocation.at(b)); });

	// Col. exchange based on permutation sequence
	
	Sort_Matrix_Row(ParityCheckbits, permutation_seq, sorted_ParityCheckbits);
}

// 將所收到的r根據前面的變化去做調整
void Sort_Function(
	vector<double>		&Rx_signal, 
	vector<size_t>		&permutation_seq, 
	vector<double>		&sorted_rx_signal)
{
	size_t counter(sorted_rx_signal.size());
	for (size_t i(0); i < counter; ++i)
		sorted_rx_signal.at(i) = Rx_signal.at(permutation_seq.at(i));
}

void Sort_Function(
	vector<__int8>		&Rx_signal,
	vector<size_t>		&permutation_seq,
	vector<__int8>		&sorted_rx_signal)
{
	size_t counter(sorted_rx_signal.size());
	for (size_t i(0); i < counter; ++i)
		sorted_rx_signal.at(i) = Rx_signal.at(permutation_seq.at(i));
}

void Sort_Parity_Length(vector<size_t>	&input_seq, vector<size_t>	&permutation_seq, vector<size_t> &output_seq) {
	for (size_t i(0); i < input_seq.size(); i++) {
		output_seq.at(i) = input_seq.at(permutation_seq.at(i));
	}
}

void Desort_Function(
	vector<size_t>		&permutation_seq, 
	vector<__int8>		&estimate_sorted_codeword,
	vector<__int8>		&estimate_codeword)
{
	size_t counter(permutation_seq.size());
	std::vector<size_t> permutation(counter,0);
	for (size_t i(0); i < counter; ++i)
		permutation.at(i) = i;
	
	// Permutation , Quick Sorting algorithm
	sort(
		permutation.begin(), 
		permutation.end(), 
		[&](const size_t& a, const size_t& b)
		{return (permutation_seq.at(a) < permutation_seq.at(b)); });

	for (size_t i(0); i < counter; ++i)
		estimate_codeword.at(i) = estimate_sorted_codeword.at(permutation.at(i));
}//end Deunsigned interleaver()

void Build_Metric_Table(vector<double> &sorted_rx_signal, MATRIX<double> &metric_table)
{
	size_t counter(sorted_rx_signal.size());
	for (size_t i(0); i < counter; ++i){
		if (sorted_rx_signal.at(i) > 0){
			metric_table._matrix[0][i] = 0;
			metric_table._matrix[1][i] = sorted_rx_signal.at(i);
		}
		else{
			metric_table._matrix[0][i] = -sorted_rx_signal.at(i);
			metric_table._matrix[1][i] = 0;
		}
	}

} //end Caculating_metric()

void Build_Fano_Metric_Table(vector<double> &sorted_rx_signal, MATRIX<double> &fano_metric_table, DECODING_INFO &decoding_info)
{
	size_t counter(sorted_rx_signal.size());
	size_t codeword_length(counter);
	size_t message_length(decoding_info.message_seq.size());
	double Es_N0; 
	Es_N0 = decoding_info.Code_Rate * pow(10, (decoding_info.SNR / 10));
	for (size_t i(0); i < message_length; ++i) {
		/*
		if (sorted_rx_signal.at(i) > 0) {
			metric_table._matrix[0][i] = 0;
			metric_table._matrix[1][i] = sorted_rx_signal.at(i);
		}
		else {
			metric_table._matrix[0][i] = -sorted_rx_signal.at(i);
			metric_table._matrix[1][i] = 0;
		}
		*/
		/*
		fano_metric_table._matrix[0][i] = log2(sqrt(Es_N0 / PI) * exp(-(Es_N0)* pow((sorted_rx_signal.at(i) - 1), 2))) + 0.5;
		fano_metric_table._matrix[1][i] = log2(sqrt(Es_N0 / PI) * exp(-(Es_N0)* pow((sorted_rx_signal.at(i) + 1), 2))) + 0.5;
		*/
		/*
		if (sorted_rx_signal.at(i) > 0) {
			fano_metric_table._matrix[0][i] = 0 - (decoding_info.Fano_Metric_Parameter * ((decoding_info.message_seq.size() - i) / (message_length / 8)));
			fano_metric_table._matrix[1][i] = sorted_rx_signal.at(i) - (decoding_info.Fano_Metric_Parameter * ((decoding_info.message_seq.size() - i) / (message_length / 8)));
		}
		else {
			fano_metric_table._matrix[0][i] = -sorted_rx_signal.at(i) - abs(decoding_info.Fano_Metric_Parameter * ((decoding_info.message_seq.size() - i) / (message_length / 8)));
			fano_metric_table._matrix[1][i] = 0 - (decoding_info.Fano_Metric_Parameter * ((decoding_info.message_seq.size() - i) / (message_length / 8)));
		}
		*/
		
		if (i < message_length) {
			if (sorted_rx_signal.at(i) > 0) {
				fano_metric_table._matrix[0][i] = 0 - decoding_info.Fano_Metric_Parameter;
				fano_metric_table._matrix[1][i] = sorted_rx_signal.at(i) - decoding_info.Fano_Metric_Parameter;
			}
			else {
				fano_metric_table._matrix[0][i] = -sorted_rx_signal.at(i) - decoding_info.Fano_Metric_Parameter;
				fano_metric_table._matrix[1][i] = 0 - decoding_info.Fano_Metric_Parameter;
			}
		}
		else {
			if (sorted_rx_signal.at(i) > 0) {
				fano_metric_table._matrix[0][i] = 0;
				fano_metric_table._matrix[1][i] = sorted_rx_signal.at(i);
			}
			else {
				fano_metric_table._matrix[0][i] = -sorted_rx_signal.at(i);
				fano_metric_table._matrix[1][i] = 0;
			}
		}
		
	}
	counter++;
}

void Place_Node(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info)//把node按照順序放入stack
{
	size_t
		position_number(Stack.size()),
		Stack_Size(Stack.size());

	for (size_t j(0); j < Stack_Size; ++j){
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).metric > child_node.metric){	//arrange stack by original metric
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

void Place_Node_ver2(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info, vector<__int8> &Record)
{
	size_t
		position_number(Stack.size()),
		Stack_Size(Stack.size()),
		constraint_i(child_node.D_z);

	if (Record.at(constraint_i) == -1) {
		cout << "?" << child_node.metric << "/" << Stack.at(Record.at(constraint_i)).metric;
		for (size_t j(Stack_Size-1); j > -1; --j) {
			++decoding_info.COM;
			// which is depend on the corresponding value of path metric to the position of new node.
			if (Stack.at(j).metric < child_node.metric) {	//arrange stack by original metric
				position_number = ++j;
				break;
			}
		}
	}
	else {
		cout << "??" << child_node.metric << "/" << Stack.at(Record.at(constraint_i)).metric;
		++decoding_info.COM;
		if (Record.at(constraint_i) == Stack_Size)position_number = Stack_Size - 1;
		else if (child_node.metric > Stack.at(Record.at(constraint_i)).metric) {
			//if()
			for (size_t j(Record.at(constraint_i)+1); j < Stack_Size; ++j) {
				++decoding_info.COM;
				// which is depend on the corresponding value of path metric to the position of new node.
				if (Stack.at(j).metric > child_node.metric) {	//arrange stack by original metric
					position_number = j;
					break;
				}
			}
		}
		else {
			for (size_t j(Record.at(constraint_i) - 1); j > -1; --j) {
				++decoding_info.COM;
				// which is depend on the corresponding value of path metric to the position of new node.
				if (Stack.at(j).metric < child_node.metric) {	//arrange stack by original metric
					position_number = ++j;
					break;
				}
			}
		}
	}
	cout << "*" << position_number << "*";
	Stack.insert(Stack.begin() + position_number, child_node);
	Record.at(constraint_i) = position_number;
	// the number of compared nodes
	//decoding_info.COM = decoding_info.COM + ((double)position_number + 1);
	//if (position_number == Stack_Size) --decoding_info.COM;

	// Stack.insert(Stack.begin() + position_number, child_node);

	// If the current capacity of stack is larger than the configured stack size
	// , then delete the last node which stores in the stack.
	if (Stack_Size > decoding_info.StackSize)
	{
		Stack.pop_back();
		if (Record.at(constraint_i) == Stack_Size) --Record.at(constraint_i);
	}
	if (decoding_info.MaxUsed_Stack < Stack_Size)
		decoding_info.MaxUsed_Stack = (unsigned int)Stack_Size;
}

void Place_Node_Fano(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info)//把node按照順序放入stack
{
	size_t
		position_number(Stack.size()),
		Stack_Size(Stack.size());

	for (size_t j(0); j < Stack_Size; ++j) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).fano_metric > child_node.fano_metric) {	//arrange stack by original metric
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

void Place_Node_multistack(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info)//把node按照順序放入multi-stack
{
	size_t
		position_number(Stack.size()),
		Stack_Size(Stack.size());

	for (size_t j(0); j < Stack_Size; ++j) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).metric > child_node.metric) {	//arrange stack by original metric
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

	if (Stack_Size > decoding_info.Multi_Stack_Size) {
		cout << "Out of Multi Stack Range ";
		system("PAUSE");
	}

}

void Pre_Procedure(
	vector<double>		&rx_signal, 
	MATRIX<__int8>		&G, 
	MATRIX<__int8>		&Sorted_G, 
	vector<size_t>		&Permutation, 
	MATRIX<double>		&metric_table)
{
	vector<double> sorting_rx_signal_seq(rx_signal);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation(rx_signal, G, Sorted_G, Permutation);
	
	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	//Sorted_R(sorting_rx_signal_seq);
}


void Pre_Procedure(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	DECODING_INFO       &decoding_info)
{
	vector<double> sorting_rx_signal_seq(rx_signal);
	//vector<double> rx_signal_seq(decoding_info.rx_signal_seq);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation(rx_signal, G, Sorted_G, Permutation);
	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	//Yin
	//vector<__int8> M = decoding_info.code_seq;
	//Sort_Function(M, Permutation, decoding_info.code_seq);
	//Yin-End
	
	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	decoding_info.Sorted_R = sorting_rx_signal_seq;
	// Yin
}

void Pre_Procedure(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	MATRIX<double>		&fano_metric_table,
	DECODING_INFO       &decoding_info)
{
	vector<double> sorting_rx_signal_seq(rx_signal);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation(rx_signal, G, Sorted_G, Permutation);

	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	//Sorted_R(sorting_rx_signal_seq);
	Build_Fano_Metric_Table(sorting_rx_signal_seq, fano_metric_table, decoding_info);
	decoding_info.Sorted_R = sorting_rx_signal_seq;
	//Yin
}

void Pre_Procedure_Fano(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&fano_metric_table,
	DECODING_INFO       &decoding_info)
{
	vector<double> sorting_rx_signal_seq(rx_signal);
	//vector<double> rx_signal_seq(decoding_info.rx_signal_seq);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation(rx_signal, G, Sorted_G, Permutation);
	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	//Yin
	//vector<__int8> M = decoding_info.code_seq;
	//Sort_Function(M, Permutation, decoding_info.code_seq);
	//Yin-End

	// To build the metric table according new received sequence.
	Build_Fano_Metric_Table(sorting_rx_signal_seq, fano_metric_table, decoding_info);
	decoding_info.Sorted_R = sorting_rx_signal_seq;
	// Yin
}

void Pre_Procedure_MultiBase(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	int					Total_Basis,
	int					Current_Basis,
	DECODING_INFO       &decoding_info)
{
	vector<double> sorting_rx_signal_seq(rx_signal);
	//vector<double> rx_signal_seq(decoding_info.rx_signal_seq);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation_MultiBase(rx_signal, G, Sorted_G, Permutation, Total_Basis,Current_Basis);
	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	//Yin
	//vector<__int8> M = decoding_info.code_seq;
	//Sort_Function(M, Permutation, decoding_info.code_seq);
	//Yin-End

	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	decoding_info.Sorted_R_Base2 = sorting_rx_signal_seq;
	// Yin
}

void Pre_Procedure_MultiBase(
	vector<double>		&rx_signal,
	MATRIX<__int8>		&G,
	MATRIX<__int8>		&Sorted_G,
	vector<size_t>		&Permutation,
	MATRIX<double>		&metric_table,
	MATRIX<double>		&fano_metric_table,
	int					Total_Basis,
	int					Current_Basis,
	DECODING_INFO       &decoding_info)
{
	vector<double> sorting_rx_signal_seq(rx_signal);
	//vector<double> rx_signal_seq(decoding_info.rx_signal_seq);

	// Determine the permutation sequence based on the received values
	// output "Sorted_G" and "Permutation"
	Determine_Permutation_MultiBase(rx_signal, G, Sorted_G, Permutation, Total_Basis, Current_Basis);
	// Reorder the received sequence according to permuation sequence
	// output new received sequence which is in order. 
	Sort_Function(rx_signal, Permutation, sorting_rx_signal_seq);
	//Yin
	//vector<__int8> M = decoding_info.code_seq;
	//Sort_Function(M, Permutation, decoding_info.code_seq);
	//Yin-End

	// To build the metric table according new received sequence.
	Build_Metric_Table(sorting_rx_signal_seq, metric_table);
	Build_Fano_Metric_Table(sorting_rx_signal_seq, fano_metric_table, decoding_info);
	decoding_info.Sorted_R_Base2 = sorting_rx_signal_seq;
	// Yin
}

void Extend_Node_Procedure(
	NODE_PATH		&Pointer,
	NODE_PATH		&child_node, 
	MATRIX<double>	&metric_table, 
	__int8			&new_bit)
{
	child_node = Pointer;
	child_node.message_bits.at(Pointer.level) = new_bit;
	++child_node.level;	
	child_node.metric += metric_table._matrix[new_bit][child_node.level - 1];
}

void Extend_Node_Procedure_Fano(
	NODE_PATH		&Pointer,
	NODE_PATH		&child_node,
	MATRIX<double>	&metric_table,
	MATRIX<double>	&fano_metric_table,
	__int8			&new_bit)
{
	child_node = Pointer;
	child_node.message_bits.at(Pointer.level) = new_bit;
	++child_node.level;
	child_node.metric += metric_table._matrix[new_bit][child_node.level - 1];
	child_node.fano_metric += fano_metric_table._matrix[new_bit][child_node.level - 1];
	//child_node.fano_metric += metric_table._matrix[new_bit][child_node.level - 1] - decoding_info.Fano_Metric_Parameter;
}


bool Update_Best_Goal_Procedure(NODE_PATH &Child_Node,NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack)
{
	if ((Child_Node.auxiliary()) < Best_Goal.metric) {
		Best_Goal = Child_Node; // Update the best candidate c odeword.
		/*unsigned int n(0), Stack_Size(Stack.size());
		while (n < Stack_Size) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) { 
				Stack.erase(Stack.begin() + n); 
				--Stack_Size;
			}
			else ++n;
		}*/

		// First, find out all indexes of deleted nodes
		size_t Stack_Size(Stack.size()), counter(0);
		vector<size_t> DeletedIndexes(Stack_Size, -1);
		for (size_t n(0); n < Stack_Size; ++n) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) {
				DeletedIndexes.at(counter) = n;
				++counter;
			}
		}
		//if (counter != 0) cout << counter << endl;
		// Second, delete those indexs from the bottom to top of stack
		while (0 < counter) {
			--counter;
			Stack.erase(Stack.begin() + DeletedIndexes.at(counter));
		}

		return TRUE;
	}
	return FALSE;
}

bool Update_Best_Goal_Procedure(NODE_PATH &Child_Node, NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack, DECODING_INFO &decoding_info)
{
	if ((Child_Node.auxiliary()) < Best_Goal.metric) {
		Best_Goal = Child_Node; // Update the best candidate c odeword.
		/*unsigned int n(0), Stack_Size(Stack.size());
		while (n < Stack_Size) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) {
				Stack.erase(Stack.begin() + n);
				--Stack_Size;
			}
			else ++n;
		}*/

		// First, find out all indexes of deleted nodes
		size_t Stack_Size(Stack.size()), counter(0);
		vector<size_t> DeletedIndexes(Stack_Size, -1);
		for (size_t n(0); n < Stack_Size; ++n) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) {
				DeletedIndexes.at(counter) = n;
				++counter;
			}
		}
		decoding_info.temp_num_Deleted_Candidate_in_Stack = 0;
		decoding_info.temp_num_Deleted_Candidate_in_Stack += counter;
		//if (counter != 0) cout << counter << endl;
		// Second, delete those indexs from the bottom to top of stack
		while (0 < counter) {
			--counter;
			Stack.erase(Stack.begin() + DeletedIndexes.at(counter));
		}

		return TRUE;
	}
	return FALSE;
}

//PoHan
bool Update_Best_Goal_Procedure(NODE_PATH &Child_Node, NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack, vector<NODE_PATH> &Multi_Stack)
{
	if ((Child_Node.auxiliary()) < Best_Goal.metric) {
		Best_Goal = Child_Node; // Update the best candidate c odeword.
		/*unsigned int n(0), Stack_Size(Stack.size());
		while (n < Stack_Size) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) {
				Stack.erase(Stack.begin() + n);
				--Stack_Size;
			}
			else ++n;
		}*/

		// First, find out all indexes of deleted nodes
		size_t Stack_Size(Stack.size()), Multi_Stack_Size(Multi_Stack.size()), counter(0);
		vector<size_t> DeletedIndexes(Stack_Size, -1), DeletedIndexes_Multi_Stack(Multi_Stack_Size, -1);
		for (size_t n(0); n < Stack_Size; ++n) {
			if (Stack.at(n).auxiliary() >= Best_Goal.metric) {
				DeletedIndexes.at(counter) = n;
				++counter;
			}
		}
		//if (counter != 0) cout << counter << endl;
		// Second, delete those indexs from the bottom to top of stack
		while (0 < counter) {
			--counter;
			Stack.erase(Stack.begin() + DeletedIndexes.at(counter));
		}

		counter = 0;
		for (size_t n(0); n < Multi_Stack_Size; ++n) {
			if (Multi_Stack.at(n).auxiliary() >= Best_Goal.metric) {
				DeletedIndexes_Multi_Stack.at(counter) = n;
				++counter;
			}
		}
		//if (counter != 0) cout << counter << endl;
		// Second, delete those indexs from the bottom to top of stack
		while (0 < counter) {
			--counter;
			Multi_Stack.erase(Multi_Stack.begin() + DeletedIndexes_Multi_Stack.at(counter));
		}

		return TRUE;
	}
	return FALSE;
}

//For A* Parity 
void Q_Function(MATRIX<__int8> &G, MATRIX<size_t>& Q_table)
{
	// scan new generator matrix 
	// Compute at least how many MRIP bits is need to uniquely determine each 
	vector<size_t> Code_Bit_Position(G.Col_number - G.Row_number, 0);
	
	/*
		Code_Bit_Position 用來儲存至少需要前幾個 bits 才能決定該位置上的 parity bits
		EX:
		G =	[1 0 0  1 1 0]
			[0 1 0  1 0 1]
			[0 0 1  0 1 1]
		->	       [1 2 2]

		位置 0 的 parity bit 至少需要 msg bit 0、1　  才能被決定
		     1					      msg bit 0、1、2
		     2					      msg bit 0、1、2
	*/
	for (size_t i(G.Row_number); i < G.Col_number; ++i)
		for (size_t j(0); j < G.Row_number; ++j)
			if (G._matrix[j][i] == 1) Code_Bit_Position.at(i - G.Row_number) = j;
	
	
	
	// Example 
	// Q_table : size : K by (N-K) 
	// 列所引值 X 代表目前已知道 bit 0 到 bit X (Q function input)
	// X 列中所存的非零元素就是已知 bit 0 ~ bit X，則能被決定的 parity bit 編號				
	//				known MRIP bits		 determined parity bits  	
	// [0 0 0 0]	bit 0			->		NULL
	// [0 0 0 0]	bit 0,1			->		NULL
	// [4 0 0 0]	bit 0,1,2		->		bit 4
	// [5 6 7 0]	bit 0,1,2,3		->		bit 4,5,6,7
	
	size_t x(0); // 存放目前的位置，為了往後增長
	for (size_t j(0); j < G.Row_number; ++j){
		x = 0;
		for (size_t i(0); i < G.Col_number - G.Row_number; ++i)
			if (Code_Bit_Position.at(i) == j){
				Q_table._matrix[j][x] = i + G.Row_number;
				++x;
			}
	}
}// End Q_Function()

// For Parity-A* decoding algorithm 
// used for encoding a part of message bits for obtaining the specific parity bit 
__int8 Bit_Encoder(MATRIX<__int8> &G, size_t parity_index, NODE_PATH &node)
{
	__int8 x(0);
	for (size_t i(0); i < node.level; ++i) {
		x ^= ((node.message_bits.at(i)) & G._matrix[i][parity_index]);
	}
	return x;
}

void Place_Node_f(vector<NODE_PATH> &Stack, NODE_PATH &child_node, DECODING_INFO &decoding_info, double alpha)
{
	size_t position_number = Stack.size();
	for (size_t j(0); j < Stack.size(); ++j) {
		// use auxiliary metric to sort stack elements
		if (Stack.at(j).auxiliary() > child_node.auxiliary()) {
			position_number = j;
			break;
		}
	}
	decoding_info.COM = decoding_info.COM + ((double)position_number + 1);
	if (position_number == Stack.size()) --decoding_info.COM;

	Stack.insert(Stack.begin() + position_number, child_node);

	if (Stack.size() > decoding_info.StackSize) {
		// 找位於後 (alpha * S) 個位置中擁有最小 level 的 node，將他們存入此 vector 中 
		vector<size_t> Delete_Index(1, (size_t)(decoding_info.StackSize*alpha));

		size_t temp_level = (decoding_info.StackSize*(1 - alpha));
		for (size_t i(decoding_info.StackSize*(1 - alpha) + 1); i < decoding_info.StackSize; ++i) {
			if (Stack.at(i).level == temp_level)
				Delete_Index.push_back(i);
			else if (Stack.at(i).level < temp_level) {
				Delete_Index.clear();
				temp_level = Stack.at(i).level;
				Delete_Index.push_back(i);
			}
		}

		// 找出 Delete_Index 當中擁有最大 metric 的 node，將之刪除
		if (Delete_Index.size() != 1) {
			size_t temp_index = Delete_Index.at(0);
			double temp_auxiliary = Stack.at(Delete_Index.at(0)).auxiliary();
			for (size_t i(0); i < Delete_Index.size(); i++) {
				if (Stack.at(Delete_Index.at(i)).auxiliary() > temp_auxiliary) {
					temp_auxiliary = Stack.at(Delete_Index.at(i)).auxiliary();
					temp_index = Delete_Index.at(i);
				}
			}
			Stack.erase(Stack.begin() + temp_index);
		}
		else {
			Stack.erase(Stack.begin() + Delete_Index.at(0));
		}
	}

	if (decoding_info.MaxUsed_Stack < Stack.size())
		decoding_info.MaxUsed_Stack = (unsigned int)Stack.size();
}

void Hard_decision_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		if (decoding_info.rx_signal_seq.at(i) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
}
//Pin
void Hard_decision_test(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		D_z1 = 0,D_z2 = 0;

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		Sorted_codeword(codeword_length, 0);


	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		if (decoding_info.rx_signal_seq.at(i) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	Sort_Function(decoding_info.estimated_codeword, Location_Index, codeword_seq);
	//計算每個位子上的錯誤
	/*
	for (int i = 0; i < codeword_length; i++) {
		if (codeword_seq.at(i) != Sorted_codeword.at(i)) {
			decoding_info.err_count[i] += 1;
		}
	}*/
	for (int i = 0; i < message_length /2; i++) {
		if (codeword_seq.at(i) != Sorted_codeword.at(i)) {
			D_z1 += 1;
		}
	}
	for (int i = message_length / 2; i < message_length; i++) {
		if (codeword_seq.at(i) != Sorted_codeword.at(i)) {
			D_z2 += 1;
		}
	}
	decoding_info.err_count[D_z1] += 1;
	decoding_info.err_count[message_length / 2 + D_z2] += 1;

}

void A_star_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	Best_Goal.metric = FLT_MAX;

	vector<NODE_PATH> Stack;
	Stack.reserve(decoding_info.StackSize + 1);
	Stack.push_back(Pointer);

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	do{
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

void A_star_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack;
	Stack.push_back(Pointer);

	Best_Goal.metric = FLT_MAX;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	//unsigned int node_counter = 0;
	do
	{
		//node_counter = 0;
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		/*if (Pounsigned inter.level == (message_length - 1)) 
			Stack.erase(Stack.begin());*/
	
		for (__int8 new_bit(0); new_bit < 2; ++new_bit){
			++decoding_info.STE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			
			/*compute heuristic value
			
			Q({u0,u1,...uX}) = { P } -> Q(X) = { P } ， { P } 為已知 u0 ~ uX 可被決定的 parity bits 所成集合
			Q_Table[x][]  where x 代表目前 parital message 有 bit 0 ~ bit x  (message bits)

			1.依照 Q_Table 中儲存的索引值
			2.將 partial message 進行編碼，得到此 bit 值
			3.將對應此所引值及其 bit 值的 bit metric 提出
			4.加總所有的 Q_Table中的所有所引值得到的 bit metric 即為此 node 的 heuristic value

			metric[ 目前 bit 值 ][ 目前 bit 位置 ]
			x -> 0 or 1
			*/

			size_t y(0);
			while (Q_table._matrix[Pointer.level][y] != 0){
				size_t x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
				++y;
				++decoding_info.STE;
			}

			if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric)){
				//
				decoding_info.CandidateCodeWord++;
				//
				Child_Node.metric = Child_Node.auxiliary();
				Child_Node.heuristic = 0;
				//
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric)) {
				
				/*if (Pointer.metric != Child_Node.metric) {
					Place_Node(Stack, Child_Node, decoding_info);
				}
				else {
					Stack.at(0) = Child_Node;
					decoding_info.COM++;
				}*/
				Place_Node(Stack, Child_Node, decoding_info);
			}
			/*else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() >= Best_Goal.metric)) {
				node_counter++;
				if (node_counter == 2) {
					Stack.erase(Stack.begin());
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

void A_star_Parity_f(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); 

	do
	{
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			decoding_info.STE++;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);

			//Compute Heuristic Value
			size_t y(0);
			while (Q_table._matrix[Pointer.level][y] != 0) {
				size_t x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];
				++y;
				++decoding_info.STE;
			}

			if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric)) {
				decoding_info.CandidateCodeWord++;
				
				Child_Node.metric = Child_Node.auxiliary();
				Child_Node.heuristic = 0;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric)) {
				Place_Node_f(Stack, Child_Node, decoding_info, decoding_info.Alpha);
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

void A_star_Fano_metric(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
	MATRIX<double> Fano_Metric_Table(2, codeword_length);


	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	Best_Goal.metric = FLT_MAX;

	vector<NODE_PATH> Stack(1, Pointer);
	decoding_info.Code_Rate = (double)message_length / (double)codeword_length;
	
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);

	do {
		Pointer = Stack.at(0);

		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure_Fano(Pointer, Child_Node, Metric_Table, Fano_Metric_Table, new_bit);
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
				Place_Node_Fano(Stack, Child_Node, decoding_info);
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

void swap_PoHan(NODE_PATH &a, NODE_PATH &b) {
	NODE_PATH temp = a;
	a = b;
	b = temp;
}
int Partition(vector<NODE_PATH> &Stack_CBC1, int front, int end, DECODING_INFO &decoding_info) {
	double pivot = Stack_CBC1.at(end).fano_metric;
	int i = front - 1;
	for (int j = front; j < end; j++) {
		if (Stack_CBC1.at(j).fano_metric < pivot) {
			i++;
			swap_PoHan(Stack_CBC1.at(i), Stack_CBC1.at(j));
		}
		decoding_info.COM++;
	}
	i++;
	swap_PoHan(Stack_CBC1.at(i), Stack_CBC1.at(end));
	return i;
}
void QuickSort(vector<NODE_PATH> &Stack_CBC1, int front, int end, DECODING_INFO &decoding_info) {
	if (front < end) {
		int pivot = Partition(Stack_CBC1, front, end, decoding_info);
		QuickSort(Stack_CBC1, front, pivot - 1, decoding_info);
		QuickSort(Stack_CBC1, pivot + 1, end, decoding_info);
	}
}
