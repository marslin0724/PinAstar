#include "AStarDecode.h"
MATRIX<double> * M;
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
	vector<NODE_PATH> S_Stack1,S_Stack2;	//C_Stack為存放達到segment_length的node
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
	//先讓C_Stack1與C_Stack2有node
	/*
	Pre_Procedure_Segment(decoding_info, Metric_Table,S_Stack1,C_Stack1, segment_length,0);
	Pre_Procedure_Segment(decoding_info, Metric_Table,S_Stack2,C_Stack2, segment_length, segment_length);
	//結合兩個node成一個candidate codeword
	Combine_Segment(C_Stack1.at(0), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
		G, Sorted_G, Metric_Table);
	if(C_Stack1.size() >= 2)
		Combine_Segment(C_Stack1.at(1), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
	*/
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

void A_star_Segment_ver2(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{
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
		
	M = &Metric_Table;

	Best_Goal.metric = FLT_MAX;

	vector<NODE_PATH> S_Stack1, S_Stack2, C_Stack1, C_Stack2;	//C_Stack為存放達到segment_length的node
	forward_list<NODE_PATH> O_Stack;					//存放結合後的codeword
	S_Stack1.reserve(decoding_info.StackSize + 1);
	S_Stack2.reserve(decoding_info.StackSize + 1);
	C_Stack1.reserve(decoding_info.StackSize + 1);
	C_Stack2.reserve(decoding_info.StackSize + 1);
	//O_Stack.reserve(decoding_info.StackSize + 1);		改成用deque
	S_Stack1.push_back(Pointer);
	Pointer.level = segment_length;
	S_Stack2.push_back(Pointer);
	//選擇Stack的參考
	size_t stack_flag = 2;
	
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	//先讓C_Stack1與C_Stack2有node
	Pre_Procedure_Segment(decoding_info, Metric_Table, S_Stack1, C_Stack1, segment_length, 0);
	Pre_Procedure_Segment(decoding_info, Metric_Table, S_Stack2, C_Stack2, segment_length, segment_length);
	//結合兩個node成一個candidate codeword
	
	Combine_Segment(C_Stack1.at(0), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
		G, Sorted_G, Metric_Table);
	
	if (C_Stack1.size() >= 2)
		Combine_Segment(C_Stack1.at(1), C_Stack1, C_Stack2, S_Stack1, S_Stack2, Best_Goal, segment_length, decoding_info,
			G, Sorted_G, Metric_Table);

	/*
	for (int i = 0; i < 4;i++) {
		cout << (int)Best_Goal.message_bits.at(i) << ",";
	}
	cout << endl;*/
	
	while ((!S_Stack1.empty()) || (!S_Stack2.empty()) || (!O_Stack.empty())) {
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
		else if(!S_Stack2.empty()){
			Pointer = S_Stack2.at(0);
			stack_flag = 2;
		}
		else {
			stack_flag = 3;
		}
		if ((!O_Stack.empty()) && ((!S_Stack2.empty()) || (!S_Stack1.empty())) ) {
			if (O_Stack.front().metric < Pointer.metric)
				stack_flag = 3;
		}
		if (stack_flag == 1) {
			if (Pointer.level == (segment_length - 1))
				S_Stack1.erase(S_Stack1.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				++decoding_info.STE;

				if ((Child_Node.level == segment_length) && (Child_Node.metric < Best_Goal.metric)) {
					Place_C_Stack(C_Stack1, Child_Node, decoding_info);
					Combine_Segment_ver2(Child_Node, C_Stack2, O_Stack,Best_Goal, segment_length, decoding_info,
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
		else if(stack_flag == 2) {
			if (Pointer.level == (message_length - 1))
				S_Stack2.erase(S_Stack2.begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				++decoding_info.STE;
				
				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric)) {
					Place_C_Stack(C_Stack2, Child_Node, decoding_info);
					Combine_Segment_ver2(Child_Node, C_Stack1, O_Stack,Best_Goal, segment_length, decoding_info,
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
		else {
			if (O_Stack.front().metric < Best_Goal.metric) {
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Systematic_Linear_Block_Code_Encoder(Sorted_G, O_Stack.front().message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					O_Stack.front().metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				//
				if (Update_Best_Goal_Segments(O_Stack.front(), Best_Goal, C_Stack1)) {
					Update_Stack(Best_Goal, C_Stack2);
					Update_Stack(Best_Goal, S_Stack1);
					Update_Stack(Best_Goal, S_Stack2);
				}
				if (O_Stack.front().segment == 2) next_extend(O_Stack.front(), C_Stack1, O_Stack, Best_Goal, segment_length, message_length, decoding_info);
				else next_extend(O_Stack.front(), C_Stack2, O_Stack, Best_Goal, segment_length, message_length, decoding_info);
			}
			else O_Stack.clear();
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
			Update_Stack( Best_Goal, C_Stack_Update);
			Update_Stack( Best_Goal, S_Stack1);
			Update_Stack( Best_Goal, S_Stack2);
			//更新bound
			metric_bound = Best_Goal.metric - Update_node.metric;
		}
		

	}
}
//改版 C_Stacl_Update使用NODE_COMB
inline void Combine_Segment(NODE_PATH &Update_node, vector<NODE_COMB>& C_Stack_Update, vector<NODE_COMB> &C_Stack, vector<NODE_PATH>& S_Stack1, vector<NODE_PATH>& S_Stack2,
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
			Update_Stack(Best_Goal, S_Stack1);
			Update_Stack(Best_Goal, S_Stack2);
			//更新bound
			metric_bound = Best_Goal.metric - Update_node.metric;
			
		}


	}
}

inline void Combine_Segment_ver2(NODE_PATH &Update_node, vector<NODE_PATH> &C_Stack, forward_list<NODE_PATH>& O_Stack,
	NODE_PATH &Best_Goal, size_t segment_length,
	DECODING_INFO& decoding_info, MATRIX<__int8>& G, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		start, end;
	NODE_PATH
		Combine_node(message_length);
	Combine_node = Update_node;
	//test(Update_node, Metric_Table);
	//Combine_node.message_bits.assign(Update_node.message_bits.begin(), Update_node.message_bits.end());
	if (Update_node.level == message_length) {
		start = 0, end = segment_length;
		Combine_node.segment = 2;
	}
	else {
		start = segment_length, end = message_length;
		Combine_node.segment = 1;
	}
	Combine_node.idx_end = C_Stack.size();
	for (int j = start; j < end; j++) {
		Combine_node.message_bits.at(j) = C_Stack.at(0).message_bits.at(j);
	}

	Combine_node.metric = Update_node.metric + C_Stack.at(0).metric;
	Combine_node.segment_metric = Update_node.metric;
	/*
	for (int i = 0; i < 4; i++) {
		cout << (int)Combine_node.message_bits.at(i) << ",";
	}
	cout << "metric = " << (double)Combine_node.metric << endl;
	*/

	Place_O_Stack(O_Stack, Combine_node, decoding_info);
	
}



inline void Place_C_Stack(vector<NODE_PATH>& Stack, NODE_PATH& child_node, DECODING_INFO& decoding_info)//把node按照順序放入stack
{
	int
		position_number(0),
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
//改版
inline void Place_C_Stack(vector<NODE_COMB>& Stack, NODE_COMB& child_node, DECODING_INFO& decoding_info)//把node按照順序放入stack
{
	int
		position_number(0),
		Stack_Size(Stack.size());

	for (int j(Stack_Size - 1); j >= 0; --j) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if (Stack.at(j).metric < child_node.metric) {	//arrange stack by original metric
			position_number = j + 1;
			break;
		}
	}

	// the number of compared nodes
	decoding_info.COM = decoding_info.COM + ((double)Stack_Size - position_number + 1);
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

inline void Place_O_Stack(forward_list<NODE_PATH>& Stack, NODE_PATH& child_node, DECODING_INFO& decoding_info)//把node按照順序放入stack
{
	forward_list<NODE_PATH>::iterator pre, cur = Stack.begin();
	int count = 0;
	if (Stack.empty() || (*cur).metric > child_node.metric) {
		Stack.push_front(child_node);
		return;
	}
	pre = cur;
	cur++;
	for (;cur != Stack.end() ; ++cur) {
		// which is depend on the corresponding value of path metric to the position of new node.
		if ((*cur).metric > child_node.metric) {	//arrange stack by original metric
			Stack.insert_after(pre, child_node);
			decoding_info.COM = decoding_info.COM + (count + 1);
			return;
		}
		count++;
		pre = cur;
	}
	Stack.insert_after(pre, child_node);
	
	decoding_info.COM = decoding_info.COM + (count + 1);
}

inline void Update_Stack(NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack)
{
		// First, find out all indexes of deleted nodes
		int Stack_Size(Stack.size());
		int n;
		for ( n = (Stack_Size - 1); n >= 0; --n) {
			if (Stack.at(n).auxiliary() <= Best_Goal.metric) {
				break;
			}
		}
		Stack.erase(Stack.begin() + n + 1,Stack.end());
}
inline void Update_Stack(NODE_PATH &Best_Goal, vector<NODE_COMB> &Stack)
{
	// First, find out all indexes of deleted nodes
	int Stack_Size(Stack.size());
	int n;
	for (n = (Stack_Size - 1); n >= 0; --n) {
		if (Stack.at(n).metric <= Best_Goal.metric) {
			break;
		}
	}
	Stack.erase(Stack.begin() + n + 1, Stack.end());
}
// polymorphism
inline void Update_Stack(NODE_PATH &Best_Goal, deque<NODE_PATH> &Stack)
{
	// First, find out all indexes of deleted nodes
	int Stack_Size(Stack.size());
	int n;
	for (n = (Stack_Size - 1); n > 0; --n) {	//不刪除 stack.at(0)
		if (Stack.at(n).auxiliary() <= Best_Goal.metric) {
			break;
		}
	}
	Stack.erase(Stack.begin() + n + 1, Stack.end());
}

inline bool Update_Best_Goal_Segments(NODE_PATH &Child_Node, NODE_PATH &Best_Goal, vector<NODE_PATH> &Stack) {
	if ((Child_Node.auxiliary()) < Best_Goal.metric) {
		Best_Goal = Child_Node; 

		// First, find out all indexes of deleted nodes
		int Stack_Size(Stack.size());
		int n;
		for (n = (Stack_Size - 1); n >= 0; --n) {
			if (Stack.at(n).auxiliary() <= Best_Goal.metric) {
				break;
			}
		}
		Stack.erase(Stack.begin() + n + 1, Stack.end());
		return TRUE;
	}
	return FALSE;
}

inline void next_extend(NODE_PATH& Node, vector<NODE_PATH> &C_Stack, forward_list<NODE_PATH>& O_Stack, NODE_PATH &Best_Goal,int segment_length,int message_length,DECODING_INFO& decoding_info) {
	Node.temp_idx++;
	int i = 1,start,end;
	forward_list<NODE_PATH>::iterator cur = O_Stack.begin(),pre;
	if (Node.temp_idx < Node.idx_end && Node.temp_idx < C_Stack.size()) {
		Node.metric = Node.segment_metric + C_Stack.at(Node.temp_idx).metric;
		if (Node.metric < Best_Goal.metric) {
			if (Node.segment == 2) {
				start = 0, end = segment_length;
			}
			else {
				start = segment_length, end = message_length;
			}
			for (int j = start; j < end; j++) {
				Node.message_bits.at(j) = C_Stack.at(Node.temp_idx).message_bits.at(j);
			}
			pre = cur;
			cur++;
			if (cur->metric > Node.metric) {
				return;
			}
			pre = cur;
			cur++;
			while(cur != O_Stack.end()) {
				i++;
				if (cur->metric > Node.metric) {
					break;
				}
				if (cur->metric == Node.metric) {
					O_Stack.pop_front();
					decoding_info.COM += i;
					return;
				}
				
				
				pre = cur;
				cur++;
			}
			O_Stack.insert_after(pre, Node);
			O_Stack.pop_front();
			decoding_info.COM += i;
			return;
		}
	}
	O_Stack.pop_front();
	
}
inline bool meg_equal(vector<__int8> meg1, vector<__int8> meg2,int len) {
	for (int i = 0; i < len; i++) {
		if (meg1.at(i) != meg2.at(i)) return false;
	}
	return true;
}

inline void test(NODE_PATH& node, MATRIX<double>& Metric_Table) {
	double metric = 0;
	static int error_count;
	
	if (node.level == 10) {
		for (int i = 5; i < 10; i++) {
			metric += Metric_Table._matrix[node.message_bits.at(i)][i];
		}
	}
	else {
		for (int i = 0; i < 5; i++) {
			metric += Metric_Table._matrix[node.message_bits.at(i)][i];
		}
	}
	/*
	for (int i = 0; i < 10; i++) {
		metric += Metric_Table._matrix[node.message_bits.at(i)][i];
	}
	*/
	if (node.metric != metric) {
		cout << "metric error,count = " << ++error_count << endl;
		for (auto c : node.message_bits) {
			cout << (int)c << ",";
		}
		cout << endl;
	}
}

//分成N段
void A_star_N_Section(MATRIX<__int8>& G, DECODING_INFO& decoding_info)
{	
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		segment_length(message_length / decoding_info.N_section);
	double
		min_metric;

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
	
	vector<vector<NODE_PATH>> S_Stack;	//C_Stack為存放達到segment_length的node
	vector<vector<NODE_COMB>> C_Stack;
	S_Stack.resize(decoding_info.N_section);
	C_Stack.resize(decoding_info.N_section);
	for (int i = 0; i < decoding_info.N_section; i++) {
		S_Stack.at(i).reserve(decoding_info.StackSize + 1);
		C_Stack.at(i).reserve(decoding_info.StackSize + 1);
		Pointer.level = segment_length * i;
		S_Stack.at(i).push_back(Pointer);
	}
	//選擇Stack的參考
	size_t stack_flag ;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	while (MutiStack_nonEmpty(S_Stack)) {
		min_metric = INT_MAX;
		for (int i = decoding_info.N_section-1; i >= 0; i--) {
			if (!S_Stack.at(i).empty()) {
				if (S_Stack.at(i).at(0).metric < min_metric) {
					stack_flag = i;
					min_metric = S_Stack.at(i).at(0).metric;
				}
			}
		}
		decoding_info.COM++;
		Pointer = S_Stack.at(stack_flag).at(0);
		if (stack_flag != decoding_info.N_section - 1) {
			if (Pointer.level == (segment_length*(stack_flag + 1) - 1))
				S_Stack.at(stack_flag).erase(S_Stack.at(stack_flag).begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				++decoding_info.STE;

				if ((Child_Node.level == segment_length * (stack_flag + 1) )&& (Child_Node.metric < Best_Goal.metric)) {
					N_Combine_Segment(Child_Node, C_Stack, S_Stack, stack_flag, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
				}
				else if ((Child_Node.level < segment_length * (stack_flag + 1)) && (Child_Node.metric < Best_Goal.metric)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack.at(stack_flag), Child_Node, decoding_info);
					else {
						S_Stack.at(stack_flag).at(0) = Child_Node;
						++decoding_info.COM;
					}
				}
			}
		}
		else {
			if ((Pointer.level == message_length - 1))
				S_Stack.at(stack_flag).erase(S_Stack.at(stack_flag).begin());

			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				++decoding_info.STE;

				if ((Child_Node.level == message_length)&& (Child_Node.metric < Best_Goal.metric)) {
					N_Combine_Segment(Child_Node, C_Stack, S_Stack, stack_flag, Best_Goal, segment_length, decoding_info,
						G, Sorted_G, Metric_Table);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric)) {
					if (Child_Node.metric != Pointer.metric)
						Place_Node(S_Stack.at(stack_flag), Child_Node, decoding_info);
					else {
						S_Stack.at(stack_flag).at(0) = Child_Node;
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

inline bool MutiStack_nonEmpty(vector<vector<NODE_PATH>>& S_Stack) {
	bool flag = 0;
	for (int i = 0; i < S_Stack.size(); i++) {
		flag |= !S_Stack.at(i).empty();
	}
	return flag;
}

//改版 C_Stacl_Update使用NODE_COMB N個Section
inline void N_Combine_Segment(NODE_PATH &Update_node, vector<vector<NODE_COMB>> &C_Stack, vector<vector<NODE_PATH>>& S_Stack,size_t& stack_flag,
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
	//double metric_bound = Best_Goal.metric - Update_node.metric;
	Combine_node.codeword_bits.assign(Update_node.message_bits.begin(), Update_node.message_bits.end());
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Combine_node.codeword_bits, Combine_node.codeword_bits);
	decoding_info.STE += (codeword_length - message_length);
	Combine_node.metric = Update_node.metric;
	Place_C_Stack(C_Stack.at(stack_flag), Combine_node, decoding_info);
	metric = Best_Goal.metric;
	recursive_combine(Combine_node, C_Stack, 0, stack_flag,
		Combine_node.metric, Best_Goal, Sorted_G, Metric_Table);
	if (metric > Best_Goal.metric) {
		for (int i = 0; i < C_Stack.size(); i++) {
			Update_Stack(Best_Goal, C_Stack.at(i));
			Update_Stack(Best_Goal, S_Stack.at(i));
		}
	}	
}

void recursive_combine(NODE_COMB Combine_node, vector<vector<NODE_COMB>> &C_Stack,size_t level,size_t& stack_flag,
	double sum_metric, NODE_PATH &Best_Goal, MATRIX<__int8>& Sorted_G, MATRIX<double>& Metric_Table) {
	NODE_COMB temp_node = Combine_node;
	if (level == C_Stack.size()) {
		for (size_t j(Sorted_G.Row_number); j < Sorted_G.Col_number; ++j)
			sum_metric += Metric_Table._matrix[Combine_node.codeword_bits.at(j)][j];
		if (sum_metric < Best_Goal.metric) {
			for (int j = 0; j < Sorted_G.Row_number; j++) {
				Best_Goal.message_bits.at(j) = Combine_node.codeword_bits.at(j);
			}
			Best_Goal.metric = sum_metric;
		}
		return;
	}
	if (level != stack_flag) {
		for (int i = 0; i < C_Stack.at(level).size(); i++) {
			if (sum_metric + C_Stack.at(level).at(i).metric < Best_Goal.metric) {
				for (int k = 0; k < Combine_node.codeword_bits.size(); k++) {
					temp_node.codeword_bits.at(k) =
						 Combine_node.codeword_bits.at(k) ^ C_Stack.at(level).at(i).codeword_bits.at(k);
					
				}
				recursive_combine(temp_node, C_Stack, level + 1, stack_flag,
					sum_metric + C_Stack.at(level).at(i).metric, Best_Goal, Sorted_G, Metric_Table);
			}
			else break;
		}
	}
	else {
		recursive_combine(temp_node, C_Stack, level + 1, stack_flag,
			sum_metric , Best_Goal, Sorted_G, Metric_Table);
	}
}
