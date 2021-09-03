#include "PinAstar/AStarDecode.h"

//Deleting machanism for parity bits
/****************************************************************************/
void ParityPath_Checking(
	MATRIX<__int8>			&G,
	DECODING_INFO			&decoding_info, 
	MATRIX<double>			&Metric_Table, 
	NODE_PATH				&Child_Node, 
	NODE_PATH				&Best_Goal, 
	vector<NODE_PATH>		&Stack, 
	vector <__int8>			&Hard_RX) {

	vector <__int8>codeword_seq(G.Col_number, 0);
	Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);

	// Coputation for the number of error in control band
	for (size_t j(G.Row_number); j < decoding_info.Control_Level; ++j) {
		if (codeword_seq.at(j) != Hard_RX.at(j)) {
			++Child_Node.D_zp;
			if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
				break;
			}
		}
	}
	Child_Node.heuristic = 0;
	// Satisfy the constraint number j so this goal node's metric can be calculated
	if ((Child_Node.D_z + Child_Node.D_zp) <= decoding_info.Constraint_j) {
		decoding_info.STE += (G.Col_number - G.Row_number);
		++decoding_info.CandidateCodeWord;
		//
		for (size_t i(G.Row_number); i < G.Col_number; ++i) {
			if (codeword_seq.at(i) != Hard_RX.at(i)) {
				Child_Node.metric += Metric_Table._matrix[codeword_seq.at(i)][i];
				if (Child_Node.metric > Best_Goal.metric) break;
			}
		}
		if(Child_Node.metric < Best_Goal.metric)
			Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
	}
	else {
		++decoding_info.DeletedNode_counter;
	} 
}

void A_star_PC_out_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack.at(0);

		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z++;
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) && 
				(Child_Node.metric < Best_Goal.metric) && 
				(Child_Node.D_z <= decoding_info.Constraint_i)) {

				ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);

			}
			else if (
				(Child_Node.level < message_length) && 
				(Child_Node.metric < Best_Goal.metric) && 
				(Child_Node.D_z == decoding_info.Constraint_i)) {

				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
					++decoding_info.STE;
				}

				ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);
			}
			else if (
				(Child_Node.level < message_length) && 
				(Child_Node.metric < Best_Goal.metric) &&
				(Child_Node.D_z < decoding_info.Constraint_i)) {
				
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.DM_COM;
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

void A_star_1_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack.at(0);

		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
			++decoding_info.STE;

			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= pc_i)) {
				// ++decoding_info.CandidateCodeWord;
				ParityPath_Checking(G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < pc_i)) {
				if (Pointer.metric != Child_Node.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == pc_i)) {

				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
					++decoding_info.STE;
				}
				// ++decoding_info.CandidateCodeWord;
				ParityPath_Checking(G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);
			}
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {
			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;

				if ((Child_Node.metric < Best_Goal.metric)){
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
			A_star_1_stack_DM_I(G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_3_stack_DM_I(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {
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
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_2_stack_DM_I(G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}

	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
	for (size_t i(0); i < codeword_length; ++i) 
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) Child_Node.D_z = Pointer.D_z + 1;  

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
			//
			Stack.erase(Stack.begin());
			//
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_DM_I(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 3, decoding_info);
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

void A_star_3_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
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

void A_star_4_stack_DM_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		else Hard_RX.at(0) = 1;
	}

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			Stack.erase(Stack.begin());
			//
			ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Pointer, Best_Goal, Stack, Hard_RX);
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

//Deleting machanism-I for parity bits + heuristic metric of parity bits
/****************************************************************************/
void A_star_PC_out_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
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

	//*** for OSD-i, Hard Decision of received sequence
	for (size_t i(0); i < codeword_length; ++i) 
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
	
	//***

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table);

	do {
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			++decoding_info.STE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);

			//*** Compute Hamming Distance : Dz = d(c,z)
			if (new_bit != Hard_RX.at(Pointer.level)) ++Child_Node.D_z;
			//***

			//*** Compute heuristic values and the Hamming distance within control level
			unsigned int y = 0, x = 0;
			while ((Q_table._matrix[Pointer.level][y] != 0)) {
				x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

				if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
					if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

					if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
						++decoding_info.DeletedNode_counter;
						break;
					}
				}
				++y;
				++decoding_info.STE;
			}
			//***

			if ((Child_Node.D_z + Child_Node.D_zp) <= decoding_info.Constraint_j) {
				if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					++decoding_info.CandidateCodeWord;
					Child_Node.metric = Child_Node.auxiliary();
					Child_Node.heuristic = 0;
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {

					for (unsigned int j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits[j] = Hard_RX[j];
						decoding_info.STE++;
					}

					Child_Node.D_zp = 0;
					Child_Node.heuristic = 0;
					ParityPath_Checking(Sorted_G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
					Place_Node(Stack, Child_Node, decoding_info);
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

void A_star_1_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table,  MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			
			//*** Compute Hamming Distance : Dz = d(c,z)
			if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z++;
			//***
			
			//*** Compute heuristic values and the Hamming distance within control level
			size_t y = 0;
			__int8 x = 0;
			while ((Q_table._matrix[Pointer.level][y] != 0)) {
				x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
				Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

				if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
					if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

					if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
						++decoding_info.DeletedNode_counter;
						break;
					}
				}
				++y;
				++decoding_info.STE;
			}

			if ((Child_Node.D_z + Child_Node.D_zp) <= decoding_info.Constraint_j) {
				if ((Child_Node.level == message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z <= pc_i)) {
					++decoding_info.CandidateCodeWord;
					Child_Node.metric = Child_Node.auxiliary();
					Child_Node.heuristic = 0;
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z < pc_i)) {
					Place_Node(Stack, Child_Node, decoding_info);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.auxiliary() < Best_Goal.metric) && (Child_Node.D_z == pc_i)) {

					for (size_t j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits[j] = Hard_RX[j];
						decoding_info.STE++;
					}
					Child_Node.D_zp = 0;
					Child_Node.heuristic = 0;
					ParityPath_Checking(G, decoding_info, Metric_Table, Child_Node, Best_Goal, Stack, Hard_RX);
				}
			}
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			++decoding_info.CandidateCodeWord;
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0;
				__int8 x = 0;
				while ((Q_table._matrix[Pointer.level][y] != 0)) {
					x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

					if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
						if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

						if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
							++decoding_info.DeletedNode_counter;
							break;
						}
					}
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) <= decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_DM_I_Parity(G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_3_stack_DM_I_Parity(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, MATRIX<size_t> &Q_table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			++decoding_info.CandidateCodeWord;
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0, x = 0;
				while ((Q_table._matrix[Pointer.level][y] != 0)) {
					x = Bit_Encoder(G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

					if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
						if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

						if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
							++decoding_info.DeletedNode_counter;
							break;
						}
					}
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) < decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
			Best_Goal_2 = Best_Goal;
			A_star_2_stack_DM_I_Parity(G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		if (Metric_Table._matrix[0][i] == 0) Hard_RX[i] = 0;
		else Hard_RX[i] = 1;
	}

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			++decoding_info.CandidateCodeWord;
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {
				
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;
				
				unsigned int y = 0, x = 0;
				while ((Q_table._matrix[Pointer.level][y] != 0)) {
					x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

					if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
						if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

						if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
							++decoding_info.DeletedNode_counter;
							break;
						}
					}
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) < decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_DM_I_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 3, decoding_info);
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

void A_star_3_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來
	
	// OSC
	double metric_thr(0);
	for (size_t i(0); i < message_length; i++) {
		//metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_thr = 0.1*metric_thr;

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
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;

				size_t y(0), x(0);
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

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) < decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			//
			Best_Goal_2 = Best_Goal;
			//
			A_star_2_stack_DM_I_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
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

void A_star_4_stack_DM_I_Parity(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		if (Metric_Table._matrix[0][i] == 0) Hard_RX[i] = 0;
		else Hard_RX[i] = 1;
	}

	MATRIX<size_t> Q_table(message_length, (codeword_length - message_length));
	Q_Function(Sorted_G, Q_table); // Q_table 儲存 Q_function 的值，一次算完儲存起來

	do {
		Pointer = Stack[0];
		Stack.erase(Stack.begin());

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			++decoding_info.CandidateCodeWord;
			Pointer.metric = Pointer.auxiliary();
			Pointer.heuristic = 0;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z = Pointer.D_z + 1;

				size_t y = 0, x = 0;
				while ((Q_table._matrix[Pointer.level][y] != 0)) {
					x = Bit_Encoder(Sorted_G, Q_table._matrix[Pointer.level][y], Child_Node);
					Child_Node.heuristic += Metric_Table._matrix[x][Q_table._matrix[Pointer.level][y]];

					if (Q_table._matrix[Pointer.level][y] <= decoding_info.Control_Level) {
						if (x != Hard_RX[Q_table._matrix[Pointer.level][y]]) Child_Node.D_zp++;

						if ((Child_Node.D_z + Child_Node.D_zp) > decoding_info.Constraint_j) {
							++decoding_info.DeletedNode_counter;
							break;
						}
					}
					++y;
					++decoding_info.STE;
				}

				if ((Child_Node.auxiliary() < Best_Goal.metric) && ((Child_Node.D_z + Child_Node.D_zp) < decoding_info.Constraint_j))
					Place_Node(Stack, Child_Node, decoding_info);
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
			Best_Goal_2 = Best_Goal;
			A_star_3_stack_DM_I_Parity(Sorted_G, Metric_Table, Q_table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
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

