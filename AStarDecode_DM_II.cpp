#include "PinAstar/AStarDecode.h"

//Deleting machanism for checking level
/****************************************************************************/
bool Deleting_Mechanisim(MATRIX<__int8> &G, NODE_PATH &Check_Node, vector<__int8> &Ref_sequence, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position;
	vector <__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Check_Node);
	Best_Goal.D_z = codeword_length;

	do {
		Pointer = Stack.at(0);
		
		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			
			Child_Node = Pointer;
			++Child_Node.level;
			Child_Node.message_bits.at(Pointer.level) = new_bit;
			if (new_bit != Ref_sequence.at(Pointer.level)) ++Child_Node.D_z;
			++decoding_info.DM_STE;

			if ((Child_Node.level == message_length) && (Child_Node.D_z < Best_Goal.D_z) && (Child_Node.D_z <= decoding_info.Constraint_i)) {

				// DM Statistic
				decoding_info.DM_STE += (decoding_info.Control_Level - message_length);
				decoding_info.CandidateCodeWord++;

				// Generate Candidate Codeword
				message_seq = Child_Node.message_bits;
				Systematic_Linear_Block_Code_Encoder(G, message_seq, codeword_seq);
				for (size_t j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Ref_sequence.at(j)) {
						Child_Node.D_z++;
					};
				}

				// Delete Bad Nodes
				if (Child_Node.D_z < Best_Goal.D_z) {
					Best_Goal = Child_Node;
					if (!Stack.empty()) {
						while ((Stack.at(Stack.size() - 1).D_z >= Best_Goal.D_z)) {
							Stack.pop_back();
							if (Stack.empty())  break;
						}
					}
				}
			}
			else if ((Child_Node.level < message_length) && (Child_Node.D_z < Best_Goal.D_z) && (Child_Node.D_z <= decoding_info.Constraint_i)) {

				if (Child_Node.D_z != Pointer.D_z) {
					position = (size_t)Stack.size();
					for (size_t j(0); j < Stack.size(); j++) {
						decoding_info.DM_COM++;
						if (Stack.at(j).D_z >= Child_Node.D_z) {
							position = (size_t)j;
							break;
						}
					}
					Stack.insert(Stack.begin() + position, Child_Node);
				}
				else {
					Stack.at(0) = Child_Node;
					decoding_info.DM_COM++;
				}
				if (Stack.size() > decoding_info.DM_II_StackSize) Stack.pop_back();
			}
		}
	} while (!Stack.empty());

	if (Best_Goal.D_z <= decoding_info.Constraint_j) return FALSE;
	else return TRUE;
}

void A_star_PC_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		H_Seq(codeword_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	
	// Hard decision
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) H_Seq.at(i) = 1;
	}
	
	bool delete_flag = FALSE;
	bool first_path_flag = FALSE;
	do {
		Pointer = Stack.at(0);
		Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) && (first_path_flag)) {
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, H_Seq, decoding_info);
		}


		if (!delete_flag) {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit){
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != H_Seq[Pointer.level]) Child_Node.D_z++;
				decoding_info.STE++;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					first_path_flag = TRUE;
					
					decoding_info.STE += (codeword_length - message_length);
					decoding_info.CandidateCodeWord++;

					message_seq = Child_Node.message_bits;
					Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, codeword_seq);
					for (size_t j(message_length); j < codeword_length; ++j)
						Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					Place_Node(Stack, Child_Node, decoding_info);
				}
			}
		}
		else {
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
		}
	} while (!Stack.empty());

	//
	message_seq = Best_Goal.message_bits;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	
	//
	decoding_info.STE /= (double)message_length;
	decoding_info.COM /= (double)message_length;

	//
	decoding_info.DM_STE /= (double)message_length;
	decoding_info.DM_COM /= (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_PCout_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number);

	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		H_Seq(codeword_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);


	// Hard decision
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) H_Seq.at(i) = 0;
		else H_Seq.at(i) = 1;
	}

	bool delete_flag = FALSE;
	bool first_path_flag = FALSE; 
	do {
		Pointer = Stack.at(0);
		//Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) && (first_path_flag)){ 
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, H_Seq, decoding_info);
		}
		else if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		if (!delete_flag) {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				//
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != H_Seq.at(Pointer.level)) Child_Node.D_z++;

				decoding_info.STE++;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					//
					decoding_info.STE += (codeword_length - message_length);
					decoding_info.CandidateCodeWord++;
					first_path_flag = TRUE;
					//
					Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
					for (size_t j(message_length); j < codeword_length; ++j)
						Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
					//
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i))
				{
					for (size_t j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits.at(j) = H_Seq.at(j);
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

				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
					if(Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						decoding_info.COM++;
					}
				}
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
		}
	} while (!Stack.empty());


	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE /= (double)message_length;
	decoding_info.COM /= (double)message_length;

	//
	decoding_info.DM_STE /= (double)message_length;
	decoding_info.DM_COM /= (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_BMA(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
	
	// Hard decision
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX[i] = 0;
		else Hard_RX[i] = 1;
	}
	//*/

	bool delete_flag = FALSE;
	bool first_path_flag = FALSE; // for the first completed path, we do not checking procedure
	do{
		Pointer = Stack[0];
		Stack.erase(Stack.begin());
		
		//Checking procedure
		if ((Pointer.level == decoding_info.Check_Level) && first_path_flag) {
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag) {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) Child_Node.D_z++;

				decoding_info.STE++;

				if ((Child_Node.level == message_length) && (Child_Node.metric <= Best_Goal.metric)) {
					first_path_flag = TRUE;
					decoding_info.STE += (codeword_length - message_length);
					decoding_info.CandidateCodeWord++;

					Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
					for (size_t j(message_length); j < codeword_length; ++j)
						Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];

					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric <= Best_Goal.metric)) {
					Place_Node(Stack, Child_Node, decoding_info);
				}
			}
		}
		else {
			decoding_info.DeletedNode_counter++;
			delete_flag = FALSE;
		}

	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	//
	decoding_info.DM_STE = decoding_info.DM_STE / (double)message_length;
	decoding_info.DM_COM = decoding_info.DM_COM / (double)message_length;
	//decoding_info.deleted_node_counter = decoding_info.deleted_node_counter / (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
}

void A_star_1_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;
	//bool first_path_flag = FALSE;
	do {
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) ){
			delete_flag = Deleting_Mechanisim(G, Pointer, Hard_RX, decoding_info);
		}
		else if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		if (!delete_flag) {
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX[Pointer.level]) ++Child_Node.D_z;
				++decoding_info.STE;

				if ((Child_Node.level == message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z <= pc_i)) {
					//
					decoding_info.STE += (codeword_length - message_length);
					decoding_info.CandidateCodeWord++;
					//
					Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);
					for (size_t j(message_length); j < codeword_length; ++j)
						Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
					//
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z < pc_i)) {
					if(Pointer.metric != Child_Node.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack[0] = Child_Node;
						decoding_info.COM++;
					}

				}
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z == pc_i)) {
					
					for (size_t j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits[j] = Hard_RX[j];
						decoding_info.STE++;
					}
					//
					decoding_info.STE += (codeword_length - message_length);
					decoding_info.CandidateCodeWord++;
					//
					Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);
					for (size_t j(message_length); j < codeword_length; ++j)
						Child_Node.metric += Metric_Table._matrix[codeword_seq[j]][j];
					//
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;

	do
	{
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Checki procedure
		if (Pointer.level == decoding_info.Check_Level) {
			delete_flag = Deleting_Mechanisim(G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag) {
			if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
				//
				Stack.erase(Stack.begin());
				//
				Systematic_Linear_Block_Code_Encoder(G, Pointer.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
				//
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
					if ((Child_Node.metric < Best_Goal.metric))
					{
						if(Pointer.metric !=Child_Node.metric)	
							Place_Node(Stack, Child_Node, decoding_info);
						else {
							Stack[0] = Child_Node;
							decoding_info.COM++;
						}
					}
				}
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
				Stack.erase(Stack.begin());
				//
				Best_Goal_2 = Best_Goal;
				A_star_1_Stack_DM_II(G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
		}

	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_3_Stack_DM_II(MATRIX<__int8> &G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;

	do
	{
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Check procedure
		if (Pointer.level == decoding_info.Check_Level) {
			delete_flag = Deleting_Mechanisim(G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag) {
			if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
				//
				Stack.erase(Stack.begin());
				//
				Systematic_Linear_Block_Code_Encoder(G, Pointer.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
				//
				decoding_info.STE += (codeword_length - message_length);
				decoding_info.CandidateCodeWord++;
				//
				Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z < pc_i)) {

				decoding_info.STE += 2;
				for (__int8 new_bit(0); new_bit < 2; new_bit++) {

					Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
					if (new_bit != Hard_RX[Pointer.level])
						Child_Node.D_z = Pointer.D_z + 1;
					if ((Child_Node.metric < Best_Goal.metric)) {
						if (Pointer.metric != Child_Node.metric)
							Place_Node(Stack, Child_Node, decoding_info);
						else {
							Stack[0] = Child_Node;
							decoding_info.COM++;
						}
					}		
				}
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z == pc_i)) {
				//
				Stack.erase(Stack.begin());
				//
				Best_Goal_2 = Best_Goal;
				A_star_2_Stack_DM_II(G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
				//
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
		}

	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_2_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;
	bool first_path_flag = FALSE;
	
	// It is very important because Multi-Stack haven't set this parameter. 
	// Therefore, the situation will cause error in using Deleting_Mechanism
	
	decoding_info.Constraint_i = 3;

	// for OSD-i
	for (size_t i(0); i < codeword_length; i++) {
		if (Metric_Table._matrix[0][i] == 0) Hard_RX[i] = 0;
		else Hard_RX[i] = 1;
	}

	do {
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) && (first_path_flag)) {
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag){
			if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
				Stack.erase(Stack.begin());
				
				first_path_flag = TRUE;
				
				// calculate the compeleted path metric 
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];

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

					if ((Child_Node.metric < Best_Goal.metric)) {
						if (Pointer.metric != Child_Node.metric)
							Place_Node(Stack, Child_Node, decoding_info);
						else {
							Stack[0] = Child_Node;
							decoding_info.COM++;
						}
					}
				}
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
				Stack.erase(Stack.begin());
				//
				Best_Goal_2 = Best_Goal;
				A_star_1_Stack_DM_II(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, (size_t)3, decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
		}
		else{
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
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

void A_star_3_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;
	bool first_path_flag = FALSE;

	// It is very important because Multi-Stack haven't set this parameter. 
	// Therefore, the situation will cause error in using Deleting_Mechanism
	decoding_info.Constraint_i = 4;

	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
	}

	do {
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) && (first_path_flag)) {
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag) {
			if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
				first_path_flag = TRUE;
				//
				Stack.erase(Stack.begin());
				 
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
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
					if (new_bit != Hard_RX[Pointer.level])
						Child_Node.D_z = Pointer.D_z + 1;

					if ((Child_Node.metric < Best_Goal.metric)) {
						if(Pointer.metric != Child_Node.metric)
							Place_Node(Stack, Child_Node, decoding_info);
						else {
							Stack[0] = Child_Node;
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
				A_star_2_Stack_DM_II(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
				//
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
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

void A_star_4_Stack_DM_II(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
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
	bool delete_flag = FALSE;
	bool first_path_flag = FALSE;

	// It is very important because Multi-Stack haven't set this parameter. 
	// Therefore, the situation will cause error in using Deleting_Mechanism

	decoding_info.Constraint_i = 5;

	// for OSD-i
	for (size_t i( 0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
	}

	do {
		Pointer = Stack[0];
		//Stack.erase(Stack.begin());

		//Check procedure
		if ((Pointer.level == decoding_info.Check_Level) && (first_path_flag)) {
			delete_flag = Deleting_Mechanisim(Sorted_G, Pointer, Hard_RX, decoding_info);
		}

		if (!delete_flag) {
			if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
				Stack.erase(Stack.begin());
				first_path_flag = TRUE;
				
				Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);
				for (size_t j(message_length); j < codeword_length; ++j)
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];

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

					if ((Child_Node.metric < Best_Goal.metric)) {
						if (Pointer.metric != Child_Node.metric)
							Place_Node(Stack, Child_Node, decoding_info);
						else {
							Stack[0] = Child_Node;
							decoding_info.COM++;
						}
					}
					
				}
			}
			else if ((Pointer.level < message_length) && (Pointer.D_z == 1)) {
				Stack.erase(Stack.begin());
				Best_Goal_2 = Best_Goal;
				A_star_3_Stack_DM_II(Sorted_G, Metric_Table, Hard_RX, Pointer, Best_Goal_2, 2, decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
		}
		else {
			Stack.erase(Stack.begin());
			delete_flag = FALSE;
			decoding_info.DeletedNode_counter++;
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
