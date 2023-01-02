#include "PinAstar/AStarDecode.h"
#include <iostream>
#include <vector>

using namespace std;

NODE_PATH Control_Band_Check_1bit(
	MATRIX <__int8>				&Sorted_G, //Sorted_G
	MATRIX <double>				&Metric_Table, //Metric_Table
	vector <__int8>				&Input_codeword_seq, //codeword_seq
	vector <__int8>				&Hard_RX, //Hard_RX
	NODE_PATH					&Input_Node, //Child_Node
	NODE_PATH					&Best_Goal_Node, //Best_Goal
	DECODING_INFO				&Decoding_info) //decoding_info
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8> codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	NODE_PATH Temp_Child_Node(message_length);
	double metric(0);
	int ADD_STE = 0;

	for (size_t bit_index(Input_Node.level); bit_index < message_length; ++bit_index) {
		Decoding_info.CBC_STE += (++ADD_STE);
		++Decoding_info.CBC_Candidate;
		++Decoding_info.Counter;
		Temp_Child_Node.D_z = Input_Node.D_z;
		Temp_Child_Node.Same_Index = Input_Node.Same_Index;
		codeword_seq = Input_codeword_seq;
		codeword_seq.at(bit_index) ^= 1;		// one-bit flipping in MRIP
		Temp_Child_Node.D_z++;
		error_counter = Input_Node.D_z + 1;		
		metric = Input_Node.metric + Metric_Table._matrix[Hard_RX.at(bit_index) ^ 1][bit_index];
		++Decoding_info.Counter;
		

		//reach control level for checking the number of errors
		for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
			//cout << bit_index << endl;
			if (Sorted_G._matrix[bit_index][j] == 1) 
				codeword_seq.at(j) ^= 1;
			
			if (codeword_seq.at(j) != Hard_RX.at(j)) { 
				++error_counter; 
				if (error_counter > Decoding_info.Constraint_j){
					++Decoding_info.DM_STE;
					break;
				}
				metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Temp_Child_Node.D_z++;
				//if (metric > best_codeword.metric) break;
			}
			else {
				Temp_Child_Node.Same_Index.push_back(j);
			}
		}
		Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);
		// codewords
		// compute the corresponding metric of codewords
		if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
			//the constraint j is satisfied, the procedure goes to calculate the rest parity bits and its metrics
			for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
				if (Sorted_G._matrix[bit_index][j] == 1) 
					codeword_seq.at(j) ^= 1;
				
				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					Temp_Child_Node.D_z++;
					if (metric > best_codeword.metric) break;
				}
				else {
					Temp_Child_Node.Same_Index.push_back(j);
				}
			}

			Decoding_info.STE += (codeword_length - message_length + 1);
			
			++Decoding_info.CandidateCodeWord;

			if (metric < best_codeword.metric) {
				best_codeword.metric = metric;
				best_codeword.D_z = Temp_Child_Node.D_z;
				best_codeword.Same_Index = Temp_Child_Node.Same_Index;
				//best_codeword.D_z = Decoding_info.Constraint_i + 1;
				for (size_t j(0); j < message_length; ++j) {
					best_codeword.message_bits.at(j) = codeword_seq.at(j);
				}
			}
			
		}
	}
	return best_codeword;
}

NODE_PATH Control_Band_Check_1bit(
	MATRIX <__int8>				&Sorted_G, //Sorted_G
	MATRIX <double>				&Metric_Table, //Metric_Table
	vector <__int8>				&Input_codeword_seq, //codeword_seq
	vector <__int8>				&Hard_RX, //Hard_RX
	NODE_PATH					&Input_Node, //Child_Node
	NODE_PATH					&Best_Goal_Node, //Best_Goal
	double						OSC_metric_thr, // OSC Metric Threshold
	double						Worst_metric, // Worst Metric
	DECODING_INFO				&Decoding_info) //decoding_info
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8> codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	double metric(0);
	int ADD_STE = 0;

	for (size_t bit_index(Input_Node.level); bit_index < message_length; ++bit_index) {
		Decoding_info.CBC_STE += (++ADD_STE);
		++Decoding_info.CBC_Candidate;
		++Decoding_info.Counter;
		codeword_seq = Input_codeword_seq;
		codeword_seq.at(bit_index) ^= 1;		// one-bit flipping in MRIP
		error_counter = Input_Node.D_z + 1;
		metric = Input_Node.metric + Metric_Table._matrix[Hard_RX.at(bit_index) ^ 1][bit_index];
		++Decoding_info.Counter;

		//reach control level for checking the number of errors
		for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
			//cout << bit_index << endl;
			if (Sorted_G._matrix[bit_index][j] == 1)
				codeword_seq.at(j) ^= 1;

			if (codeword_seq.at(j) != Hard_RX.at(j)) {
				++error_counter;
				if (error_counter > Decoding_info.Constraint_j) {
					++Decoding_info.DM_STE;
					break;
				}
				metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				if (metric > best_codeword.metric) break;
			}
		}
		Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);
		// codewords
		// compute the corresponding metric of codewords
		if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
			//the constraint j is satisfied, the procedure goes to calculate the rest parity bits and its metrics
			for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
				if (Sorted_G._matrix[bit_index][j] == 1)
					codeword_seq.at(j) ^= 1;

				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (metric > best_codeword.metric) break;
				}
			}

			Decoding_info.STE += (codeword_length - message_length + 1);
			++Decoding_info.CandidateCodeWord;

			if (metric < best_codeword.metric) {
				best_codeword.metric = metric;
				best_codeword.D_z = Decoding_info.Constraint_i + 1;
				for (size_t j(0); j < message_length; ++j) {
					best_codeword.message_bits.at(j) = codeword_seq.at(j);
				}
			}
			if (best_codeword.metric < OSC_metric_thr) break;
		}
	}
	return best_codeword;
}

// flip two bits
NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
	int                          CBC1)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8>
		codeword_seq(Input_codeword_seq),
		best_codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	double metric(0);

	size_t bit2_index;

	for (size_t bit1_index(Input_Node.level); bit1_index < message_length; ++bit1_index) {
		//cout << bit1_index << ": ";
		if (CBC1 != -1 && CBC1 > bit1_index) {
			bit2_index = CBC1;
		}
		else bit2_index = bit1_index;
		
		
		/*
		if (Decoding_info.CBC1 != -1 && Decoding_info.CBC1 > bit1_index) {
			bit2_index = Decoding_info.CBC1;
		}
		else bit2_index = bit1_index;*/
		for (bit2_index; bit2_index < message_length; ++bit2_index) {
			++Decoding_info.Counter;
			//cout << bit2_index << " ";
	// flip bits starting from the last
	//for (size_t bit1_index(message_length - 1); bit1_index >= Input_Node.level; --bit1_index) {
	//	for (size_t bit2_index(message_length - 1); bit2_index >= bit1_index; --bit2_index) {
			
	        // generate message bits
			codeword_seq = Input_codeword_seq;
			if (bit1_index != bit2_index) {  // bit flipping CBC-2
				codeword_seq.at(bit1_index) ^= 1;
				codeword_seq.at(bit2_index) ^= 1;
				error_counter = Input_Node.D_z + 2;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
					+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index];
			}
			else { // CBC-1
				codeword_seq.at(bit1_index) ^= 1;
				error_counter = Input_Node.D_z + 1;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index];
			}

			// Encode / Deleting mechanism: 算到Control level的長度  (k--Control Level)
			for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {	
				if (bit1_index != bit2_index) {
					if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
						codeword_seq.at(j) ^= 1;
				}
				else if (bit1_index == bit2_index) {
					if (Sorted_G._matrix[bit1_index][j] == 1)
						codeword_seq.at(j) ^= 1;
				}

				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					++error_counter;
					if (error_counter > Decoding_info.Constraint_j) break;

					metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (metric > best_codeword.metric) break;
				}
			}

			Decoding_info.Binary_STE += (Decoding_info.Control_Level - Input_Node.level);

			// Check the number of errors within control level   (Control Level--n)
			if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
				for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
					if (bit1_index != bit2_index) {
						if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
							codeword_seq.at(j) ^= 1;
					}
					else if (bit1_index == bit2_index) {
						if (Sorted_G._matrix[bit1_index][j] == 1)
							codeword_seq.at(j) ^= 1;
					}
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (metric > best_codeword.metric) break;
					}
				}

				Decoding_info.STE += (codeword_length - message_length + 2); // for 2 bit flipping
				if (bit1_index == bit2_index)--Decoding_info.STE; // for 1 bit flipping
				++Decoding_info.CandidateCodeWord;

				if (metric < best_codeword.metric) {
					best_codeword.metric = metric;
					best_codeword_seq = codeword_seq;
				}
			}
		}
	}
	//cout << "END" << endl << endl;
	best_codeword.message_bits.assign(
		best_codeword_seq.begin(), 
		best_codeword_seq.begin() + message_length);
	return best_codeword;
}

NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8>
		codeword_seq(Input_codeword_seq),
		best_codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	NODE_PATH Temp_Child_Node(message_length);
	double metric(0);

	int ADD_STE = 0,ADD_STE_CYCLE2;
	for (size_t bit1_index(Input_Node.level); bit1_index < message_length; ++bit1_index) {
		if (bit1_index < message_length / 2) bit1_index = message_length / 2;
		Decoding_info.CBC_STE += (++ADD_STE);
		ADD_STE_CYCLE2 = ADD_STE;
		while (--ADD_STE_CYCLE2 > 0) {
			Decoding_info.CBC_STE += ADD_STE_CYCLE2;
		}
		//cout << bit1_index << ": ";
		for (size_t bit2_index(bit1_index); bit2_index < message_length; ++bit2_index) {
			//if (Decoding_info.Cancelled_Candidate_i != 0 && )
			//Decoding_info.CBC_STE += ADD_STE;
			++Decoding_info.CBC_Candidate;
			++Decoding_info.Counter;
			Temp_Child_Node.D_z = Input_Node.D_z;
			Temp_Child_Node.Same_Index = Input_Node.Same_Index;
			//cout << bit2_index << " ";
	// flip bits starting from the last
	//for (size_t bit1_index(message_length - 1); bit1_index >= Input_Node.level; --bit1_index) {
	//	for (size_t bit2_index(message_length - 1); bit2_index >= bit1_index; --bit2_index) {

			// generate message bits
			codeword_seq = Input_codeword_seq;
			if (bit1_index != bit2_index) {  // bit flipping CBC-2
				codeword_seq.at(bit1_index) ^= 1;
				codeword_seq.at(bit2_index) ^= 1;
				error_counter = Input_Node.D_z + 2;
				Temp_Child_Node.D_z += 2;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
					+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index];
			}
			else { // CBC-1
				codeword_seq.at(bit1_index) ^= 1;
				error_counter = Input_Node.D_z + 1;
				Temp_Child_Node.D_z += 1;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index];
			}

			// Encode / Deleting mechanism: 算到Control level的長度  (k--Control Level)
			for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
				if (bit1_index != bit2_index) {
					if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
						codeword_seq.at(j) ^= 1;
				}
				else if (bit1_index == bit2_index) {
					if (Sorted_G._matrix[bit1_index][j] == 1)
						codeword_seq.at(j) ^= 1;
				}

				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					++error_counter;
					if (error_counter > Decoding_info.Constraint_j) {
						++Decoding_info.DM_STE;
						break;
					}
					metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					Temp_Child_Node.D_z++;
					if (metric > best_codeword.metric) break;
				}
				else {
					Temp_Child_Node.Same_Index.push_back(j);
				}

			}

			Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);

			// Check the number of errors within control level   (Control Level--n)
			if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
				for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
					if (bit1_index != bit2_index) {
						if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
							codeword_seq.at(j) ^= 1;
					}
					else if (bit1_index == bit2_index) {
						if (Sorted_G._matrix[bit1_index][j] == 1)
							codeword_seq.at(j) ^= 1;
					}
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						Temp_Child_Node.D_z++;
						if (metric > best_codeword.metric) break;
					}
					else {
						Temp_Child_Node.Same_Index.push_back(j);
					}
				}

				Decoding_info.STE += (codeword_length - message_length + 2); // for 2 bit flipping
				if (bit1_index == bit2_index)--Decoding_info.STE; // for 1 bit flipping
				++Decoding_info.CandidateCodeWord;

				if (metric < best_codeword.metric) {
					best_codeword.metric = metric;
					best_codeword_seq = codeword_seq;
					best_codeword.D_z = Temp_Child_Node.D_z;
					best_codeword.Same_Index = Temp_Child_Node.Same_Index;
					//Yin
					/*
					if (bit1_index != bit2_index) best_codeword.D_z = Decoding_info.Constraint_i + 2;
					else best_codeword.D_z = Decoding_info.Constraint_i + 1;
					*/
				}
			}
		}
	}
	//cout << "END" << endl << endl;
	best_codeword.message_bits.assign(
		best_codeword_seq.begin(),
		best_codeword_seq.begin() + message_length);
	return best_codeword;
}

NODE_PATH Control_Band_Check_2bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	double						OSC_metric_thr,
	DECODING_INFO				&Decoding_info)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8>
		codeword_seq(Input_codeword_seq),
		best_codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	double metric(0);

	int ADD_STE = 0, ADD_STE_CYCLE2;
	for (size_t bit1_index(Input_Node.level); bit1_index < message_length; ++bit1_index) {
		if (bit1_index < message_length / 2) bit1_index = message_length / 2;
		Decoding_info.CBC_STE += (++ADD_STE);
		ADD_STE_CYCLE2 = ADD_STE;
		while (--ADD_STE_CYCLE2 > 0) {
			Decoding_info.CBC_STE += ADD_STE_CYCLE2;
		}
		//cout << bit1_index << ": ";
		for (size_t bit2_index(bit1_index); bit2_index < message_length; ++bit2_index) {
			//if (Decoding_info.Cancelled_Candidate_i != 0 && )
			//Decoding_info.CBC_STE += ADD_STE;
			++Decoding_info.CBC_Candidate;
			++Decoding_info.Counter;
			//cout << bit2_index << " ";
	// flip bits starting from the last
	//for (size_t bit1_index(message_length - 1); bit1_index >= Input_Node.level; --bit1_index) {
	//	for (size_t bit2_index(message_length - 1); bit2_index >= bit1_index; --bit2_index) {

			// generate message bits
			codeword_seq = Input_codeword_seq;
			if (bit1_index != bit2_index) {  // bit flipping CBC-2
				codeword_seq.at(bit1_index) ^= 1;
				codeword_seq.at(bit2_index) ^= 1;
				error_counter = Input_Node.D_z + 2;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
					+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index];
			}
			else { // CBC-1
				codeword_seq.at(bit1_index) ^= 1;
				error_counter = Input_Node.D_z + 1;
				metric =
					Input_Node.metric
					+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index];
			}

			// Encode / Deleting mechanism: 算到Control level的長度  (k--Control Level)
			for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
				if (bit1_index != bit2_index) {
					if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
						codeword_seq.at(j) ^= 1;
				}
				else if (bit1_index == bit2_index) {
					if (Sorted_G._matrix[bit1_index][j] == 1)
						codeword_seq.at(j) ^= 1;
				}

				if (codeword_seq.at(j) != Hard_RX.at(j)) {
					++error_counter;
					if (error_counter > Decoding_info.Constraint_j) {
						++Decoding_info.DM_STE;
						break;
					}
					metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (metric > best_codeword.metric) break;
				}
			}

			Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);

			// Check the number of errors within control level   (Control Level--n)
			if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
				for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
					if (bit1_index != bit2_index) {
						if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
							codeword_seq.at(j) ^= 1;
					}
					else if (bit1_index == bit2_index) {
						if (Sorted_G._matrix[bit1_index][j] == 1)
							codeword_seq.at(j) ^= 1;
					}
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (metric > best_codeword.metric) break;
					}
				}

				Decoding_info.STE += (codeword_length - message_length + 2); // for 2 bit flipping
				if (bit1_index == bit2_index)--Decoding_info.STE; // for 1 bit flipping
				++Decoding_info.CandidateCodeWord;

				if (metric < best_codeword.metric) {
					best_codeword.metric = metric;
					best_codeword_seq = codeword_seq;
					//Yin
					if (bit1_index != bit2_index) best_codeword.D_z = Decoding_info.Constraint_i + 2;
					else best_codeword.D_z = Decoding_info.Constraint_i + 1;

				}
				if (best_codeword.metric < OSC_metric_thr)break;
			}
		}
		if (best_codeword.metric < OSC_metric_thr)break;
	}
	//cout << "END" << endl << endl;
	best_codeword.message_bits.assign(
		best_codeword_seq.begin(),
		best_codeword_seq.begin() + message_length);
	return best_codeword;
}

// flip two bits, limited CBC
NODE_PATH Control_Band_Check_2bits_Limited(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
	size_t                       limited_CBC_length)
	{
		size_t
			message_length(Sorted_G.Row_number),
			codeword_length(Sorted_G.Col_number),
			error_counter(0);
		vector<__int8>
			codeword_seq(Input_codeword_seq),
			best_codeword_seq(Input_codeword_seq);
		NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
		double metric(0);
		
		if (Input_Node.level < limited_CBC_length) {
			//cout << Input_Node.level << "," << limited_CBC_length << endl;
			Input_Node.level = limited_CBC_length;
		}
		for (size_t bit1_index(Input_Node.level); bit1_index < message_length; ++bit1_index) {
			//cout << bit1_index << ": ";
			for (size_t bit2_index(bit1_index); bit2_index < message_length; ++bit2_index) {
				++Decoding_info.Counter;
				//cout << bit2_index << " ";
		// flip bits starting from the last
		//for (size_t bit1_index(message_length - 1); bit1_index >= Input_Node.level; --bit1_index) {
		//	for (size_t bit2_index(message_length - 1); bit2_index >= bit1_index; --bit2_index) {

				// generate message bits
				codeword_seq = Input_codeword_seq;
				if (bit1_index != bit2_index) {  // bit flipping CBC-2
					codeword_seq.at(bit1_index) ^= 1;
					codeword_seq.at(bit2_index) ^= 1;
					error_counter = Input_Node.D_z + 2;
					metric =
						Input_Node.metric
						+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
						+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index];
				}
				else { // CBC-1
					codeword_seq.at(bit1_index) ^= 1;
					error_counter = Input_Node.D_z + 1;
					metric =
						Input_Node.metric
						+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index];
				}

				// Encode / Deleting mechanism: 算到Control level的長度  (k--Control Level)
				for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
					if (bit1_index != bit2_index) {
						if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
							codeword_seq.at(j) ^= 1;
					}
					else if (bit1_index == bit2_index) {
						if (Sorted_G._matrix[bit1_index][j] == 1)
							codeword_seq.at(j) ^= 1;
					}

					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						++error_counter;
						if (error_counter > Decoding_info.Constraint_j) break;

						metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (metric > best_codeword.metric) break;
					}
				}

				Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);

				// Check the number of errors within control level   (Control Level--n)
				if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
					for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
						if (bit1_index != bit2_index) {
							if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
								codeword_seq.at(j) ^= 1;
						}
						else if (bit1_index == bit2_index) {
							if (Sorted_G._matrix[bit1_index][j] == 1)
								codeword_seq.at(j) ^= 1;
						}
						if (codeword_seq.at(j) != Hard_RX.at(j)) {
							metric += Metric_Table._matrix[codeword_seq.at(j)][j];
							if (metric > best_codeword.metric) break;
						}
					}

					Decoding_info.STE += (codeword_length - message_length + 2); // for 2 bit flipping
					if (bit1_index == bit2_index)--Decoding_info.STE; // for 1 bit flipping					
					++Decoding_info.CandidateCodeWord;

					if (metric < best_codeword.metric) {
						best_codeword.metric = metric;
						best_codeword_seq = codeword_seq;
					}
				}
			}
		}
		//cout << "END" << endl << endl;
		best_codeword.message_bits.assign(
			best_codeword_seq.begin(),
			best_codeword_seq.begin() + message_length);
		return best_codeword;
	}

// flip three bits
NODE_PATH Control_Band_Check_3bits(
	MATRIX <__int8>				&Sorted_G,
	MATRIX <double>				&Metric_Table,
	vector <__int8>				&Input_codeword_seq,
	vector <__int8>				&Hard_RX,
	NODE_PATH					&Input_Node,
	NODE_PATH					&Best_Goal_Node,
	DECODING_INFO				&Decoding_info,
	double                       limited_ratio)
{
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0);
	vector<__int8>
		codeword_seq(Input_codeword_seq),
		best_codeword_seq(Input_codeword_seq);
	NODE_PATH best_codeword((size_t)message_length); best_codeword.metric = Best_Goal_Node.metric;
	double metric(0);

	if (limited_ratio != 1) {
		Input_Node.level = limited_ratio * message_length;
	}
	int ADD_STE = 0, ADD_STE_CYCLE2, ADD_STE_CYCLE3;
	//cout << "Start";
	for (size_t bit1_index(Input_Node.level); bit1_index < message_length; ++bit1_index) {
		if (bit1_index < message_length / 2) bit1_index = message_length / 2;
		Decoding_info.CBC_STE += (++ADD_STE);
		ADD_STE_CYCLE2 = ADD_STE;
		while (--ADD_STE_CYCLE2 > 0) {
			Decoding_info.CBC_STE += ADD_STE_CYCLE2;
			ADD_STE_CYCLE3 = ADD_STE_CYCLE2;
			while (--ADD_STE_CYCLE3 > 0) {
				Decoding_info.CBC_STE += ADD_STE_CYCLE3;
			}
		}
		for (size_t bit2_index(bit1_index); bit2_index < message_length; ++bit2_index) {
			for (size_t bit3_index(bit2_index); bit3_index < message_length; ++bit3_index) {
				// flip bits starting from the last
				if ((bit3_index != bit2_index) && (bit1_index == bit2_index)) continue;
				++Decoding_info.Counter;
				++Decoding_info.CBC_Candidate;
				Decoding_info.CBC_STE += ADD_STE;
				//cout << "(" << bit1_index << "," << bit2_index << "," << bit3_index << ")";
				codeword_seq = Input_codeword_seq;
				if ((bit1_index == bit2_index) && (bit2_index == bit3_index)) { // 3 same
					//cout << "0";
					codeword_seq.at(bit1_index) ^= 1;
					error_counter = Input_Node.D_z + 1;
					metric =
						Input_Node.metric
						+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index];
				}
				else if ((bit1_index != bit2_index) && (bit2_index == bit3_index)) {    // 2 same 1 different
					//cout << "0";
					codeword_seq.at(bit1_index) ^= 1;
					codeword_seq.at(bit2_index) ^= 1;
					error_counter = Input_Node.D_z + 2;
					metric =
						Input_Node.metric
						+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
						+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index];
				}
				else { 
					//cout << "1";
					// 3 different
					//if (bit1_index == bit2_index) cout << "a";
					//if (bit1_index == bit3_index) cout << "b";
					//if (bit3_index == bit2_index) cout << "c";

					codeword_seq.at(bit1_index) ^= 1;
					codeword_seq.at(bit2_index) ^= 1;
					codeword_seq.at(bit3_index) ^= 1;
					error_counter = Input_Node.D_z + 3;
					metric =
						Input_Node.metric
						+ Metric_Table._matrix[Hard_RX.at(bit1_index) ^ 1][bit1_index]
						+ Metric_Table._matrix[Hard_RX.at(bit2_index) ^ 1][bit2_index] + Metric_Table._matrix[Hard_RX.at(bit3_index) ^ 1][bit3_index];
				}
				//cout << endl;
				// Encode
				for (size_t j(message_length); j < Decoding_info.Control_Level; ++j) {
					if ((bit1_index == bit2_index) && (bit3_index == bit2_index)) {
						if (Sorted_G._matrix[bit1_index][j] == 1)
							codeword_seq.at(j) ^= 1;
					}
					else if ((bit1_index != bit2_index) && (bit2_index == bit3_index)) {
						if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
							codeword_seq.at(j) ^= 1;
					}
					else {
						int temp = Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j];
						temp ^= Sorted_G._matrix[bit3_index][j];
						if (temp == 1) codeword_seq.at(j) ^= 1;
					}

					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						++error_counter;
						if (error_counter > Decoding_info.Constraint_j) {
							++Decoding_info.DM_STE;
							break;
						}
						metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (metric > best_codeword.metric) break;
					}
				}

				Decoding_info.Binary_STE += (Decoding_info.Control_Level - message_length);

				// Check the number of errors within control level
				if ((error_counter <= Decoding_info.Constraint_j) && (metric <= best_codeword.metric)) {
					for (size_t j(Decoding_info.Control_Level); j < codeword_length; ++j) {
						if ((bit1_index == bit2_index) && (bit3_index == bit2_index)) {
							if (Sorted_G._matrix[bit1_index][j] == 1)
								codeword_seq.at(j) ^= 1;
						}
						else if ((bit1_index != bit2_index) && (bit2_index == bit3_index)) {
							if ((Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j]) == 1)
								codeword_seq.at(j) ^= 1;
						}
						else {
							int temp = Sorted_G._matrix[bit1_index][j] ^ Sorted_G._matrix[bit2_index][j];
							temp ^= Sorted_G._matrix[bit3_index][j];
							if (temp == 1)
								codeword_seq.at(j) ^= 1;
						}

						if (codeword_seq.at(j) != Hard_RX.at(j)) {
							metric += Metric_Table._matrix[codeword_seq.at(j)][j];
							if (metric > best_codeword.metric) 
								break;
						}
					}

					Decoding_info.STE += (codeword_length - Decoding_info.Control_Level); // for 3 bit flipping
					if (bit3_index == bit2_index) --Decoding_info.STE; // for 2 bit flipping
					if (bit1_index == bit2_index)--Decoding_info.STE; // for 1 bit flipping
					++Decoding_info.CandidateCodeWord;

					if (metric < best_codeword.metric) {
						best_codeword.metric = metric;
						best_codeword_seq = codeword_seq;
					}
				}

			}
			//cout << endl;
		}
	}
	//cout << " End" << endl << endl;
	best_codeword.message_bits.assign(
		best_codeword_seq.begin(),
		best_codeword_seq.begin() + message_length);
	return best_codeword;
}

void A_star_PC_out_CBC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0), // rx_signal hard decision 結果
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);
	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length),
		temp_Node(message_length); //CBC 回傳的best node
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	decoding_info.First_nonzero_metric = 0;
	decoding_info.skip_Flag = FALSE;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double PNC = 0;
	bool Update_Flag;
	size_t Update_Num(0), Non_Update_Num(0);
	for (size_t i(message_length- decoding_info.Constraint_i- decoding_info.CBC_FlippingBit-1); i < message_length; ++i) {
		PNC += abs(decoding_info.rx_signal_seq.at(i));
		//cout << "a";
	}
	//cout << endl;
	// for OSD-i
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
	}
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	decoding_info.Counter = 0;
	do {
		Pointer = Stack.front();
		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Update_Flag = FALSE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) && 
				(Child_Node.metric < Best_Goal.metric) && 
				(Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}

				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}

				
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}
				
				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;
				
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				
			}
			else if ((Child_Node.level < message_length) && 
				(Child_Node.metric < Best_Goal.metric) && 
				(Child_Node.D_z == decoding_info.Constraint_i)){
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				
				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G, 
							Metric_Table, 
							codeword_seq, 
							Hard_RX, 
							Child_Node, 
							Best_Goal, 
							decoding_info);
				}
				else if(decoding_info.CBC_FlippingBit == 2){
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G, 
							Metric_Table, 
							codeword_seq, 
							Hard_RX, 
							Child_Node, 
							Best_Goal, 
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
						    1);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;

				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				
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
		//if (Best_Goal.metric < PNC) break;
	} while (!Stack.empty());

	
	Non_Update_Num = 0;
	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
		decoding_info.skip_Flag = TRUE;
	}
	decoding_info.New_OSC_Alpha = Best_Goal.metric / decoding_info.First_nonzero_metric; //PoHan
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

}

void A_star_PC_out_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	double metric_total(0); //PoHan
	decoding_info.First_nonzero_metric = 0;

	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	/*
	double K_amount = 0;
	for (int i = message_length - 5; i < message_length; ++i) {
		K_amount += abs(decoding_info.Sorted_R.at(i));
	}
	K_amount /= (5);
	double K1 = abs(decoding_info.Sorted_R.at(message_length / 2));
	double K2 = abs(decoding_info.Sorted_R.at(message_length));
	double K3 = abs(decoding_info.Sorted_R.at(message_length*1.5));
	double parameterK = K2 * (K1 - K2)*(K2 - K3);

	double Mean = abs(decoding_info.Sorted_R.at(message_length - 1));
	double Ave_LRIP = 0;
	for (int i = message_length; i < codeword_length; ++i) {
		Ave_LRIP += abs(decoding_info.Sorted_R.at(i));
	}
	Ave_LRIP /= message_length;

	double MRIP = 0;
	double LRIP = 0;
	*/
	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}

	int CRC_Period = StartPeriod;
	int CRC_Period_Add = PeriodLength;
	vector<__int8> MM;
	double Metric_Temp = 0;

	// 開始 Tree Search
	do {
		if (CRC_Check == TRUE) {
			if (decoding_info.Counter > CRC_Period) {
				//cout << Best_Goal.metric / OSC_metric_thr << endl;
				//cout << "o";
				if (MM != Best_Goal.message_bits) {
					//vector <size_t> Temp = Location_Index;
					//cout << "o";

					MM = Best_Goal.message_bits;
					Systematic_Linear_Block_Code_Encoder(Sorted_G, MM, codeword_seq);
					Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

					//if (Location_Index != Temp) cout << "GG";

					vector<__int8> Est_M(decoding_info.estimated_codeword.begin(), decoding_info.estimated_codeword.begin() + MM.size());

					if (decoding_info.CRC_Examination(Est_M, decoding_info.CRC) == TRUE) {
						break;
					}
				}
				CRC_Period += CRC_Period_Add;
			}
		}

		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;
				/*
				// Reach Control level to check hamming distance
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {

					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (size_t j(message_length); j < decoding_info.Control_Level; ++j) {
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {

					}
				}
				*/
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {

					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//cout << "a2";
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}

				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//cout << "p";
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				//cout << Child_Node.metric << "/" << Pointer.metric << " ";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
				//cout << "(" << Stack.size() << ") " << endl;
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;
	} while (!Stack.empty());
	/*
	if (Best_Goal.D_z == 0) ++decoding_info.Dz_0_Number;
	else if (Best_Goal.D_z == 1) ++decoding_info.Dz_1_Number;
	else if (Best_Goal.D_z == 2) ++decoding_info.Dz_2_Number;
	else if (Best_Goal.D_z == 3) ++decoding_info.Dz_3_Number;
	else if (Best_Goal.D_z == 4) ++decoding_info.Dz_4_Number;
	else if (Best_Goal.D_z == 5) ++decoding_info.Dz_5_Number;
	else {

	}
	vector<__int8> Position;
	for (int i = 0; i < message_length; ++i) {
		if (Best_Goal.message_bits.at(i) != Hard_RX.at(i)) Position.push_back(i);
	}
	for (int i = 0; i < Position.size(); ++i) {
		if (i == 0 && Position.at(i) < decoding_info.D1)decoding_info.D1 = Position.at(i);
		else if (i == 1 && Position.at(i) < decoding_info.D2)decoding_info.D2 = Position.at(i);
		else if (i == 2 && Position.at(i) < decoding_info.D3)decoding_info.D3 = Position.at(i);
		else if (i == 3 && Position.at(i) < decoding_info.D4)decoding_info.D4 = Position.at(i);
		else if (i == 4 && Position.at(i) < decoding_info.D5)decoding_info.D5 = Position.at(i);
	}*/

	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//PoHan
	
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	decoding_info.New_OSC_Alpha = Best_Goal.metric / decoding_info.First_nonzero_metric; //PoHan


	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void A_star_PC_out_CBC_OSC_Verified(MATRIX<__int8> &G, DECODING_INFO &decoding_info)  // 修正部分狀況沒有套用DM-I的問題
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序
	//system("pause");

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	bool Update_Flag;
	size_t Update_Num(0), Non_Update_Num(0);
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	
	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//decoding_info.DoubleDecoder = TRUE;
	// 開始 Tree Search
	//vector<__int8> BestGoalTemp;
	//decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());

	do {
		
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Update_Flag = FALSE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;
				
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, decoding_info);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
					if (Best_Goal.metric >= OSC_metric_thr) {
						decoding_info.num_Best_Node_Update += decoding_info.temp_num_Deleted_Candidate_in_Stack;
					}
				}
				
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							OSC_metric_thr,
							1,
							decoding_info);
					/*
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
					*/
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							OSC_metric_thr,
							decoding_info);
					/*
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
					*/
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
				}
				//cout << "p";
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
			//cout << "C";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Place_Node(Stack, Child_Node, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;
		
	} while (!Stack.empty());
	//cout << "Counter: " << decoding_info.Counter << endl;

	
	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	
	//cout << decoding_info.code_seq.size() << "," << decoding_info.Sorted_R.size() << endl;
	//system("pause");

	
	//Error Test
	/*
	for (int i = 0; i < 192; ++i) {
		//cout << decoding_info.Sorted_R.at(i) << ",";
		if (decoding_info.Sorted_R.at(i) < 0 && decoding_info.code_seq.at(i) == 0) decoding_info.Error_Accumulation.at(i)++;
		else if (decoding_info.Sorted_R.at(i) > 0 && decoding_info.code_seq.at(i) == 1) decoding_info.Error_Accumulation.at(i)++;
	}*/
	//cout << endl;

	/*
	int Total_Error = 0;
	int MRIP_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size() / 2; ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++MRIP_Error;
			//cout << "index: " << i << ", ";
		}
	}
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++Total_Error;
			//cout << "index: " << i << ", ";
		}
	}
	//cout << MRIP_Error << endl;
	*/
	/*
	int Total_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != codeword_seq.at(i)) {
			++Total_Error;
			cout << "index: (" << i << ", " << decoding_info.Sorted_R.at(i) << ") / ";
		}
	}*/

	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;
	
	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void  A_star_PC_out_CBC_OSC_Verified_Function(MATRIX<__int8> &G, DECODING_INFO &decoding_info, vector<double> &Rx, vector<__int8> &outoputcodeword, long double &metric, double &Average_LLR)  // 修正部分狀況沒有套用DM-I的問題
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(Rx, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	//cout << "M:" << abs(decoding_info.Sorted_R.at(message_length)) << endl;
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序
	//system("pause");

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(Rx.at(i));
	}
	Average_LLR = OSC_metric_thr / codeword_length;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//decoding_info.DoubleDecoder = TRUE;
	// 開始 Tree Search
	//vector<__int8> BestGoalTemp;
	decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());

	do {

		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//cout << "p";
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;

	} while (!Stack.empty());

	metric = Best_Goal.metric;

	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);

	Desort_Function(Location_Index, codeword_seq, outoputcodeword);

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}


void A_star_PC_out_CBC_OSC_MultiDecoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while (Minus_i-- != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}


}
//A_star_PC_out_CBC_OSCf_Adaptive_i

void A_star_PC_out_CBC_OSC_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {

	Adaptive_I_Decode_Info Adaptive_Info;
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Adaptive_Info.Hard_RX = Hard_RX;
	Adaptive_Info.MRIP_codeword = MRIP_codeword;
	Adaptive_Info.Sorted_G = Sorted_G;
	Adaptive_Info.Metric_Table = Metric_Table;
	Adaptive_Info.Best_Goal = Best_Goal;
	Adaptive_Info.Best_Goal.metric = DBL_MAX;

	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Adaptive_Info.Sorted_G, Location_Index, Adaptive_Info.Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	Adaptive_Info.OSC_metric_thr = 0;

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Adaptive_Info.Metric_Table._matrix[0][i] != 0) Adaptive_Info.Hard_RX.at(i) = 1;

		// for OSC threshold
		Adaptive_Info.OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq.assign(Adaptive_Info.Hard_RX.begin(), Adaptive_Info.Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, message_seq, Adaptive_Info.MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	Adaptive_Info.OSC_metric_thr = decoding_info.OSC_Alpha*Adaptive_Info.OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		Adaptive_Info.OSC_metric_thr += Adaptive_Info.Metric_Table._matrix[0][i] + Adaptive_Info.Metric_Table._matrix[1][i];
	}
	decoding_info.DoubleDecoder = TRUE;

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;

	//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
	decoding_info.phase = 1;
	A_star_PC_out_CBC_OSC_Block(G, decoding_info, Adaptive_Info);  //   1st
	//cout << "(1):" << decoding_info.Counter << endl;
	//cout << "(1):" << decoding_info.Counter << endl;
	if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
		//cout << "(X2)" << endl;
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
		decoding_info.phase = 2;
		A_star_PC_out_CBC_OSC_Block(G, decoding_info, Adaptive_Info); //   2nd
		//cout << "(2):" << decoding_info.Counter << endl;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder3_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder3_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder2_i;
			//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
			decoding_info.phase = 3;
			A_star_PC_out_CBC_OSC_Block(G, decoding_info, Adaptive_Info);  //    3rd
			if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder4_i != 0)) {
				decoding_info.Constraint_i = Temp_i;
				decoding_info.CBC_FlippingBit = Temp_CBC;
				Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder4_i;
				while ((Minus_i--) != 0) {
					if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
					else --decoding_info.CBC_FlippingBit;
				}
				decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder3_i;
				//Adaptive_Info.OSC_metric_thr *= Adaptive_i_Parameter;
				A_star_PC_out_CBC_OSC_Block(G, decoding_info, Adaptive_Info);  //    4th
			}
		}
	}
	decoding_info.CBC_FlippingBit = Temp_CBC;
	decoding_info.Constraint_i = Temp_i;


	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Adaptive_Info.Sorted_G, Adaptive_Info.Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	decoding_info.estimated_codeword_66 = decoding_info.estimated_codeword;

	decoding_info.STE = decoding_info.STE / (double)message_length;
	/*
	decoding_info.STE_1 = decoding_info.STE_1 / (double)message_length;
	decoding_info.STE_2 = decoding_info.STE_2 / (double)message_length;
	decoding_info.STE_3 = decoding_info.STE_3 / (double)message_length;
	*/
	decoding_info.STE_1 = decoding_info.STE;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.COM_1 = decoding_info.COM;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

}

void A_star_PC_out_CBC_OSC_Block(MATRIX<__int8> &G, DECODING_INFO &decoding_info, Adaptive_I_Decode_Info &Adaptive_info)  // A_star_PC_out_CBC_OSC_DoubleDecoder的block
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);

	
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	if (decoding_info.Cancelled_Candidate_i != 0) {
		Operater_Deletion = TRUE;
		Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
		Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
	}
	//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
	// 開始 Tree Search
	do {
		//cout << Stack.size() << endl;
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
			//cout << "C:" << Stack.size() << endl;
			// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
			if ((Operater_Deletion == TRUE) && (Pointer.level == Level_k_previous) && (Pointer.D_z <= Difference)) {
				if (Pointer.level != (message_length - 1)) Stack.erase(Stack.begin());
				//cout << "!";
				break;
			}
			// End
			Extend_Node_Procedure(Pointer, Child_Node, Adaptive_info.Metric_Table, new_bit);
			if (new_bit != Adaptive_info.Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now
			//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "(A1)";
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if (Operater_Deletion == TRUE && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
				// End
				++decoding_info.Counter;
				codeword_seq = Adaptive_info.MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
						Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				if (decoding_info.phase == 1) {
					decoding_info.STE_1 += (codeword_length - message_length);
				}
				else if (decoding_info.phase == 2) {
					decoding_info.STE_2 += (codeword_length - message_length);
				}
				else if (decoding_info.phase == 3) {
					decoding_info.STE_3 += (codeword_length - message_length);
				}
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Adaptive_info.Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = Adaptive_info.MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Adaptive_info.Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				if (decoding_info.phase == 1) {
					decoding_info.STE_1 += (codeword_length - Child_Node.level);
				}
				else if (decoding_info.phase == 2) {
					decoding_info.STE_2 += (codeword_length - Child_Node.level);
				}
				else if (decoding_info.phase == 3) {
					decoding_info.STE_3 += (codeword_length - Child_Node.level);
				}
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Adaptive_info.Sorted_G,
							Adaptive_info.Metric_Table,
							codeword_seq,
							Adaptive_info.Hard_RX,
							Child_Node,
							Adaptive_info.Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				if (Operater_Deletion == FALSE || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Adaptive_info.Hard_RX.at(j)) {
							Child_Node.metric += Adaptive_info.Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Adaptive_info.Best_Goal.metric) break;
						}
					}
					//cout << Child_Node.metric << endl;
					if (temp_Node.metric < Child_Node.metric)Child_Node = temp_Node;
				}
				else Child_Node = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Child_Node, Adaptive_info.Best_Goal, Stack);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Adaptive_info.Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
			//cout << "(A3)";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}
		if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr) {
			//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
			break;
		}
	} while (!Stack.empty());
	//cout << decoding_info.Counter << endl;
	if (Adaptive_info.Best_Goal.metric < Adaptive_info.OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	//cout << endl << endl;
}

/*
void A_star_PC_out_CBC_OSC_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);


	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	
	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	decoding_info.Cancelled_Candidate_i = 0;
	int message_2 = message_length - 2;
	__int8 Constraint_i_choice = 0;
	Best_Goal.metric = FLT_MAX;

	for (__int8 DeTime = 0; DeTime < decoding_info.Adaptive_i_Iteration; DeTime++) {
		//DeTime = 1;
		if (DeTime == 0) {
			if (decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i == 2) {
				decoding_info.Constraint_i--;
				decoding_info.CBC_FlippingBit--;
				Constraint_i_choice = 2;
			}
			else if (decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i == 1)
			{
				decoding_info.Constraint_i--;
				Constraint_i_choice = 1;
			}
			else {
				Constraint_i_choice = 0;
			}
		}
		if (DeTime == 1) {
			if((decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i == 1 )&&(decoding_info.Adaptive_i_Iteration - DeTime >= 2)) decoding_info.Constraint_i--;
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		}
		if (DeTime == 2) decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder2_i;
		int Difference = decoding_info.Cancelled_Candidate_i - decoding_info.CBC_FlippingBit;
		//cout << DeTime << ": " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
		NODE_PATH
			Pointer(message_length),
			Child_Node(message_length),
			temp_Node(message_length);
		vector<NODE_PATH> Stack(1, Pointer);

		//Best_Goal.metric = FLT_MAX;

		int CRC_Period = StartPeriod;
		int CRC_Period_Add = PeriodLength;
		vector<__int8> MM;
		double Metric_Temp = 0;
		//cout << endl << endl << endl;
		// 開始 Tree Search
		//cout << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
		do {
		    //if(Pointer)
			Pointer = Stack.at(0);
			if (Pointer.level == (message_length - 1)) {
				Stack.erase(Stack.begin());
			}
			//cout << Stack.size() << " ";
			//++decoding_info.Counter;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
				if ((DeTime != 0) && (Pointer.level == message_2) && (Pointer.D_z <= Difference)) {
					Stack.erase(Stack.begin());
					break;
				}
				// End

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}

				++decoding_info.STE;
				// Child_Node: The node we are examining now

				// Reach level k
				if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if (decoding_info.Cancelled_Candidate_i != 0 && Child_Node.D_z <= decoding_info.Cancelled_Candidate_i) continue;
					// End
					++decoding_info.Counter;
					codeword_seq = MRIP_codeword;
					// DM-I: Reach Control level to check hamming distance
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					error_counter = Child_Node.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
					}
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX.at(j)) {
							Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
							if (Child_Node.metric > Best_Goal.metric) break;
						}
					}

					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
					//cout << "B";
					//if (decoding_info.Cancelled_Candidate_i == 0 || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
					++decoding_info.Counter;
					for (__int16 j(Child_Node.level); j < message_length; ++j) {
						Child_Node.message_bits.at(j) = Hard_RX.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

					codeword_seq = MRIP_codeword;
					for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node.level);
					//}
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Sorted_G,
								Metric_Table,
								codeword_seq,
								Hard_RX,
								Child_Node,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Sorted_G,
								Metric_Table,
								codeword_seq,
								Hard_RX,
								Child_Node,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Sorted_G,
								Metric_Table,
								codeword_seq,
								Hard_RX,
								Child_Node,
								Best_Goal,
								decoding_info
								, 1);
					}

					else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
					//if (decoding_info.Cancelled_Candidate_i == 0 || Child_Node.D_z > decoding_info.Cancelled_Candidate_i) {
					if (DeTime == 0) {
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX.at(j)) {
								Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
								if (Child_Node.metric > Best_Goal.metric) break;
							}
						}
						if (temp_Node.metric < Child_Node.metric)
							Child_Node = temp_Node;
					}
					else {
						Child_Node = temp_Node;
					}
					//}
					Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
					//cout << "p";
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
					//cout << "C";
					//cout << Child_Node.metric << "/" << Pointer.metric << " ";
					if (Child_Node.metric != Pointer.metric)
						Place_Node(Stack, Child_Node, decoding_info);
					else {
						Stack.at(0) = Child_Node;
						++decoding_info.COM;
					}
					//cout << "(" << Stack.size() << ") " << endl;
				}
			}
			//if((Cancelled_Candidate_i!=0) && (Stack.at(0).level== message_1) && (Stack.at(0).D_z <= ))
			if (Best_Goal.metric < OSC_metric_thr) break;
		} while (!Stack.empty());

		if (DeTime == 0) {
			if (Constraint_i_choice == 2) {
				decoding_info.Constraint_i++;
				decoding_info.CBC_FlippingBit++;
			}
			else if (Constraint_i_choice == 1)
			{
				decoding_info.Constraint_i++;
			}
			else {

			}
		}
		if (DeTime == 1 && decoding_info.Adaptive_i_Iteration - DeTime >= 2) {
			cout << "??" << endl;
			if (decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i == 0) decoding_info.Constraint_i++;
		}
		decoding_info.OSC_Ratio = Best_Goal.metric / OSC_metric_thr;
		if (decoding_info.OSC_Ratio < Adaptive_i_Parameter) break;

	}
	//cout << decoding_info.Adaptive_i_Iteration<<": "<< decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}
*/
void A_star_PC_out_CBC_OSC_Dynamic_I(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	/*
	if (Dynamic_I_Constraint == TRUE) {
		decoding_info.bit1 = -1;
		decoding_info.bit2 = -1;
		decoding_info.bit3 = -1;
		decoding_info.CBC1 = -1;
	}*/
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;
   
	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	decoding_info.Counter = 0;

	// 開始 Tree Search
	//cout << decoding_info.bit1 << "," << decoding_info.bit2 << "," << decoding_info.bit3 << "," << decoding_info.CBC1 << endl;
	do {
		decoding_info.Counter++;
		
		Pointer = Stack.at(0);
		//cout << Stack.size() << " ";
		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
	
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			/*
			if (decoding_info.Counter > message_length) {
				if (Child_Node.D_z == 1 && Child_Node.level <= decoding_info.bit1 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
				if (Child_Node.D_z == 2 && Child_Node.level <= decoding_info.bit2 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
				if (Child_Node.D_z == 3 && Child_Node.level <= decoding_info.bit3 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
			}*/
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

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
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info); // decoding_info.CBC1
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
				//cout << "(" << Stack.size() << ") " << endl;
			}
		}

		if (Best_Goal.metric < OSC_metric_thr) break;
	} while (!Stack.empty());

	decoding_info.TotalCounter += decoding_info.Counter;
	//cout << decoding_info.Counter << endl;
	


	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void A_star_PC_out_CBC_OSC_Test2(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;

	/*
	for (int i = 0; i < codeword_length; ++i) {
		if(abs(decoding_info.rx_signal_seq.at(i))>1) decoding_info.rx_signal_seq.at(i)*=1.3;
	}
	*/
	//cout << decoding_info.rx_signal_seq.size();
	/*
	if (Adapting_I_Constraint == TRUE) {
		decoding_info.bit1 = -1;
		decoding_info.bit2 = -1;
		decoding_info.bit3 = -1;
		double delta = 0.3;
		decoding_info.bit1Threshold = sqrt(decoding_info.var) * (5 - 0) - 1;
		decoding_info.bit2Threshold = sqrt(decoding_info.var) * (4 - 0) - 1;
		decoding_info.bit3Threshold = sqrt(decoding_info.var) * (3.5 - delta) - 1;
		decoding_info.CBC1 = -1;
		decoding_info.CBC1Threshold = sqrt(decoding_info.var)*(3 - delta) - 1;
	}
	*/
	/*
	if (Dynamic_I_Constraint == TRUE) {
		decoding_info.bit1 = -1;
		decoding_info.bit2 = -1;
		decoding_info.bit3 = -1;
		decoding_info.CBC1 = -1;
	}
	*/
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}


	//Yin
	//decoding_info.Unsorted_Est.assign(Hard_RX.begin(), Hard_RX.end());

	//double parameter = decoding_info.var*sqrt(decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;
	//OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr*log10(message_length)*parameter /12;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	decoding_info.Counter = 0;
	
	double temp = 0;
	int Max = 679120;   // 192:64593560 128:679120
	int Period_Start = Max / 20;
	int Period_increase = Max / 20;
	double Threshold_Add = OSC_metric_thr * 0.012;
	//int index = 0;
	//double beta = 1.25;
	//double beta_increase = 0.01;
	
	//cout << decoding_info.bit1 << "," << decoding_info.bit2 << "," << decoding_info.bit3 << "," << decoding_info.CBC1 << endl;

	// 開始 Tree Search
	do {

		//cout << decoding_info.var << endl;
		//if (Stack.size() > 1) temp = Stack.at(1).metric;
		//if ((Stack.size() > 1 )&& (Stack.at(Stack.size() - 1).metric > Best_Goal.metric) && (decoding_info.Counter > message_length)) cout << "+";
		
		if (decoding_info.Counter == message_length + 1) temp = Best_Goal.metric;
		if (decoding_info.Counter > Period_Start) {
			OSC_metric_thr += Threshold_Add;
			Period_Start += Period_increase;
		}
	
		


		// 自己想的
		/*
		if (decoding_info.Counter == message_length + 1)  temp = Best_Goal.metric;
		if (decoding_info.Counter > Period_Start) {
			if ((temp / OSC_metric_thr) < (beta)) {
				break;
			}
			//if ((temp / OSC_metric_thr) > 2.5 && Period_Start/ Max > 9) break;
			temp = Best_Goal.metric;
			beta += beta_increase;
			Period_Start += Period_increase;
		}
		*/
		/*
		if (decoding_info.Counter == message_length+1)  temp = Best_Goal.metric;
		if (decoding_info.Counter > Period_Start) {
			if (temp != Best_Goal.metric) {
				//cout << temp / Best_Goal.metric << endl;
				int j;
				if ((temp / OSC_metric_thr) - 1 >= 2) j = 19;
				else j = (int)(((temp / OSC_metric_thr) - 1) * 10);
				decoding_info.FlipRecorder[j][index]++;
				temp = Best_Goal.metric;
			}
			index++;
			Period_Start += Period_increase;
		}
		*/
		decoding_info.Counter++;

		Pointer = Stack.at(0);
		//cout << Stack.size() << " ";
		if (Pointer.level == (message_length - 1)) {
			Stack.erase(Stack.begin());
		}

		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			//cout << "(" << Stack.size() << "," << Pointer.level << ")";
			//Yin
			//cout << "(" << Pointer.D_z << ")";

			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			//cout << (int)new_bit << ":" << Stack.size() << " ";
			++decoding_info.STE;
			/*
			if (decoding_info.Counter > message_length) {
				if (Child_Node.D_z == 1 && Child_Node.level <= decoding_info.bit1 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
				if (Child_Node.D_z == 2 && Child_Node.level <= decoding_info.bit2 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
				if (Child_Node.D_z == 3 && Child_Node.level <= decoding_info.bit3 && Metric_Table._matrix[new_bit][Child_Node.level - 1] != 0) continue;
			}
			//
			*/
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

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
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info
						); // decoding_info.CBC1
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//cout << "p";
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << Child_Node.metric << "/" << Pointer.metric << " ";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
				//cout << "(" << Stack.size() << ") " << endl;
			}
		}

		if (Best_Goal.metric < OSC_metric_thr) break;
	} while (!Stack.empty());
	decoding_info.TotalCounter += decoding_info.Counter;

	//if(decoding_info.Counter > 679120) cout << decoding_info.Counter << endl <<endl;
	//cout << decoding_info.Counter << endl << endl;
	//cout << decoding_info.var << endl;
	//cout << temp << endl;
	/*
	if (temp == 0) { temp = Best_Goal.metric; }
	else
	{
		cout << temp / Best_Goal.metric << " / " << temp << "," << Best_Goal.metric << endl;
	}
	//cout << decoding_info.Counter << ": " << Best_Goal.metric << ", " << temp << ", " << OSC_metric_thr << "{" << temp / OSC_metric_thr << "}" << endl;
	//int node = (int)(Best_Goal.metric * 100 / temp);
	//if (node > 100) node = 100;
	//cout << decoding_info.Accumulate_Break.size() << "/" <<node << "(" << ;
	//system("pause");
	//decoding_info.Accumulate_Break.at(node)++;
	//cout << decoding_info.Counter << endl;
	//if (decoding_info.Counter > Max / 100) cout << "," << setw(50) << "F:" << " (" << decoding_info.Counter / (Max / 100) << "%)" << setprecision(3) << Best_Goal.metric / OSC_metric_thr << endl <<endl;
	//cout << endl;
	/*
	if (decoding_info.Counter == message_length) {
		decoding_info.Accumulate_Break.at(0)++;
		++decoding_info.Frame_Counter;
		decoding_info.OSC_Ratio;
	}
	else {

	}
	*/
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

// (Yin) Do not start CBC before level reached limieted length
void A_star_PC_out_Limited_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		limited_CBC_length((G.Row_number*(3 /4) ));  // (Yin) Do not start CBC before level reached limieted length
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
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;

	//cout << decoding_info.rx_signal_seq.size();
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//cout << "!";
	//decoding_info.Counter = 0;
	// 開始 Tree Search
	do {
		//++decoding_info.Counter;

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
			
			++decoding_info.STE;
			
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				// Pointer 到達 level k 直接做extension 
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
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				// Pointer 未到達 level k, 但已到達i的上限 -> 做CBC
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits_Limited(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							limited_CBC_length);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							0.5);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				// Pointer 既未到達level k,也未到達i的上限 -> 繼續做extension (general case)
	
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
		}

		if (Best_Goal.metric < OSC_metric_thr) break;
	} while (!Stack.empty());
	//cout << decoding_info.Counter<<endl;
	
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

////////////////////////////////////////////////////////////////////////////////////////// 正在做的/////////////////////////////////////////////////////////////////////////////////

// Yin

void A_star_PC_out_CBC_OSC_Test(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		MRIP_codeword(codeword_length, 0),
		Level_Record(decoding_info.Constraint_i + decoding_info.CBC_FlippingBit + 1, -1);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;

	//cout << decoding_info.rx_signal_seq.size();
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	//Yin
	//decoding_info.Unsorted_Est.assign(Hard_RX.begin(), Hard_RX.end());
	//Yin


	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}

	int Saturation_Maxsize = 64593560;  // 128:679120  196:64593560
	int Saturation_Sec_Maxsize = 3321960; // 128:41664 196:3321960
	int Saturation_Start = Saturation_Sec_Maxsize / (2 * decoding_info.var);
	int Saturation_Maxsize_Period = Saturation_Start / 3; //(Saturation_Maxsize * decoding_info.snr) / MRIP_codeword.size();
	//double alpha = OSC_metric_thr * 1.3;
	//double alpha = OSC_metric_thr * (1 + (2 / pow((double)decoding_info.snr, 2)));
	double alpha = OSC_metric_thr * (1.0 + (decoding_info.var*decoding_info.var));
	double beta = OSC_metric_thr * (1 + (2.5 / pow((double)decoding_info.snr, 2)));
	//cout << alpha << endl;
	//system("pause");
	//int counter = 0;
	//bool check = false;
	double Temp_Best = 0;
	//cout << Saturation_Maxsize_Period << "/" << alpha << "/" << OSC_metric_thr<< endl;

	//double Start, End;
	//Start = clock();
	decoding_info.Counter = 0;
	//double Lowest_Metric[3] = { 0 };
	//int Lowest_Level[3] = { 0 };
	//double TestMetric = 0;
	//int adder = 0;
	double TempSave = 0;
	int accumulate = 0;

	//cout << endl << endl;
	// 開始 Tree Search
	do {
		++decoding_info.Counter;
		cout << decoding_info.Counter << ": ";
		/*
		for (int i = 0; i < Stack.size(); ++i) {
			cout << setprecision(2) << Stack.at(i).metric << "(" << Stack.at(i).D_z << ") ";
		}
		*/
		cout << "Stacksize: "<<Stack.size();
		cout << endl <<endl;

		if ((decoding_info.Counter > Saturation_Start)) {
			Saturation_Start += Saturation_Maxsize_Period;
			//if(Stack.size()>3) cout << decoding_info.Counter << ": (" << Stack.at(1).metric << ", " << Stack.at(2).metric << ", " << Stack.at(3).metric << "), B: " << Best_Goal.metric << ", T: " << OSC_metric_thr << endl ;
			//if (Temp_Best == Best_Goal.metric && Best_Goal.metric < alpha) break;
			//else Temp_Best = Best_Goal.metric;
			//if (Stack.at(1).metric * 2.5 > OSC_metric_thr && decoding_info.Counter> Saturation_Maxsize/2) break;
			//if (Saturation_Start > Saturation_Maxsize / 8) Saturation_Maxsize_Period = Saturation_Start * 10;
			//if (decoding_info.Counter > Saturation_Maxsize / 4 && Best_Goal.metric < beta) break;
		}

		/*
		if (Stack.size() == 10) {
			for (int i = 0; i < 10; ++i) {
				cout << setprecision(4) << Stack.at(i).metric << " ";
			}
			cout << endl;
		}

		++counter;
		if (counter == 500) {
			Lowest_Metric[0] = Stack.at(1).metric;
			Lowest_Metric[1] = Stack.at(2).metric;
			Lowest_Metric[2] = Stack.at(3).metric;
			Lowest_Level[0] = Stack.at(1).level;
			Lowest_Level[1] = Stack.at(2).level;
			Lowest_Level[2] = Stack.at(3).level;
		}
		if ((counter % 1000 == 0) && (adder!=5)) {
			bool check = false;
			if ((Lowest_Metric[0] == Stack.at(1).metric) && (Lowest_Metric[1] == Stack.at(2).metric) && (Lowest_Metric[2] == Stack.at(3).metric)) ++adder;
			else {
				Lowest_Metric[0] = Stack.at(1).metric;
				Lowest_Metric[1] = Stack.at(2).metric;
				Lowest_Metric[2] = Stack.at(3).metric;
				Lowest_Level[0] = Stack.at(1).level;
				Lowest_Level[1] = Stack.at(2).level;
				Lowest_Level[2] = Stack.at(3).level;
				adder = 0;
				cout << " GGG" << counter << "GGG G:" << Best_Goal.metric << ", Dmin: " << Best_Goal.D_z << endl;
				cout << Lowest_Metric[0] << "," << Lowest_Level[0]<< " / " << Lowest_Metric[1] << "," << Lowest_Level[1] << " / " << Lowest_Metric[2] << ","<< Lowest_Level[2]  << endl << endl;
			}
			if (adder == 5) {
				TestMetric = Best_Goal.metric;
			}
		}

		//cout << Stack.size() << " ";
		*/
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
			//cout << (int)new_bit << ":" << Stack.size() << " ";
			++decoding_info.STE;
			//
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				cout << "a";

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
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				decoding_info.STE += (codeword_length - message_length);
				++decoding_info.CandidateCodeWord;

				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				// Yin
				if (Level_Record.at(1) > Stack.size()) Level_Record.at(1) = Stack.size();
				if (Level_Record.at(2) > Stack.size()) Level_Record.at(2) = Stack.size();
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//
				cout << "b";
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				//cout << "p";

				//Yin
				if (Level_Record.at(1) > Stack.size()) Level_Record.at(1) = Stack.size();
				if (Level_Record.at(2) > Stack.size()) Level_Record.at(2) = Stack.size();
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << Child_Node.metric << "/" << Pointer.metric << " ";
				cout << "c:(";
				cout << Stack.size() << "," << (int)Level_Record.at(1) << "," << (int)Level_Record.at(2) << "," << (int)Level_Record.at(3) <<")";
				cout << Child_Node.D_z << ")";
				if (Child_Node.metric != Pointer.metric) {
					if(decoding_info.Counter<= message_length)	Place_Node(Stack, Child_Node, decoding_info);
					else {
						cout << "&";
						Place_Node_ver2(Stack, Child_Node, decoding_info, Level_Record);
					}
				}
				else {
					cout << "d";
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
				//cout << "(" << Stack.size() << ") " << endl;
				cout << "t" << (int)Level_Record.at(1) << "t" << (int)(Level_Record.at(2)) << "t";
				if ((int)Level_Record.at(1) > (int)Stack.size()) Level_Record.at(1) = Stack.size();
				if ((int)Level_Record.at(2) > (int)Stack.size()) Level_Record.at(2) = Stack.size();
			}
		}

		if (Best_Goal.metric < OSC_metric_thr) break;
	} while (!Stack.empty());
	//End = clock();
	//cout << decoding_info.Counter << endl;
	//if (counter > 1000) cout << "Bad!: " << counter << ", " << Best_Goal.metric / OSC_metric_thr << ", B:" << Best_Goal.metric << " / OSC:" << OSC_metric_thr << endl;
	/*
	if (counter > 600) cout << endl;
	if ((End - Start) > 50) {
		cout << counter << endl;
		cout << "Total: " << counter << ", Metric: " << Best_Goal.metric << ", OSC: " << OSC_metric_thr;
	}
	*/
	//if (decoding_info.Counter > 1000000) cout << "Bad!: " << decoding_info.Counter << ", **" << Best_Goal.metric / OSC_metric_thr << "**, B:" << Best_Goal.metric<< " / OSC:" << OSC_metric_thr << endl;
	//else cout << "Good: " << decoding_info.Counter << ", " << Best_Goal.metric / OSC_metric_thr << endl;
	//cout << "Total: " << counter << ", Metric: " << Best_Goal.metric << ", OSC: " << OSC_metric_thr << endl;
	/*
	End = clock();
	if (End - Start > 500) cout << "***";

	cout << endl;
	//if (Best_Goal.metric < OSC_metric_thr) cout << "*";
	*/
	//cout << counter << " ";
	if (Best_Goal.D_z == 5) cout << " ! ";
	//if (counter == 147536) cout << Best_Goal.D_z << "/" << Best_Goal.D_zp << "/" << Best_Goal.metric << "/" << OSC_metric_thr << endl;
	//if (Best_Goal.metric < OSC_metric_thr) cout << Best_Goal.metric << " / " << OSC_metric_thr << endl;
	//else cout << "*" << Best_Goal.metric << " / " << OSC_metric_thr << endl;
	//cout << Best_Goal.metric << " / "<< OSC_metric_thr << endl;
	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	//Yin
	/*
	int ErrorAmount = 0;
	for (int i = 0; i < Best_Goal.message_bits.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++ErrorAmount;
	}
	if (ErrorAmount != 0) cout << endl << "Error:" << ErrorAmount << " ! " << endl << endl;
	cout << endl;
	*/

	//Yin

	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

//////////////////////////////////////////////////////////////////////////////////////////// END//////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////// 博瀚正在做的/////////////////////////////////////////////////////////////////////////////////

void A_star_1_stack_CBC_OSC(
	MATRIX<__int8>			&Sorted_G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Node, 
	NODE_PATH				&Pre_Best_Goal, 
	size_t					pc_i,
	double					metric_thr,
	DECODING_INFO			&decoding_info)
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
		Child_Node(message_length),
		temp_Node(message_length);;
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
				++decoding_info.CandidateCodeWord;
				//
				//Systematic_Linear_Block_Code_Encoder(G, Child_Node.message_bits, codeword_seq);
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
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z < pc_i)) {
				if (Pointer.metric != Child_Node.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Stack.at(0) = Child_Node;
					++decoding_info.COM;
				}
			}
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Pre_Best_Goal.metric) && (Child_Node.D_z == pc_i)) {
				//
				for (size_t j(Child_Node.level); j < message_length; ++j) 
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				
				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				
				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;
				//
				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				/*
				else {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);*/
			}
			//if (Best_Goal.metric < metric_thr) break;
		}
	} while (!Stack.empty());
	Pre_Best_Goal = Best_Goal;
}

void A_star_2_stack_CBC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	decoding_info.First_nonzero_metric = 0;
	double metric_total(0); //PoHan
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// for OSD
	for (size_t i(0); i < codeword_length; ++i) {
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
		metric_total += abs(decoding_info.rx_signal_seq.at(i));
	}

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

			for (size_t j(message_length); j < codeword_length; ++j) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				if (Pointer.metric > Best_Goal.metric) break;
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Child_Node.metric;
			}
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {
			decoding_info.STE += 2;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}

				if (Child_Node.metric < Best_Goal.metric) {
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
			A_star_1_stack_CBC_OSC(
				Sorted_G,
				Metric_Table,
				Hard_RX,
				MRIP_codeword,
				Pointer,
				Best_Goal_2,
				(size_t)(1 + 2),
				DBL_MAX,
				decoding_info);

			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Pointer.metric;
			}
		}
	} while (!Stack.empty());

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	decoding_info.New_OSC_Alpha = Best_Goal.metric / decoding_info.First_nonzero_metric; //PoHan
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;
	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void A_star_2_stack_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	decoding_info.First_nonzero_metric = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double metric_thr(0);
	double metric_total(0); //PoHan
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSD
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
		// for OSC
		metric_thr += abs(decoding_info.rx_signal_seq.at(i)); // for OSC
	}
	metric_total = metric_thr;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

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

			for (size_t j(message_length); j < codeword_length; ++j) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				if (Pointer.metric > Best_Goal.metric) break;
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Child_Node.metric;
			}
			if (Best_Goal.metric < metric_thr)	break;
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

				if (Child_Node.metric < Best_Goal.metric) {
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
			A_star_1_stack_CBC_OSC(
				Sorted_G,
				Metric_Table,
				Hard_RX,
				MRIP_codeword,
				Pointer,
				Best_Goal_2,
				(size_t)(1 + 2),
				metric_thr,
				decoding_info);

			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Pointer.metric;
			}
			if (Best_Goal.metric < metric_thr)	break;
		}
	} while (!Stack.empty());

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	
	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;
	

	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	decoding_info.New_OSC_Alpha = Best_Goal.metric / decoding_info.First_nonzero_metric; //PoHan
	if (isnan(decoding_info.New_OSC_Alpha) == 1) {
		cout << "\n Best_Goal metric: " << Best_Goal.metric;
		cout << "\n First nonzero metric: " << decoding_info.First_nonzero_metric;
		system("PAUSE");
	}
	
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void A_star_2_stack_CBC_OSC(
	MATRIX<__int8>			&G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Node,
	NODE_PATH				&Pre_Best_Goal,
	size_t					pc_i,
	double					metric_thr,
	DECODING_INFO			&decoding_info
	) {
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
	
	do{
		Pointer = Stack.at(0);
		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
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
			++decoding_info.CandidateCodeWord;
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			//if (Best_Goal.metric < metric_thr) break;
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
			Best_Goal_2 = Best_Goal;
			A_star_1_stack_CBC_OSC(
				G, 
				Metric_Table, 
				Hard_RX, 
				MRIP_codeword, 
				Pointer, 
				Best_Goal_2, 
				pc_i+2,
				metric_thr,
				decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			//if (Best_Goal.metric < metric_thr) break;
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}


void A_star_3_stack_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
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
		Best_Goal_2(message_length),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);
	Best_Goal.metric = FLT_MAX;
	decoding_info.First_nonzero_metric = 0;
	decoding_info.skip_Flag = FALSE;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	double metric_thr(0);
	double metric_total(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSD
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
		// for OSC
		metric_thr += abs(decoding_info.rx_signal_seq.at(i)); // for OSC
	}
	metric_total = metric_thr;

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	metric_thr = decoding_info.OSC_Alpha*metric_thr;

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
			//
			Stack.erase(Stack.begin());
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
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Pointer.metric;
			}
			if (Best_Goal.metric < metric_thr)	break;
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

				if (Child_Node.metric < Best_Goal.metric) {
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
			A_star_2_stack_CBC_OSC(
				Sorted_G,
				Metric_Table,
				Hard_RX,
				MRIP_codeword,
				Pointer,
				Best_Goal_2,
				1+1,
				metric_thr,
				decoding_info);

			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			if (Best_Goal.metric == Pointer.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Pointer.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Pointer.metric;
			}
			if (Best_Goal.metric < metric_thr)	break;
		}
	} while (!Stack.empty());

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);
	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;

	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
		decoding_info.skip_Flag = TRUE;
	}
	decoding_info.New_OSC_Alpha = Best_Goal.metric / decoding_info.First_nonzero_metric; //PoHan
	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

}


void A_star_PC_out_CBC_OSC_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info)  // 修正部分狀況沒有套用DM-I的問題
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);
	MATRIX<double> Fano_Metric_Table(2, message_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, Fano_Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序
	//system("pause");

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	bool Update_Flag;
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	
	//decoding_info.DoubleDecoder = TRUE;
	// 開始 Tree Search
	//vector<__int8> BestGoalTemp;
	//decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());
	do {

		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure_Fano(Pointer, Child_Node, Metric_Table, Fano_Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
					
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				/*
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack); 

				//cout << "p";
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				Place_Node_Fano(Stack, Child_Node, decoding_info);
				/*
				if (Child_Node.metric != Pointer.metric)
					Place_Node_Fano(Stack, Child_Node, decoding_info);
				else {
					Place_Node(Stack, Child_Node, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
				*/
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;

	} while (!Stack.empty());
	
	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);

	//cout << decoding_info.code_seq.size() << "," << decoding_info.Sorted_R.size() << endl;
	//system("pause");


	//Error Test
	/*
	for (int i = 0; i < 192; ++i) {
		//cout << decoding_info.Sorted_R.at(i) << ",";
		if (decoding_info.Sorted_R.at(i) < 0 && decoding_info.code_seq.at(i) == 0) decoding_info.Error_Accumulation.at(i)++;
		else if (decoding_info.Sorted_R.at(i) > 0 && decoding_info.code_seq.at(i) == 1) decoding_info.Error_Accumulation.at(i)++;
	}*/
	//cout << endl;

	/*
	int Total_Error = 0;
	int MRIP_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size() / 2; ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++MRIP_Error;
			//cout << "index: " << i << ", ";
		}
	}
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++Total_Error;
			//cout << "index: " << i << ", ";
		}
	}
	//cout << MRIP_Error << endl;
	*/
	/*
	int Total_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != codeword_seq.at(i)) {
			++Total_Error;
			cout << "index: (" << i << ", " << decoding_info.Sorted_R.at(i) << ") / ";
		}
	}*/

	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
} 

void A_star_multiple_basis_CBC_OSC(
	MATRIX<__int8>			&Sorted_G,
	MATRIX<double>			&Metric_Table,
	vector<__int8>			&Hard_RX,
	vector<__int8>			&MRIP_codeword,
	NODE_PATH				&Best_Goal,
	vector<NODE_PATH>		&Stack_CBC,
	double					OSC_metric_thr,
	int						Basis,
	DECODING_INFO			&decoding_info
	) {
	size_t
		message_length(Sorted_G.Row_number),
		codeword_length(Sorted_G.Col_number),
		error_counter(0),
		position_number(0);

	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);

	vector<NODE_PATH> Stack(1, Pointer);
	
	do {
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (Basis == 1) {
				Child_Node.base = Sorted_Base_1;
			}
			if (Basis == 2) {
				Child_Node.base = Sorted_Base_2;
			}
			if (Basis == 3) {
				Child_Node.base = Sorted_Base_3;
			}
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, Stack_CBC);
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				
				if (Basis == 1) {//為了最後一次的decode作準備
					Stack_CBC.insert(Stack_CBC.begin() + position_number, Child_Node);
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				if (Basis == 1) {
					temp_Node.base = Sorted_Base_1;
				}
				if (Basis == 2) {
					temp_Node.base = Sorted_Base_2;
				}
				if (Basis == 3) {
					temp_Node.base = Sorted_Base_3;
				}
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, Stack_CBC);

				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Place_Node(Stack, Child_Node, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) {
			decoding_info.DoubleDecoder = FALSE;
			break;
		}
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
			//}

	} while (!Stack.empty());
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
}

void A_star_2_Base_PC_out_CBC_OSC(MATRIX<__int8> &G, DECODING_INFO &decoding_info) 
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8> Sorted_G_Base1(G);
	MATRIX<__int8> Sorted_G_Base2(G);
	MATRIX<double> Metric_Table_Base1(2, codeword_length);
	MATRIX<double> Metric_Table_Base2(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC(1, Pointer_Base1);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, 2, 2, decoding_info);
	
	/*
	for (size_t i(0); i < Sorted_G_Base1.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base1.Col_number; j++) {
			cout << (int)Sorted_G_Base1._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	for (size_t i(0); i < Sorted_G_Base2.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base2.Col_number; j++) {
			cout << (int)Sorted_G_Base2._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	*/
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	size_t Update_Num(0), Non_Update_Num(0);
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword

	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	decoding_info.DoubleDecoder = TRUE;
	A_star_multiple_basis_CBC_OSC(
		Sorted_G_Base1,
		Metric_Table_Base1,
		Hard_RX_Base1,
		MRIP_codeword_Base1,
		Best_Goal,
		Stack_CBC,
		OSC_metric_thr,
		1,
		decoding_info
	);
	if (decoding_info.DoubleDecoder == TRUE) {
		A_star_multiple_basis_CBC_OSC(
			Sorted_G_Base2,
			Metric_Table_Base2,
			Hard_RX_Base2,
			MRIP_codeword_Base2,
			Best_Goal,
			Stack_CBC,
			OSC_metric_thr,
			2,
			decoding_info
		);
		/*
		if (Best_Goal.metric < (OSC_metric_thr * 1.5)) {
			A_star_multiple_basis_CBC_OSC(
				Sorted_G_Base2,
				Metric_Table_Base2,
				Hard_RX_Base2,
				MRIP_codeword_Base2,
				Best_Goal,
				Stack_CBC,
				OSC_metric_thr,
				2,
				decoding_info
			);
		}
		*/
		if (decoding_info.DoubleDecoder == TRUE) {
			decoding_info.CBC_FlippingBit = 2;
			// Compute permutation by using "quick sorting algorithm" 
			
			do {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Child_Node_CBC = Stack_CBC.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_CBC.erase(Stack_CBC.begin());
				//++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_CBC.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_CBC.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_CBC.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_CBC.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_CBC.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_CBC.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_CBC.metric)
					Child_Node_CBC = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_CBC, Best_Goal, Stack_CBC);
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					decoding_info.CBC_FlippingBit = 1;
					break;
				}
				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/

			} while (!Stack_CBC.empty());
			decoding_info.CBC_FlippingBit = 1;
		}
		
	}




	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	/*
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	*/
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

	
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	
}

void A_star_2_Base_PC_out_CBC_OSC_Parallel(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8> Sorted_G_Base1(G);
	MATRIX<__int8> Sorted_G_Base2(G);
	MATRIX<double> Metric_Table_Base1(2, codeword_length);
	MATRIX<double> Metric_Table_Base2(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC(1, Pointer_Base1);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, 2, 2, decoding_info);

	/*
	for (size_t i(0); i < Sorted_G_Base1.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base1.Col_number; j++) {
			cout << (int)Sorted_G_Base1._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	for (size_t i(0); i < Sorted_G_Base2.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base2.Col_number; j++) {
			cout << (int)Sorted_G_Base2._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	*/
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	size_t Next_Flag = Sorted_Base_1;
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword

	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	decoding_info.DoubleDecoder = TRUE;
	
	do {
		if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base1 = Stack_Base1.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base1.erase(Stack_Base1.begin());
			/*
			if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
				Stack_Base2.erase(Stack_Base2.begin());
			}
			*/
			//++decoding_info.Counter;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
				Child_Node_Base1.base = Sorted_Base_1;
				if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
					++Child_Node_Base1.D_z;
					Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
				}
				++decoding_info.STE;
				++decoding_info.Binary_STE;
				// Child_Node: The node we are examining now

				// Reach level k
				if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
					//cout << "A";
					++decoding_info.Counter;
					codeword_seq = MRIP_codeword_Base1;
					Next_Flag = Sorted_Base_2;
					// DM-I: Reach Control level to check hamming distance

					for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
						}
					}
					error_counter = Child_Node_Base1.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
					}
					decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
						//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
							Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base1.metric > Best_Goal.metric) break;
						}
					}
					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC);
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
					//cout << "B";
					++decoding_info.Counter;
					if (Child_Node_Base1.level < message_length - Multiple_Basis_Bits) {
						Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					}
					for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
						Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

					if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
						if (position_number == 0) {
							Stack_CBC.at(0) = Child_Node_Base1;
						}
						else {
							Stack_CBC.insert(Stack_CBC.begin() + position_number, Child_Node_Base1);
						}
						position_number++;
						//Place_Node(Stack_CBC, Child_Node, decoding_info);
					}

					codeword_seq = MRIP_codeword_Base1;
					for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node_Base1.level);
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Sorted_G_Base1,
								Metric_Table_Base1,
								codeword_seq,
								Hard_RX_Base1,
								Child_Node_Base1,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Sorted_G_Base1,
								Metric_Table_Base1,
								codeword_seq,
								Hard_RX_Base1,
								Child_Node_Base1,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Sorted_G_Base1,
								Metric_Table_Base1,
								codeword_seq,
								Hard_RX_Base1,
								Child_Node_Base1,
								Best_Goal,
								decoding_info,
								1);
					}
					else temp_Node.metric = DBL_MAX;
					temp_Node.base = Sorted_Base_1;
					Next_Flag = Sorted_Base_2;
					//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
							Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base1.metric > Best_Goal.metric) break;
						}
					}

					if (temp_Node.metric < Child_Node_Base1.metric)
						Child_Node_Base1 = temp_Node;
					Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC);

					//cout << "p";
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
					//cout << "C";
					if (Child_Node_Base1.metric != Pointer_Base1.metric)
						Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
					else {
						Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
						//Stack.at(0) = Child_Node;
						//++decoding_info.COM;
					}
					
					if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						else {
							//Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
							Stack_Base2.at(0) = Child_Node_Base1;
							++decoding_info.COM;
						}
					}
					
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		}
		else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
			Next_Flag = Sorted_Base_2;
		}
		else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			//++decoding_info.Counter;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
				Child_Node_Base2.base = Sorted_Base_2;
				if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
					++Child_Node_Base2.D_z;
					Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
				}
				++decoding_info.STE;
				++decoding_info.Binary_STE;
				// Child_Node: The node we are examining now

				// Reach level k
				if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
					//cout << "A";
					++decoding_info.Counter;
					codeword_seq = MRIP_codeword_Base2;
					Next_Flag = Sorted_Base_1;
					// DM-I: Reach Control level to check hamming distance

					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
						}
					}
					error_counter = Child_Node_Base2.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
					}
					decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
							Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base2.metric > Best_Goal.metric) break;
						}
					}
					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC);
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
					//cout << "B";
					++decoding_info.Counter;
					for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
						Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
					codeword_seq = MRIP_codeword_Base2;
					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node_Base2.level);
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info,
								1);
					}
					else temp_Node.metric = DBL_MAX;
					temp_Node.base = Sorted_Base_2;
					Next_Flag = Sorted_Base_1;
					//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
							Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base2.metric > Best_Goal.metric) break;
						}
					}

					if (temp_Node.metric < Child_Node_Base2.metric)
						Child_Node_Base2 = temp_Node;
					Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC);

					//cout << "p";
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
					//cout << "C";
					if (Child_Node_Base2.metric != Pointer_Base2.metric)
						Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
					else {
						Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
						//Stack.at(0) = Child_Node;
						//++decoding_info.COM;
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		}
		else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
			Next_Flag = Sorted_Base_1;
		}
	} while (!Stack_Base1.empty() || !Stack_Base2.empty());
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		decoding_info.CBC_FlippingBit = 2;
		// Compute permutation by using "quick sorting algorithm" 
		do {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Child_Node_CBC = Stack_CBC.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_CBC.erase(Stack_CBC.begin());
			//++decoding_info.Counter;
			codeword_seq = MRIP_codeword_Base1;
			for (size_t index(0); index < Child_Node_CBC.Diff_Index.size(); ++index) {
				codeword_seq.at(Child_Node_CBC.Diff_Index.at(index)) ^= 1;
				for (__int16 j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_CBC.Diff_Index.at(index)][j];
				}
			}
			//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

			decoding_info.STE += (codeword_length - Child_Node_CBC.level);
			//++decoding_info.CandidateCodeWord;

			// CBC
			if (decoding_info.CBC_FlippingBit == 1) {
				temp_Node
					= Control_Band_Check_1bit(
						Sorted_G_Base1,
						Metric_Table_Base1,
						codeword_seq,
						Hard_RX_Base1,
						Child_Node_CBC,
						Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 2) {
				temp_Node
					= Control_Band_Check_2bits(
						Sorted_G_Base1,
						Metric_Table_Base1,
						codeword_seq,
						Hard_RX_Base1,
						Child_Node_CBC,
						Best_Goal,
						decoding_info);
			}
			else if (decoding_info.CBC_FlippingBit == 3) {
				temp_Node
					= Control_Band_Check_3bits(
						Sorted_G_Base1,
						Metric_Table_Base1,
						codeword_seq,
						Hard_RX_Base1,
						Child_Node_CBC,
						Best_Goal,
						decoding_info,
						1);
			}
			else temp_Node.metric = DBL_MAX;
			temp_Node.base = Sorted_Base_1;
			//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
			for (size_t j(message_length); j < codeword_length; ++j) {
				if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
					Child_Node_CBC.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
					if (Child_Node_CBC.metric > Best_Goal.metric) break;
				}
			}

			if (temp_Node.metric < Child_Node_CBC.metric)
				Child_Node_CBC = temp_Node;
			Update_Best_Goal_Procedure(Child_Node_CBC, Best_Goal, Stack_CBC);
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				decoding_info.CBC_FlippingBit = 1;
				break;
			}
			//cout << "p";
			/*
			if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
				decoding_info.First_nonzero_metric = Best_Goal.metric;
			}
			else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
				decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
			}
			*/

		} while (!Stack_CBC.empty());
		decoding_info.CBC_FlippingBit = 1;
	}

	




	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	/*
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	*/
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;


	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);


}

void A_star_2_Base_PC_out_CBC_OSC_Latest(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8> Sorted_G_Base1(G);
	MATRIX<__int8> Sorted_G_Base2(G);
	MATRIX<double> Metric_Table_Base1(2, codeword_length);
	MATRIX<double> Metric_Table_Base2(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC(1, Pointer_Base1);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, 2, 2, decoding_info);

	/*
	for (size_t i(0); i < Sorted_G_Base1.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base1.Col_number; j++) {
			cout << (int)Sorted_G_Base1._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	for (size_t i(0); i < Sorted_G_Base2.Row_number; i++) {
		for (size_t j(0); j < Sorted_G_Base2.Col_number; j++) {
			cout << (int)Sorted_G_Base2._matrix[i][j];
		}
		cout << "\n";
	}
	cout << "\n";
	*/
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	size_t Update_Num(0), Non_Update_Num(0);
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword

	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	decoding_info.DoubleDecoder = TRUE;
	do {
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer_Base1 = Stack_Base1.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack_Base1.erase(Stack_Base1.begin());
		if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
			Stack_Base2.erase(Stack_Base2.begin());
		}

		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
			Child_Node_Base1.base = Sorted_Base_1;
			if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
				++Child_Node_Base1.D_z;
				Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node_Base1.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC);
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
					Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}

				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC.at(0) = Child_Node_Base1;
					}
					else {
						Stack_CBC.insert(Stack_CBC.begin() + position_number, Child_Node_Base1);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}

				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_Base1.metric)
					Child_Node_Base1 = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC);

				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node_Base1.metric != Pointer_Base1.metric)
					Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
				else {
					Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}

				if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
					if (Child_Node_Base1.metric != Pointer_Base1.metric)
						Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					else {
						Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						//Stack_Base2.at(0) = Child_Node_Base2;
						//++decoding_info.COM;
					}
				}

			}
		}
		if (Best_Goal.metric < OSC_metric_thr) {
			decoding_info.DoubleDecoder = FALSE;
			break;
		}
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
			//}
	} while (!Stack_Base1.empty());
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		do {
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			//++decoding_info.Counter;
			for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
				Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
				Child_Node_Base2.base = Sorted_Base_2;
				if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
					++Child_Node_Base2.D_z;
					Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
				}
				++decoding_info.STE;
				++decoding_info.Binary_STE;
				// Child_Node: The node we are examining now

				// Reach level k
				if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
					//cout << "A";
					++decoding_info.Counter;
					codeword_seq = MRIP_codeword_Base2;
					// DM-I: Reach Control level to check hamming distance

					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
						}
					}
					error_counter = Child_Node_Base2.D_z;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
					}
					decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
					if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
						++decoding_info.DM_STE;
						//cout << decoding_info.DM_STE <<" ";
						continue;
					}
					// if DM-I condition did not fit, then continue
					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
						for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
							//cout << "o";
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
						}
					}
					for (__int16 j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
							Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base2.metric > Best_Goal.metric) break;
						}
					}
					decoding_info.STE += (codeword_length - message_length);
					//++decoding_info.CandidateCodeWord;
					//cout << "a3";
					Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC);
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Did not reach level k, but reach i errors (compared with hard decision result)
				else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
					//cout << "B";
					++decoding_info.Counter;
					for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
						Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
					}
					//
					//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
					codeword_seq = MRIP_codeword_Base2;
					for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
						codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
						for (__int16 j(message_length); j < codeword_length; ++j) {
							codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
						}
					}
					//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

					decoding_info.STE += (codeword_length - Child_Node_Base2.level);
					//++decoding_info.CandidateCodeWord;

					// CBC
					if (decoding_info.CBC_FlippingBit == 1) {
						temp_Node
							= Control_Band_Check_1bit(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 2) {
						temp_Node
							= Control_Band_Check_2bits(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info);
					}
					else if (decoding_info.CBC_FlippingBit == 3) {
						temp_Node
							= Control_Band_Check_3bits(
								Sorted_G_Base2,
								Metric_Table_Base2,
								codeword_seq,
								Hard_RX_Base2,
								Child_Node_Base2,
								Best_Goal,
								decoding_info,
								1);
					}
					else temp_Node.metric = DBL_MAX;
					temp_Node.base = Sorted_Base_2;
					//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
					for (size_t j(message_length); j < codeword_length; ++j) {
						if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
							Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
							if (Child_Node_Base2.metric > Best_Goal.metric) break;
						}
					}

					if (temp_Node.metric < Child_Node_Base2.metric)
						Child_Node_Base2 = temp_Node;
					Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC);

					//cout << "p";
					/*
					if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
						decoding_info.First_nonzero_metric = Best_Goal.metric;
					}
					else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
						decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
					}
					*/
				}
				// Neither reach level k nor reach i error
				else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
					//cout << "C";
					if (Child_Node_Base2.metric != Pointer_Base2.metric)
						Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
					else {
						Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
						//Stack.at(0) = Child_Node;
						//++decoding_info.COM;
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		} while (!Stack_Base2.empty());
		if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
		if (decoding_info.DoubleDecoder == TRUE) {
			decoding_info.CBC_FlippingBit = 2;
			// Compute permutation by using "quick sorting algorithm" 

			do {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Child_Node_CBC = Stack_CBC.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_CBC.erase(Stack_CBC.begin());
				//++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_CBC.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_CBC.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_CBC.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_CBC.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_CBC,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_CBC.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_CBC.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_CBC.metric)
					Child_Node_CBC = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_CBC, Best_Goal, Stack_CBC);
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					decoding_info.CBC_FlippingBit = 1;
					break;
				}
				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/

			} while (!Stack_CBC.empty());
			decoding_info.CBC_FlippingBit = 1;
		}
	}
	



	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}
	

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	/*
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}
	*/
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

}


void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i(MATRIX<__int8> &G, DECODING_INFO &decoding_info) { //除了第一層的multibase 其他層都是各自長各自的
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length),
		Initial_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);
	

	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer_Base1 = Stack_Base1.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack_Base1.erase(Stack_Base1.begin());
		/*
		if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
			Stack_Base2.erase(Stack_Base2.begin());
		}
		*/
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
			Child_Node_Base1.base = Sorted_Base_1;
			if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
				++Child_Node_Base1.D_z;
				Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node_Base1.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				/*
				if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					Child_Node_Base1.base = Sorted_Base_1;
				}
				*/
				if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base1.at(0) = Child_Node_Base1;
					}
					else {
						//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
						Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
					Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_Base1.metric)
					Child_Node_Base1 = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node_Base1.metric != Pointer_Base1.metric)
					Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
				else {
					Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
				if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
					if (Stack_Base2.at(0).level == 0) {
						Stack_Base2.at(0) = Child_Node_Base1;
						++decoding_info.COM;
					}
					else {
						Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					}
					Child_Node_Base1.base = Sorted_Base_1;
				}
				
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) {
			decoding_info.DoubleDecoder = FALSE;
			break;
		}
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
			//}
	} while (!Stack_Base1.empty()); // phase 1
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		position_number = 0;
		do {
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			//++decoding_info.Counter;
			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Pointer_Base2);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;

						//
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
							}
							position_number++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}

						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		} while (!Stack_Base2.empty());
		if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
			Stack_Base1 = Stack_CBC_Base1;
			Stack_CBC_Base1.clear();
			/***開始decode****/
			position_number = 0;
			Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
			if (decoding_info.Cancelled_Candidate_i != 0) {
				Operater_Deletion = TRUE;
				Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
				Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
			}
			// 開始 Tree Search
			do {
				//cout << Stack.size() << endl;
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				if (Pointer_Base1.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack_Base1.erase(Stack_Base1.begin());
				}
				
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
					//cout << "C:" << Stack.size() << endl;
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now
					//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "(A1)";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;

						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "(A2)";
						//cout << "B";
						++decoding_info.Counter;
						
						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
							}
							position_number++;
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
						//cout << Adaptive_info.Best_Goal.metric << endl;
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "(A3)";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Stack_Base1.at(0) = Child_Node_Base1;
							++decoding_info.COM;
						}
						
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
					break;
				}
			} while (!Stack_Base1.empty()); // phase 2
			if (decoding_info.DoubleDecoder == TRUE) {
				Stack_Base2 = Stack_CBC_Base2;
				Stack_CBC_Base2.clear();
				/***開始decode****/
				position_number = 0;
				Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
				do {
					//cout << Stack.size() << endl;
					// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
					Pointer_Base2 = Stack_Base2.at(0);
					if (Pointer_Base2.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
						Stack_Base2.erase(Stack_Base2.begin());
					}
					//++decoding_info.Counter;
					for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
						//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
						//cout << "C:" << Stack.size() << endl;
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
							if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
							//cout << "!";
							break;
						}
						// End
						Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
						Child_Node_Base2.base = Sorted_Base_2;
						if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
							++Child_Node_Base2.D_z;
							Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
						}
						++decoding_info.STE;
						++decoding_info.Binary_STE;
						// Child_Node: The node we are examining now
						//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
						// Reach level k
						if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
							//cout << "(A1)";
							// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
							if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
							// End
							++decoding_info.Counter;
							codeword_seq = MRIP_codeword_Base2;

							// DM-I: Reach Control level to check hamming distance
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							error_counter = Child_Node_Base2.D_z;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
							}
							decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
							if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
								++decoding_info.DM_STE;
								//cout << decoding_info.DM_STE <<" ";
								continue;
							}
							// if DM-I condition did not fit, then continue
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							for (__int16 j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							decoding_info.STE += (codeword_length - message_length);
							//++decoding_info.CandidateCodeWord;
							//cout << "a3";
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						}
						// Did not reach level k, but reach i errors (compared with hard decision result)
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
							//cout << "(A2)";
							//cout << "B";
							++decoding_info.Counter;
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}

							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
								for (size_t j(message_length); j < codeword_length; ++j) {
									if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
										Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
										if (Child_Node_Base2.metric > Best_Goal.metric) break;
									}
								}
								//cout << Child_Node.metric << endl;
								if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
							}
							else Child_Node_Base2 = temp_Node;
							//cout << Adaptive_info.Best_Goal.metric << endl;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
							//cout << Adaptive_info.Best_Goal.metric << endl;
						}
						// Neither reach level k nor reach i error
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
							//cout << "(A3)";
							if (Child_Node_Base2.metric != Pointer_Base2.metric)
								Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
							else {
								Stack_Base2.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
						}
						if (Best_Goal.metric < OSC_metric_thr) {
							decoding_info.DoubleDecoder = FALSE;
							//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
							break;
						}
					}
				} while (!Stack_Base2.empty());
				if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
			}
		}
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;

	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

}

void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number_Base1(0),
		position_number_Base2(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Initial_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);

	
	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr (0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	bool First_Flag_Base1 = TRUE;
	bool First_Flag_Base2 = TRUE;
	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base1 = Stack_Base1.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base1.erase(Stack_Base1.begin());
			
			/*
			if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
				Stack_Base2.erase(Stack_Base2.begin());
			}
			*/
			//可能會刪到不該刪的Node 在第二次傳回來的時候
			
			//++decoding_info.Counter;
			if ((Pointer_Base1.level < message_length) && (Pointer_Base1.metric < Best_Goal.metric) && (Pointer_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Pointer_Base1.level); j < message_length; ++j) {
					Pointer_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Pointer_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Pointer_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base1.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base1.metric)Pointer_Base1 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base1, Best_Goal, Stack_Base1);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base2 = Child_Node_Base1;
						Child_Node_Base2.base = Sorted_Base_2;
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base1.metric)
							Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base2.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
								if (First_Flag_Base2 == TRUE) {
									Stack_CBC_Base2.at(0) = Child_Node_Base2;
									First_Flag_Base2 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
								}
								position_number_Base2++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base2.metric)
								Child_Node_Base2 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
						
						if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
							if (Stack_Base2.at(0).level == 0) {
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
							else {
								Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
							}
							Child_Node_Base1.base = Sorted_Base_1;
						}
						
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		}
		else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
			Next_Flag = Sorted_Base_2;
		}
		else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			
			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (First_Flag_Base2 == TRUE) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
						First_Flag_Base2 = FALSE;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Pointer_Base2);
					}
					position_number_Base2++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base2.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
						}
						*/
						Child_Node_Base1 = Child_Node_Base2;
						Child_Node_Base1.base = Sorted_Base_1;
						Child_Node_Base2.base = Sorted_Base_2;
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//

						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
								if (First_Flag_Base1 == TRUE) {
									Stack_CBC_Base1.at(0) = Child_Node_Base1;
									First_Flag_Base1 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
								}
								position_number_Base1++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
								Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base1;
							for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base1.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_1;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base1.metric)
								Child_Node_Base1 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
						
						if (Child_Node_Base2.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							if (Stack_Base1.empty())Stack_Base1.insert(Stack_Base1.begin(), Initial_Node);
							if (Stack_Base1.at(0).level == 0) {
								Stack_Base1.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
							else {
								Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
							}
							Child_Node_Base2.base = Sorted_Base_2;

						}
						
					}
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
		}
		else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
			Next_Flag = Sorted_Base_1;
		}
	} while (!Stack_Base1.empty() || !Stack_Base2.empty());

	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		if (decoding_info.Cancelled_Candidate_i != 0) {
			Operater_Deletion = TRUE;
			Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
			Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
		}
		Stack_Base1 = Stack_CBC_Base1;
		Stack_Base2 = Stack_CBC_Base2;
		Stack_CBC_Base1.clear();
		Stack_CBC_Base2.clear();
		Stack_CBC_Base1.insert(Stack_CBC_Base1.begin(), Initial_Node);
		Stack_CBC_Base2.insert(Stack_CBC_Base2.begin(), Initial_Node);
		First_Flag_Base1 = TRUE;
		First_Flag_Base2 = TRUE;
		position_number_Base1 = 0;
		position_number_Base2 = 0;
		do {
			if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base1.erase(Stack_Base1.begin());
				/*
				if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
					Stack_Base2.erase(Stack_Base2.begin());
				}
				//可能會刪到不該刪的Node 在第二次傳回來的時候
				*/
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						

						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Child_Node_Base1.metric != Pointer_Base1.metric)
								Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
							else {
								//Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
						}
						*/
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
				Next_Flag = Sorted_Base_2;
			}
			else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base2 = Stack_Base2.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base2.erase(Stack_Base2.begin());

				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
						if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						//
						if (Child_Node_Base2.base == 2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}

						
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						Next_Flag = Sorted_Base_1;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
						}
						else Child_Node_Base2 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						//cout << Adaptive_info.Best_Goal.metric << endl;

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
				Next_Flag = Sorted_Base_1;
			}
		} while (!Stack_Base1.empty() || !Stack_Base2.empty());
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;

	
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;


}

void A_star_PC_out_CBC_OSC_WorstMetric(MATRIX<__int8> &G, DECODING_INFO &decoding_info)  
{
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		error_counter(0);
	vector <size_t>
		Location_Index(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq(message_length, 0),
		Hard_RX(codeword_length, 0),
		MRIP_codeword(codeword_length, 0);
	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH Best_Goal(message_length);

	NODE_PATH
		Pointer(message_length),
		Child_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序
	//system("pause");

	double OSC_metric_thr(0);
	double metric_total(0); //PoHan
	double Worst_Metric(0);
	bool Update_Flag;
	size_t Update_Num(0), Non_Update_Num(0);
	decoding_info.First_nonzero_metric = 0;
	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table._matrix[0][i] != 0) Hard_RX.at(i) = 1;
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	metric_total = OSC_metric_thr;
	//decoding_info.Ave_LLR += ((OSC_metric_thr / codeword_length) * 2 / decoding_info.var);

	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold
	Worst_Metric = decoding_info.Worst_Metric_Ratio * metric_total;

	for (size_t i(codeword_length - decoding_info.Number_of_the_last_symbols); i < codeword_length; ++i) {
		OSC_metric_thr += Metric_Table._matrix[0][i] + Metric_Table._matrix[1][i];
	}
	//decoding_info.DoubleDecoder = TRUE;
	// 開始 Tree Search
	//vector<__int8> BestGoalTemp;
	//decoding_info.code_seq.erase(decoding_info.code_seq.begin() + (decoding_info.code_seq.size() / 2), decoding_info.code_seq.end());

	do {

		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
		//}
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer = Stack.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack.erase(Stack.begin());
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Update_Flag = FALSE;
			Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
			if (new_bit != Hard_RX.at(Pointer.level)) {
				++Child_Node.D_z;
				Child_Node.Diff_Index.push_back(Pointer.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now
			if (Child_Node.metric > Worst_Metric) continue;

			// Reach level k
			if ((Child_Node.level == message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword;

				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack, decoding_info);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
					if (Best_Goal.metric >= OSC_metric_thr) {
						decoding_info.num_Best_Node_Update += decoding_info.temp_num_Deleted_Candidate_in_Stack;
					}
				}

				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Child_Node.level); j < message_length; ++j) {
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G,
							Metric_Table,
							codeword_seq,
							Hard_RX,
							Child_Node,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;

				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;

				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX.at(j)) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node.metric)
					Child_Node = temp_Node;
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Update_Flag == TRUE) {
					decoding_info.num_Best_Node_Update++;
				}
				//cout << "p";
				if (Best_Goal.metric == Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Child_Node.metric;
				}
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node.level < message_length) && (Child_Node.metric < Best_Goal.metric) && (Child_Node.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node.metric != Pointer.metric)
					Place_Node(Stack, Child_Node, decoding_info);
				else {
					Place_Node(Stack, Child_Node, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
			}
		}
		if (Best_Goal.metric < OSC_metric_thr) break;

	} while (!Stack.empty());
	//cout << "Counter: " << decoding_info.Counter << endl;

	//PoHan
	if ((Best_Goal.metric / metric_total) > decoding_info.Worst_OSC_Ratio) {
		decoding_info.Worst_OSC_Ratio = (Best_Goal.metric / metric_total);
	}
	//end
	decoding_info.TotalCounter += decoding_info.Counter;
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);

	//cout << decoding_info.code_seq.size() << "," << decoding_info.Sorted_R.size() << endl;
	//system("pause");


	//Error Test
	/*
	for (int i = 0; i < 192; ++i) {
		//cout << decoding_info.Sorted_R.at(i) << ",";
		if (decoding_info.Sorted_R.at(i) < 0 && decoding_info.code_seq.at(i) == 0) decoding_info.Error_Accumulation.at(i)++;
		else if (decoding_info.Sorted_R.at(i) > 0 && decoding_info.code_seq.at(i) == 1) decoding_info.Error_Accumulation.at(i)++;
	}*/
	//cout << endl;

	/*
	int Total_Error = 0;
	int MRIP_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size() / 2; ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++MRIP_Error;
			//cout << "index: " << i << ", ";
		}
	}
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
			++Total_Error;
			//cout << "index: " << i << ", ";
		}
	}
	//cout << MRIP_Error << endl;
	*/
	/*
	int Total_Error = 0;
	for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
		if (decoding_info.code_seq.at(i) != codeword_seq.at(i)) {
			++Total_Error;
			cout << "index: (" << i << ", " << decoding_info.Sorted_R.at(i) << ") / ";
		}
	}*/

	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);

	/*
	int Error = 0;
	for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) ++Error;
	}
	if (Error != 0) {
		cout << "MRIP Error: " << MRIP_Error << endl;
		cout << "Wrong Index: ";
		for (int i = 0; i < decoding_info.code_seq.size(); ++i) {
			if (decoding_info.code_seq.at(i) != Hard_RX.at(i)) {
				cout << i << ", " << decoding_info.rx_signal_seq.at(i) << " | ";
			}
		}
		cout << endl;
		cout << "Total Error: " << Error << endl;
		//cout << "Error Index: ";
		//for (int i = 0; i < decoding_info.message_seq.size(); ++i) {
			//if (decoding_info.estimated_codeword.at(i) != decoding_info.message_seq.at(i)) cout << "(" << i << "), ";
		//}
		cout << endl <<endl;
	}*/

	//cout << endl << "Check: ";
	//system("pause");
	//cout << "a";


	//cout << "b";
	if (decoding_info.First_nonzero_metric == 0) {
		decoding_info.First_nonzero_metric = 1.0;
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;
}

void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano_Sufficient_Condition(MATRIX<__int8> &G, DECODING_INFO &decoding_info) { //除了第一層的multibase 其他層都是各自長各自的
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length),
		Fano_Metric_Table_Base1(2, message_length),
		Fano_Metric_Table_Base2(2, message_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length),
		Initial_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);


	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer_Base1 = Stack_Base1.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack_Base1.erase(Stack_Base1.begin());
		/*
		if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
			Stack_Base2.erase(Stack_Base2.begin());
		}
		*/
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
			Child_Node_Base1.base = Sorted_Base_1;
			if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
				++Child_Node_Base1.D_z;
				Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node_Base1.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;

				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				/*
				if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					Child_Node_Base1.base = Sorted_Base_1;
				}
				*/
				if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base1.at(0) = Child_Node_Base1;
					}
					else {
						//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
						Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
					Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_Base1.metric)
					Child_Node_Base1 = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node_Base1.metric != Pointer_Base1.metric)
					Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
				else {
					Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
				if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
					if (Stack_Base2.at(0).level == 0) {
						Stack_Base2.at(0) = Child_Node_Base1;
						++decoding_info.COM;
					}
					else {
						Place_Node_Fano(Stack_Base2, Child_Node_Base1, decoding_info);
					}
					Child_Node_Base1.base = Sorted_Base_1;
				}

			}
		}
		if (Best_Goal.metric < OSC_metric_thr) {
			decoding_info.DoubleDecoder = FALSE;
			break;
		}
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
			//}
	} while (!Stack_Base1.empty()); // phase 1
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		position_number = 0;
		do {
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			//++decoding_info.Counter;
			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Pointer_Base2);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;

						//
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
							}
							position_number++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}

						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		} while (!Stack_Base2.empty());
		if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
			Stack_Base1 = Stack_CBC_Base1;
			Stack_CBC_Base1.clear();
			/***開始decode****/
			position_number = 0;
			Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
			if (decoding_info.Cancelled_Candidate_i != 0) {
				Operater_Deletion = TRUE;
				Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
				Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
			}
			// 開始 Tree Search
			do {
				//cout << Stack.size() << endl;
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				if (Pointer_Base1.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack_Base1.erase(Stack_Base1.begin());
				}

				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
					//cout << "C:" << Stack.size() << endl;
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now
					//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "(A1)";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;

						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "(A2)";
						//cout << "B";
						++decoding_info.Counter;

						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
							}
							position_number++;
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
						//cout << Adaptive_info.Best_Goal.metric << endl;
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "(A3)";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Stack_Base1.at(0) = Child_Node_Base1;
							++decoding_info.COM;
						}

					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
					break;
				}
			} while (!Stack_Base1.empty()); // phase 2
			if (decoding_info.DoubleDecoder == TRUE) {
				Stack_Base2 = Stack_CBC_Base2;
				Stack_CBC_Base2.clear();
				/***開始decode****/
				position_number = 0;
				Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
				do {
					//cout << Stack.size() << endl;
					// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
					Pointer_Base2 = Stack_Base2.at(0);
					if (Pointer_Base2.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
						Stack_Base2.erase(Stack_Base2.begin());
					}
					//++decoding_info.Counter;
					for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
						//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
						//cout << "C:" << Stack.size() << endl;
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
							if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
							//cout << "!";
							break;
						}
						// End
						Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
						Child_Node_Base2.base = Sorted_Base_2;
						if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
							++Child_Node_Base2.D_z;
							Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
						}
						++decoding_info.STE;
						++decoding_info.Binary_STE;
						// Child_Node: The node we are examining now
						//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
						// Reach level k
						if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
							//cout << "(A1)";
							// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
							if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
							// End
							++decoding_info.Counter;
							codeword_seq = MRIP_codeword_Base2;

							// DM-I: Reach Control level to check hamming distance
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							error_counter = Child_Node_Base2.D_z;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
							}
							decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
							if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
								++decoding_info.DM_STE;
								//cout << decoding_info.DM_STE <<" ";
								continue;
							}
							// if DM-I condition did not fit, then continue
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							for (__int16 j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							decoding_info.STE += (codeword_length - message_length);
							//++decoding_info.CandidateCodeWord;
							//cout << "a3";
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						}
						// Did not reach level k, but reach i errors (compared with hard decision result)
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
							//cout << "(A2)";
							//cout << "B";
							++decoding_info.Counter;
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}

							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
								for (size_t j(message_length); j < codeword_length; ++j) {
									if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
										Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
										if (Child_Node_Base2.metric > Best_Goal.metric) break;
									}
								}
								//cout << Child_Node.metric << endl;
								if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
							}
							else Child_Node_Base2 = temp_Node;
							//cout << Adaptive_info.Best_Goal.metric << endl;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
							//cout << Adaptive_info.Best_Goal.metric << endl;
						}
						// Neither reach level k nor reach i error
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
							//cout << "(A3)";
							if (Child_Node_Base2.metric != Pointer_Base2.metric)
								Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							else {
								Stack_Base2.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
						}
						if (Best_Goal.metric < OSC_metric_thr) {
							decoding_info.DoubleDecoder = FALSE;
							//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
							break;
						}
					}
				} while (!Stack_Base2.empty());
				if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
			}
		}
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;

	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

}

void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano_Sufficient_Condition(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number_Base1(0),
		position_number_Base2(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length),
		Fano_Metric_Table_Base1(2, message_length),
		Fano_Metric_Table_Base2(2, message_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Initial_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);


	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	bool First_Flag_Base1 = TRUE;
	bool First_Flag_Base2 = TRUE;
	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base1 = Stack_Base1.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base1.erase(Stack_Base1.begin());

			/*
			if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
				Stack_Base2.erase(Stack_Base2.begin());
			}
			*/
			//可能會刪到不該刪的Node 在第二次傳回來的時候

			//++decoding_info.Counter;
			if ((Pointer_Base1.level < message_length) && (Pointer_Base1.metric < Best_Goal.metric) && (Pointer_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Pointer_Base1.level); j < message_length; ++j) {
					Pointer_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Pointer_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Pointer_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base1.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base1.metric)Pointer_Base1 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base1, Best_Goal, Stack_Base1);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base2 = Child_Node_Base1;
						Child_Node_Base2.base = Sorted_Base_2;
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base1.metric)
							Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base2.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
								if (First_Flag_Base2 == TRUE) {
									Stack_CBC_Base2.at(0) = Child_Node_Base2;
									First_Flag_Base2 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
								}
								position_number_Base2++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base2.metric)
								Child_Node_Base2 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}

						if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
							if (Stack_Base2.at(0).level == 0) {
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
							else {
								Place_Node_Fano(Stack_Base2, Child_Node_Base1, decoding_info);
							}
							Child_Node_Base1.base = Sorted_Base_1;
						}

					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		}
		else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
			Next_Flag = Sorted_Base_2;
		}
		else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());

			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (First_Flag_Base2 == TRUE) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
						First_Flag_Base2 = FALSE;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Pointer_Base2);
					}
					position_number_Base2++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base2.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
						}
						*/
						Child_Node_Base1 = Child_Node_Base2;
						Child_Node_Base1.base = Sorted_Base_1;
						Child_Node_Base2.base = Sorted_Base_2;
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//

						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
								if (First_Flag_Base1 == TRUE) {
									Stack_CBC_Base1.at(0) = Child_Node_Base1;
									First_Flag_Base1 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
								}
								position_number_Base1++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
								Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base1;
							for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base1.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_1;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base1.metric)
								Child_Node_Base1 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}

						if (Child_Node_Base2.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							if (Stack_Base1.empty())Stack_Base1.insert(Stack_Base1.begin(), Initial_Node);
							if (Stack_Base1.at(0).level == 0) {
								Stack_Base1.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
							else {
								Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
							}
							Child_Node_Base2.base = Sorted_Base_2;

						}

					}
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
		}
		else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
			Next_Flag = Sorted_Base_1;
		}
	} while (!Stack_Base1.empty() || !Stack_Base2.empty());

	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		if (decoding_info.Cancelled_Candidate_i != 0) {
			Operater_Deletion = TRUE;
			Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
			Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
		}
		Stack_Base1 = Stack_CBC_Base1;
		Stack_Base2 = Stack_CBC_Base2;
		Stack_CBC_Base1.clear();
		Stack_CBC_Base2.clear();
		Stack_CBC_Base1.insert(Stack_CBC_Base1.begin(), Initial_Node);
		Stack_CBC_Base2.insert(Stack_CBC_Base2.begin(), Initial_Node);
		First_Flag_Base1 = TRUE;
		First_Flag_Base2 = TRUE;
		position_number_Base1 = 0;
		position_number_Base2 = 0;
		do {
			if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base1.erase(Stack_Base1.begin());
				/*
				if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
					Stack_Base2.erase(Stack_Base2.begin());
				}
				//可能會刪到不該刪的Node 在第二次傳回來的時候
				*/
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}


						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Child_Node_Base1.metric != Pointer_Base1.metric)
								Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
							else {
								//Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
						}
						*/
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
				Next_Flag = Sorted_Base_2;
			}
			else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base2 = Stack_Base2.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base2.erase(Stack_Base2.begin());

				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
						if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						//
						if (Child_Node_Base2.base == 2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}


						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						Next_Flag = Sorted_Base_1;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
						}
						else Child_Node_Base2 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						//cout << Adaptive_info.Best_Goal.metric << endl;

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
				Next_Flag = Sorted_Base_1;
			}
		} while (!Stack_Base1.empty() || !Stack_Base2.empty());
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;


	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;


}

void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info) { //除了第一層的multibase 其他層都是各自長各自的
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length),
		Fano_Metric_Table_Base1(2, message_length),
		Fano_Metric_Table_Base2(2, message_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Pointer_CBC(message_length),
		Child_Node_CBC(message_length),
		temp_Node(message_length),
		Initial_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);


	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除

	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
		Pointer_Base1 = Stack_Base1.at(0);
		/*
		if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
			Stack.erase(Stack.begin());
		}
		*/
		Stack_Base1.erase(Stack_Base1.begin());
		/*
		if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
			Stack_Base2.erase(Stack_Base2.begin());
		}
		*/
		//++decoding_info.Counter;
		for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
			Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
			Child_Node_Base1.base = Sorted_Base_1;
			if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
				++Child_Node_Base1.D_z;
				Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
			}
			++decoding_info.STE;
			++decoding_info.Binary_STE;
			// Child_Node: The node we are examining now

			// Reach level k
			if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
				//cout << "A";
				++decoding_info.Counter;
				codeword_seq = MRIP_codeword_Base1;
				// DM-I: Reach Control level to check hamming distance

				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
					}
				}
				error_counter = Child_Node_Base1.D_z;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
				}
				decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
				if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
					++decoding_info.DM_STE;
					//cout << decoding_info.DM_STE <<" ";
					continue;
				}
				// if DM-I condition did not fit, then continue
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
						//cout << "o";
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				for (__int16 j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}
				decoding_info.STE += (codeword_length - message_length);
				//++decoding_info.CandidateCodeWord;
				//cout << "a3";
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Did not reach level k, but reach i errors (compared with hard decision result)
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "B";
				++decoding_info.Counter;

				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
				/*
				if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
					Child_Node_Base1.base = Sorted_Base_1;
				}
				*/
				if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base1.at(0) = Child_Node_Base1;
					}
					else {
						//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
						Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
					Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Child_Node_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Child_Node_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				temp_Node.base = Sorted_Base_1;
				//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Child_Node_Base1.metric > Best_Goal.metric) break;
					}
				}

				if (temp_Node.metric < Child_Node_Base1.metric)
					Child_Node_Base1 = temp_Node;
				Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

				//cout << "p";
				/*
				if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
					decoding_info.First_nonzero_metric = Best_Goal.metric;
				}
				else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
					decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
				}
				*/
			}
			// Neither reach level k nor reach i error
			else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
				//cout << "C";
				if (Child_Node_Base1.metric != Pointer_Base1.metric)
					Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
				else {
					Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
					//Stack.at(0) = Child_Node;
					//++decoding_info.COM;
				}
				if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
					Child_Node_Base1.base = Sorted_Base_2;
					if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
					if (Stack_Base2.at(0).level == 0) {
						Stack_Base2.at(0) = Child_Node_Base1;
						++decoding_info.COM;
					}
					else {
						Place_Node_Fano(Stack_Base2, Child_Node_Base1, decoding_info);
					}
					Child_Node_Base1.base = Sorted_Base_1;
				}

			}
		}
		if (Best_Goal.metric < OSC_metric_thr) {
			decoding_info.DoubleDecoder = FALSE;
			break;
		}
		//if (BestGoalTemp != Best_Goal.message_bits) {
			//BestGoalTemp = Best_Goal.message_bits;
			//if (BestGoalTemp == decoding_info.code_seq) break;
			//}
	} while (!Stack_Base1.empty()); // phase 1
	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		position_number = 0;
		do {
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());
			//++decoding_info.Counter;
			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (position_number == 0) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Pointer_Base2);
					}
					position_number++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;

						//
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
							}
							position_number++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}

						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		} while (!Stack_Base2.empty());
		if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
		if ((decoding_info.DoubleDecoder == TRUE) && (Adaptive_i_Decoder2_i != 0)) {
			decoding_info.Constraint_i = Temp_i;
			decoding_info.CBC_FlippingBit = Temp_CBC;
			Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
			while ((Minus_i--) != 0) {
				if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
				else --decoding_info.CBC_FlippingBit;
			}
			decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
			Stack_Base1 = Stack_CBC_Base1;
			Stack_CBC_Base1.clear();
			/***開始decode****/
			position_number = 0;
			Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
			if (decoding_info.Cancelled_Candidate_i != 0) {
				Operater_Deletion = TRUE;
				Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
				Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
			}
			// 開始 Tree Search
			do {
				//cout << Stack.size() << endl;
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				if (Pointer_Base1.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack_Base1.erase(Stack_Base1.begin());
				}

				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
					//cout << "C:" << Stack.size() << endl;
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now
					//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "(A1)";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;

						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "(A2)";
						//cout << "B";
						++decoding_info.Counter;

						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (position_number == 0) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number, Child_Node_Base1);
							}
							position_number++;
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1);
						//cout << Adaptive_info.Best_Goal.metric << endl;
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "(A3)";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Stack_Base1.at(0) = Child_Node_Base1;
							++decoding_info.COM;
						}

					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
					break;
				}
			} while (!Stack_Base1.empty()); // phase 2
			if (decoding_info.DoubleDecoder == TRUE) {
				Stack_Base2 = Stack_CBC_Base2;
				Stack_CBC_Base2.clear();
				/***開始decode****/
				position_number = 0;
				Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number, Child_Node_Base2);
				do {
					//cout << Stack.size() << endl;
					// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
					Pointer_Base2 = Stack_Base2.at(0);
					if (Pointer_Base2.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
						Stack_Base2.erase(Stack_Base2.begin());
					}
					//++decoding_info.Counter;
					for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
						//cout << "B: " << Adaptive_info.Best_Goal.metric << endl;
						//cout << "C:" << Stack.size() << endl;
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
							if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
							//cout << "!";
							break;
						}
						// End
						Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
						Child_Node_Base2.base = Sorted_Base_2;
						if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
							++Child_Node_Base2.D_z;
							Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
						}
						++decoding_info.STE;
						++decoding_info.Binary_STE;
						// Child_Node: The node we are examining now
						//cout << Child_Node.level << "," << Child_Node.metric << ","<<Child_Node.D_z << endl;
						// Reach level k
						if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
							//cout << "(A1)";
							// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
							if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
							// End
							++decoding_info.Counter;
							codeword_seq = MRIP_codeword_Base2;

							// DM-I: Reach Control level to check hamming distance
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							error_counter = Child_Node_Base2.D_z;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
							}
							decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
							if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
								++decoding_info.DM_STE;
								//cout << decoding_info.DM_STE <<" ";
								continue;
							}
							// if DM-I condition did not fit, then continue
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
									//cout << "o";
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							for (__int16 j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							decoding_info.STE += (codeword_length - message_length);
							//++decoding_info.CandidateCodeWord;
							//cout << "a3";
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						}
						// Did not reach level k, but reach i errors (compared with hard decision result)
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
							//cout << "(A2)";
							//cout << "B";
							++decoding_info.Counter;
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}

							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
								for (size_t j(message_length); j < codeword_length; ++j) {
									if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
										Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
										if (Child_Node_Base2.metric > Best_Goal.metric) break;
									}
								}
								//cout << Child_Node.metric << endl;
								if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
							}
							else Child_Node_Base2 = temp_Node;
							//cout << Adaptive_info.Best_Goal.metric << endl;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
							//cout << Adaptive_info.Best_Goal.metric << endl;
						}
						// Neither reach level k nor reach i error
						else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
							//cout << "(A3)";
							if (Child_Node_Base2.metric != Pointer_Base2.metric)
								Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							else {
								Stack_Base2.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
						}
						if (Best_Goal.metric < OSC_metric_thr) {
							decoding_info.DoubleDecoder = FALSE;
							//cout << Adaptive_info.Best_Goal.metric <<"," << Adaptive_info.OSC_metric_thr <<endl;
							break;
						}
					}
				} while (!Stack_Base2.empty());
				if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
			}
		}
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;

	decoding_info.TotalCounter += decoding_info.Counter;
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;

}

void A_star_2_Base_PC_out_CBC_OSC_Adaptive_i_Parallel_Fano(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	size_t
		message_length(G.Row_number),
		codeword_length(G.Col_number),
		position_number_Base1(0),
		position_number_Base2(0),
		error_counter(0);
	vector <size_t>
		Location_Index_Base1(G.Col_number, 0),
		Location_Index_Base2(G.Col_number, 0);
	vector<__int8>
		codeword_seq(codeword_length, 0),
		message_seq_Base1(message_length, 0),
		Hard_RX_Base1(codeword_length, 0),
		MRIP_codeword_Base1(codeword_length, 0),
		message_seq_Base2(message_length, 0),
		Hard_RX_Base2(codeword_length, 0),
		MRIP_codeword_Base2(codeword_length, 0);
	MATRIX<__int8>
		Sorted_G_Base1(G),
		Sorted_G_Base2(G);
	MATRIX<double>
		Metric_Table_Base1(2, codeword_length),
		Metric_Table_Base2(2, codeword_length),
		Fano_Metric_Table_Base1(2, message_length),
		Fano_Metric_Table_Base2(2, message_length);

	NODE_PATH Best_Goal(message_length);
	Best_Goal.metric = FLT_MAX;

	NODE_PATH
		Pointer_Base1(message_length),
		Child_Node_Base1(message_length),
		Pointer_Base2(message_length),
		Child_Node_Base2(message_length),
		Initial_Node(message_length),
		temp_Node(message_length);
	vector<NODE_PATH> Stack_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_Base2(1, Pointer_Base2);
	vector<NODE_PATH> Stack_CBC_Base1(1, Pointer_Base1);
	vector<NODE_PATH> Stack_CBC_Base2(1, Pointer_Base2);


	decoding_info.Counter = 0;

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G_Base1, Location_Index_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, decoding_info);
	Pre_Procedure_MultiBase(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, 2, 2, decoding_info);	// sorting_rx_signal_seq 為sorting rx結果
	// Location_index 紀錄排序

	double OSC_metric_thr(0);

	for (size_t i(0); i < codeword_length; ++i) {
		// for MRIP constraint
		if (Metric_Table_Base1._matrix[0][i] != 0) Hard_RX_Base1.at(i) = 1;
		if (Metric_Table_Base2._matrix[0][i] != 0) Hard_RX_Base2.at(i) = 1;

		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
		//cout << 
	}
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);  // MRIP_codeword: MRIP message sequence所算出的codeword
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);  // MRIP_codeword: MRIP message sequence所算出的codeword
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;  // 算出OSC threshold

	size_t Next_Flag = Sorted_Base_1;
	decoding_info.DoubleDecoder = TRUE;
	int Level_k_previous, Difference;
	bool Operater_Deletion = FALSE; //第二次的tree search把之前搜尋過的刪除
	bool First_Flag_Base1 = TRUE;
	bool First_Flag_Base2 = TRUE;
	//cout << "A";
	// Decoder(i'-2) -> Early Termination -> Decoder(i')

	size_t Temp_i = decoding_info.Constraint_i, Temp_CBC = decoding_info.CBC_FlippingBit;    // 紀錄一開始的i, CBC
	int Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder1_i;
	while ((Minus_i--) != 0) {
		if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
		else --decoding_info.CBC_FlippingBit;
	}
	//cout << "(X1): " << decoding_info.Constraint_i << "," << decoding_info.CBC_FlippingBit << endl;
	decoding_info.Cancelled_Candidate_i = 0;
	do {
		if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base1 = Stack_Base1.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base1.erase(Stack_Base1.begin());

			/*
			if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
				Stack_Base2.erase(Stack_Base2.begin());
			}
			*/
			//可能會刪到不該刪的Node 在第二次傳回來的時候

			//++decoding_info.Counter;
			if ((Pointer_Base1.level < message_length) && (Pointer_Base1.metric < Best_Goal.metric) && (Pointer_Base1.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				for (__int16 j(Pointer_Base1.level); j < message_length; ++j) {
					Pointer_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base1;
				for (size_t index(0); index < Pointer_Base1.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base1.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Pointer_Base1.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base1.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base1,
							Metric_Table_Base1,
							codeword_seq,
							Hard_RX_Base1,
							Pointer_Base1,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
						Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base1.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base1.metric)Pointer_Base1 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base1, Best_Goal, Stack_Base1);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base2 = Child_Node_Base1;
						Child_Node_Base2.base = Sorted_Base_2;
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}
						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base1.metric)
							Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base2.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
								if (First_Flag_Base2 == TRUE) {
									Stack_CBC_Base2.at(0) = Child_Node_Base2;
									First_Flag_Base2 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
								}
								position_number_Base2++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
								Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base2;
							for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base2.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base2,
										Metric_Table_Base2,
										codeword_seq,
										Hard_RX_Base2,
										Child_Node_Base2,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_2;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base2.metric)
								Child_Node_Base2 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}

						if (Child_Node_Base1.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Stack_Base2.empty())Stack_Base2.insert(Stack_Base2.begin(), Initial_Node);
							if (Stack_Base2.at(0).level == 0) {
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
							else {
								Place_Node_Fano(Stack_Base2, Child_Node_Base1, decoding_info);
							}
							Child_Node_Base1.base = Sorted_Base_1;
						}

					}
				}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
			//if (BestGoalTemp != Best_Goal.message_bits) {
				//BestGoalTemp = Best_Goal.message_bits;
				//if (BestGoalTemp == decoding_info.code_seq) break;
				//}
		}
		else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
			Next_Flag = Sorted_Base_2;
		}
		else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
			// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
			Pointer_Base2 = Stack_Base2.at(0);
			/*
			if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
				Stack.erase(Stack.begin());
			}
			*/
			Stack_Base2.erase(Stack_Base2.begin());

			if ((Pointer_Base2.level < message_length) && (Pointer_Base2.metric < Best_Goal.metric) && (Pointer_Base2.D_z == decoding_info.Constraint_i)) {
				//cout << "(A2)";
				//cout << "B";
				++decoding_info.Counter;
				Pointer_Base2.base == Sorted_Base_2;
				if (Pointer_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
					if (First_Flag_Base2 == TRUE) {
						Stack_CBC_Base2.at(0) = Pointer_Base2;
						First_Flag_Base2 = FALSE;
					}
					else {
						//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
						Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Pointer_Base2);
					}
					position_number_Base2++;
					//Place_Node(Stack_CBC, Child_Node, decoding_info);
				}
				for (__int16 j(Pointer_Base2.level); j < message_length; ++j) {
					Pointer_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
				}
				//
				//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

				codeword_seq = MRIP_codeword_Base2;
				for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
					codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
					for (__int16 j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
					}
				}
				//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

				decoding_info.STE += (codeword_length - Pointer_Base2.level);
				//++decoding_info.CandidateCodeWord;

				// CBC
				if (decoding_info.CBC_FlippingBit == 1) {
					temp_Node
						= Control_Band_Check_1bit(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 2) {
					temp_Node
						= Control_Band_Check_2bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info);
				}
				else if (decoding_info.CBC_FlippingBit == 3) {
					temp_Node
						= Control_Band_Check_3bits(
							Sorted_G_Base2,
							Metric_Table_Base2,
							codeword_seq,
							Hard_RX_Base2,
							Pointer_Base2,
							Best_Goal,
							decoding_info,
							1);
				}
				else temp_Node.metric = DBL_MAX;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
						Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
						if (Pointer_Base2.metric > Best_Goal.metric) break;
					}
				}
				//cout << Pointer.metric << endl;
				if (temp_Node.metric < Pointer_Base2.metric)Pointer_Base2 = temp_Node;
				//cout << Adaptive_info.Best_Goal.metric << endl;
				Update_Best_Goal_Procedure(Pointer_Base2, Best_Goal, Stack_Base2);
				//cout << Adaptive_info.Best_Goal.metric << endl;
			}
			else {
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base2.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
						}
						*/
						Child_Node_Base1 = Child_Node_Base2;
						Child_Node_Base1.base = Sorted_Base_1;
						Child_Node_Base2.base = Sorted_Base_2;
						if (Child_Node_Base2.base == Sorted_Base_2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}
						//

						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						for (size_t j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}

						if (temp_Node.metric < Child_Node_Base2.metric)
							Child_Node_Base2 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						if ((Best_Goal.metric > OSC_metric_thr) && (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits)) {
							if (Child_Node_Base1.base == Sorted_Base_1) {//為了最後一次的decode作準備
								if (First_Flag_Base1 == TRUE) {
									Stack_CBC_Base1.at(0) = Child_Node_Base1;
									First_Flag_Base1 = FALSE;
								}
								else {
									//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
									Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
								}
								position_number_Base1++;
								//Place_Node(Stack_CBC, Child_Node, decoding_info);
							}
							for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
								Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
							}
							//
							//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);				

							codeword_seq = MRIP_codeword_Base1;
							for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
								codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
								for (__int16 j(message_length); j < codeword_length; ++j) {
									codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
								}
							}
							//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

							decoding_info.STE += (codeword_length - Child_Node_Base1.level);
							//++decoding_info.CandidateCodeWord;

							// CBC
							if (decoding_info.CBC_FlippingBit == 1) {
								temp_Node
									= Control_Band_Check_1bit(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 2) {
								temp_Node
									= Control_Band_Check_2bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info);
							}
							else if (decoding_info.CBC_FlippingBit == 3) {
								temp_Node
									= Control_Band_Check_3bits(
										Sorted_G_Base1,
										Metric_Table_Base1,
										codeword_seq,
										Hard_RX_Base1,
										Child_Node_Base1,
										Best_Goal,
										decoding_info,
										1);
							}
							else temp_Node.metric = DBL_MAX;
							temp_Node.base = Sorted_Base_1;
							//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}

							if (temp_Node.metric < Child_Node_Base1.metric)
								Child_Node_Base1 = temp_Node;
							Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						}
						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}

						if (Child_Node_Base2.level == message_length - Multiple_Basis_Bits) {
							Child_Node_Base2.base = Sorted_Base_1;
							if (Stack_Base1.empty())Stack_Base1.insert(Stack_Base1.begin(), Initial_Node);
							if (Stack_Base1.at(0).level == 0) {
								Stack_Base1.at(0) = Child_Node_Base2;
								++decoding_info.COM;
							}
							else {
								Place_Node(Stack_Base1, Child_Node_Base2, decoding_info);
							}
							Child_Node_Base2.base = Sorted_Base_2;

						}

					}
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			if (Best_Goal.metric < OSC_metric_thr) {
				decoding_info.DoubleDecoder = FALSE;
				break;
			}
		}
		else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
			Next_Flag = Sorted_Base_1;
		}
	} while (!Stack_Base1.empty() || !Stack_Base2.empty());

	if (Best_Goal.metric < OSC_metric_thr*Adaptive_i_Parameter) decoding_info.DoubleDecoder = FALSE;
	if (decoding_info.DoubleDecoder == TRUE) {
		decoding_info.Constraint_i = Temp_i;
		decoding_info.CBC_FlippingBit = Temp_CBC;
		Minus_i = decoding_info.Constraint_i + decoding_info.CBC_FlippingBit - Adaptive_i_Decoder2_i;
		while ((Minus_i--) != 0) {
			if (decoding_info.Constraint_i > decoding_info.CBC_FlippingBit) --decoding_info.Constraint_i;
			else --decoding_info.CBC_FlippingBit;
		}
		decoding_info.Cancelled_Candidate_i = Adaptive_i_Decoder1_i;
		if (decoding_info.Cancelled_Candidate_i != 0) {
			Operater_Deletion = TRUE;
			Level_k_previous = message_length - (decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i);
			Difference = decoding_info.Constraint_i - decoding_info.Cancelled_Candidate_i;
		}
		Stack_Base1 = Stack_CBC_Base1;
		Stack_Base2 = Stack_CBC_Base2;
		Stack_CBC_Base1.clear();
		Stack_CBC_Base2.clear();
		Stack_CBC_Base1.insert(Stack_CBC_Base1.begin(), Initial_Node);
		Stack_CBC_Base2.insert(Stack_CBC_Base2.begin(), Initial_Node);
		First_Flag_Base1 = TRUE;
		First_Flag_Base2 = TRUE;
		position_number_Base1 = 0;
		position_number_Base2 = 0;
		do {
			if ((Next_Flag == Sorted_Base_1) && (!Stack_Base1.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base1 = Stack_Base1.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base1.erase(Stack_Base1.begin());
				/*
				if (Pointer_Base1.level < message_length - Multiple_Basis_Bits) {
					Stack_Base2.erase(Stack_Base2.begin());
				}
				//可能會刪到不該刪的Node 在第二次傳回來的時候
				*/
				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base1.level == Level_k_previous) && (Pointer_Base1.D_z <= Difference)) {
						if (Pointer_Base1.level != (message_length - 1)) Stack_Base1.erase(Stack_Base1.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base1, Child_Node_Base1, Metric_Table_Base1, Fano_Metric_Table_Base1, new_bit);
					Child_Node_Base1.base = Sorted_Base_1;
					if (new_bit != Hard_RX_Base1.at(Pointer_Base1.level)) {
						++Child_Node_Base1.D_z;
						Child_Node_Base1.Diff_Index.push_back(Pointer_Base1.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base1.level == message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base1.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base1;
						Next_Flag = Sorted_Base_2;
						// DM-I: Reach Control level to check hamming distance

						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base1.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
								Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base1.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						/*
						if (Child_Node_Base1.level < message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
						}
						*/
						Child_Node_Base1.base = Sorted_Base_1;
						if (Child_Node_Base1.base == 1) {//為了最後一次的decode作準備
							if (First_Flag_Base1 == TRUE) {
								Stack_CBC_Base1.at(0) = Child_Node_Base1;
								First_Flag_Base1 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base1, Child_Node_Base1, decoding_info);
								Stack_CBC_Base1.insert(Stack_CBC_Base1.begin() + position_number_Base1, Child_Node_Base1);
							}
							position_number_Base1++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base1.level); j < message_length; ++j) {
							Child_Node_Base1.message_bits.at(j) = Hard_RX_Base1.at(j);
						}


						//
						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);

						codeword_seq = MRIP_codeword_Base1;
						for (size_t index(0); index < Child_Node_Base1.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base1.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Child_Node_Base1.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base1.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base1,
									Metric_Table_Base1,
									codeword_seq,
									Hard_RX_Base1,
									Child_Node_Base1,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_1;
						Next_Flag = Sorted_Base_2;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base1.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) {
									Child_Node_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base1.metric > Best_Goal.metric) break;
								}
							}
							if (temp_Node.metric < Child_Node_Base1.metric)Child_Node_Base1 = temp_Node;
						}
						else Child_Node_Base1 = temp_Node;
						Update_Best_Goal_Procedure(Child_Node_Base1, Best_Goal, Stack_Base1, Stack_CBC_Base1);

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base1.level < message_length) && (Child_Node_Base1.metric < Best_Goal.metric) && (Child_Node_Base1.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base1.metric != Pointer_Base1.metric)
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
						else {
							Place_Node_Fano(Stack_Base1, Child_Node_Base1, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
						/*
						if (Child_Node_Base1.level <= message_length - Multiple_Basis_Bits) {
							Child_Node_Base1.base = Sorted_Base_2;
							if (Child_Node_Base1.metric != Pointer_Base1.metric)
								Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
							else {
								//Place_Node(Stack_Base2, Child_Node_Base1, decoding_info);
								Stack_Base2.at(0) = Child_Node_Base1;
								++decoding_info.COM;
							}
						}
						*/
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_1) && (Stack_Base1.empty())) {
				Next_Flag = Sorted_Base_2;
			}
			else if ((Next_Flag == Sorted_Base_2) && (!Stack_Base2.empty())) {
				// 這裡的pointer不是真的pointer, 只是用pointer去等於Stack
				Pointer_Base2 = Stack_Base2.at(0);
				/*
				if (Pointer.level == (message_length - 1)) { // 當pointer的level達到k-1之後, 接下來的兩個child node都會是k, 因此在這個步驟把stack最上面的值給削掉(用pointer存資料之後就能pop掉的概念)
					Stack.erase(Stack.begin());
				}
				*/
				Stack_Base2.erase(Stack_Base2.begin());

				//++decoding_info.Counter;
				for (__int8 new_bit(0); new_bit < 2; ++new_bit) {
					// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
					if ((Operater_Deletion == TRUE) && (Pointer_Base2.level == Level_k_previous) && (Pointer_Base2.D_z <= Difference)) {
						if (Pointer_Base2.level != (message_length - 1)) Stack_Base2.erase(Stack_Base2.begin());
						//cout << "!";
						break;
					}
					// End
					Extend_Node_Procedure_Fano(Pointer_Base2, Child_Node_Base2, Metric_Table_Base2, Fano_Metric_Table_Base2, new_bit);
					Child_Node_Base2.base = Sorted_Base_2;
					if (new_bit != Hard_RX_Base2.at(Pointer_Base2.level)) {
						++Child_Node_Base2.D_z;
						Child_Node_Base2.Diff_Index.push_back(Pointer_Base2.level);
					}
					++decoding_info.STE;
					++decoding_info.Binary_STE;
					// Child_Node: The node we are examining now

					// Reach level k
					if ((Child_Node_Base2.level == message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z <= decoding_info.Constraint_i)) {
						//cout << "A";
						// ***  "和原本的Pcout-CBC-OSC不同的地方"  ***
						if (Operater_Deletion == TRUE && Child_Node_Base2.D_z <= decoding_info.Cancelled_Candidate_i) continue;
						// End
						++decoding_info.Counter;
						codeword_seq = MRIP_codeword_Base2;
						Next_Flag = Sorted_Base_1;
						// DM-I: Reach Control level to check hamming distance
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j]; //這裡只有算到control而已所以沒有用systemetic encoder
							}
						}
						error_counter = Child_Node_Base2.D_z;
						for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
						}
						decoding_info.Binary_STE += (decoding_info.Control_Level - message_length);
						if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
							++decoding_info.DM_STE;
							//cout << decoding_info.DM_STE <<" ";
							continue;
						}
						// if DM-I condition did not fit, then continue
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							//codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
							for (__int16 j(decoding_info.Control_Level); j < codeword_length; ++j) {
								//cout << "o";
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						for (__int16 j(message_length); j < codeword_length; ++j) {
							if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
								Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
								if (Child_Node_Base2.metric > Best_Goal.metric) break;
							}
						}
						decoding_info.STE += (codeword_length - message_length);
						//++decoding_info.CandidateCodeWord;
						//cout << "a3";
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2, Stack_CBC_Base2);
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Did not reach level k, but reach i errors (compared with hard decision result)
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z == decoding_info.Constraint_i)) {
						//cout << "B";
						++decoding_info.Counter;
						//
						if (Child_Node_Base2.base == 2) {//為了最後一次的decode作準備
							if (First_Flag_Base2 == TRUE) {
								Stack_CBC_Base2.at(0) = Child_Node_Base2;
								First_Flag_Base2 = FALSE;
							}
							else {
								//Place_Node(Stack_CBC_Base2, Child_Node_Base2, decoding_info);
								Stack_CBC_Base2.insert(Stack_CBC_Base2.begin() + position_number_Base2, Child_Node_Base2);
							}
							position_number_Base2++;
							//Place_Node(Stack_CBC, Child_Node, decoding_info);
						}
						for (__int16 j(Child_Node_Base2.level); j < message_length; ++j) {
							Child_Node_Base2.message_bits.at(j) = Hard_RX_Base2.at(j);
						}


						//Systematic_Linear_Block_Code_Encoder(Sorted_G, Child_Node.message_bits, codeword_seq);
						codeword_seq = MRIP_codeword_Base2;
						for (size_t index(0); index < Child_Node_Base2.Diff_Index.size(); ++index) {
							codeword_seq.at(Child_Node_Base2.Diff_Index.at(index)) ^= 1;
							for (__int16 j(message_length); j < codeword_length; ++j) {
								codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Child_Node_Base2.Diff_Index.at(index)][j];
							}
						}
						//codeword_seq為到達i的message直接做extension到level k, 並延伸到level n的seqeunce

						decoding_info.STE += (codeword_length - Child_Node_Base2.level);
						//++decoding_info.CandidateCodeWord;

						// CBC
						if (decoding_info.CBC_FlippingBit == 1) {
							temp_Node
								= Control_Band_Check_1bit(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 2) {
							temp_Node
								= Control_Band_Check_2bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info);
						}
						else if (decoding_info.CBC_FlippingBit == 3) {
							temp_Node
								= Control_Band_Check_3bits(
									Sorted_G_Base2,
									Metric_Table_Base2,
									codeword_seq,
									Hard_RX_Base2,
									Child_Node_Base2,
									Best_Goal,
									decoding_info,
									1);
						}
						else temp_Node.metric = DBL_MAX;
						temp_Node.base = Sorted_Base_2;
						Next_Flag = Sorted_Base_1;
						//else cout << endl << "CBC should be equal to or smaller than 3! Please reset the system !" << endl;
						if (Operater_Deletion == FALSE || Child_Node_Base2.D_z > decoding_info.Cancelled_Candidate_i) {
							for (size_t j(message_length); j < codeword_length; ++j) {
								if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) {
									Child_Node_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
									if (Child_Node_Base2.metric > Best_Goal.metric) break;
								}
							}
							//cout << Child_Node.metric << endl;
							if (temp_Node.metric < Child_Node_Base2.metric)Child_Node_Base2 = temp_Node;
						}
						else Child_Node_Base2 = temp_Node;
						//cout << Adaptive_info.Best_Goal.metric << endl;
						Update_Best_Goal_Procedure(Child_Node_Base2, Best_Goal, Stack_Base2);
						//cout << Adaptive_info.Best_Goal.metric << endl;

						//cout << "p";
						/*
						if (Best_Goal.metric == Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {		//PoHan
							decoding_info.First_nonzero_metric = Best_Goal.metric;
						}
						else if (Best_Goal.metric != Sorted_Child_Node.metric && decoding_info.First_nonzero_metric == 0) {
							decoding_info.First_nonzero_metric = Sorted_Child_Node.metric;
						}
						*/
					}
					// Neither reach level k nor reach i error
					else if ((Child_Node_Base2.level < message_length) && (Child_Node_Base2.metric < Best_Goal.metric) && (Child_Node_Base2.D_z < decoding_info.Constraint_i)) {
						//cout << "C";
						if (Child_Node_Base2.metric != Pointer_Base2.metric)
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
						else {
							Place_Node_Fano(Stack_Base2, Child_Node_Base2, decoding_info);
							//Stack.at(0) = Child_Node;
							//++decoding_info.COM;
						}
					}
				}
				if (Best_Goal.metric < OSC_metric_thr) {
					decoding_info.DoubleDecoder = FALSE;
					break;
				}
				//if (BestGoalTemp != Best_Goal.message_bits) {
					//BestGoalTemp = Best_Goal.message_bits;
					//if (BestGoalTemp == decoding_info.code_seq) break;
					//}
			}
			else if ((Next_Flag == Sorted_Base_2) && (Stack_Base2.empty())) {
				Next_Flag = Sorted_Base_1;
			}
		} while (!Stack_Base1.empty() || !Stack_Base2.empty());
	}
	decoding_info.Constraint_i = Temp_i;
	decoding_info.CBC_FlippingBit = Temp_CBC;


	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else if (Best_Goal.base == Sorted_Base_2) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	// BESTONE 
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;

	if (decoding_info.CandidateCodeWord > decoding_info.Worst_Case_Candidate)
		decoding_info.Worst_Case_Candidate = decoding_info.CandidateCodeWord;


}
