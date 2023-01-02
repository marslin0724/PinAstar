#include "AStarDecode.h"

// use Best Graph Algorithm (BGA)
void ABP(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	vector<size_t> permutation_seq(G.Col_number);
	//initial permutation_seq
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	MATRIX<__int8> H,sorted_H,tmp_H;
	tmp_H.Building_Empty_Matrix(G.Col_number-G.Row_number,G.Col_number);
	tmp_H._matrix = G._H;
	sorted_H = tmp_H;
	vector<double> orig_rx(G.Col_number);
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++)
		orig_rx.at(i) = decoding_info.rx_signal_seq.at(i) * (2 / decoding_info.var);
	/*
	for (int i = 0; i < G.Col_number; i++) {
		cout << orig_rx.at(i) << ",";
	}
	cout << endl;*/
	
	ABP_Permutation(decoding_info.rx_signal_seq, tmp_H, sorted_H, permutation_seq);
	//
	// Change into Channel LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) 
		decoding_info.rx_signal_seq.at(i) *= (2 / decoding_info.var);
	vector < vector <double> > M, Eji,tmp_M;
	vector <double> Li(tmp_H.Col_number), zi;
	int I_count = 0, Col = tmp_H.Col_number, Row = tmp_H.Row_number;
	double m_temp, temp;
	vector<__int8> codeword_seq(G.Col_number);
	/*
	for (int i = 0; i < Col; i++) {
		cout << decoding_info.rx_signal_seq.at(i) << ",";
	}
	cout << endl;
	for (int i = 0; i < Col; i++) {
		cout << permutation_seq.at(i) << ",";
	}
	cout << endl;*/

	vector<double> vect(decoding_info.rx_signal_seq.size(), 0);   // for template container
	for (int i = 0; i < Col; i++) {
		vect.at(i) = (decoding_info.rx_signal_seq.at(i));
	}
	for (int j = 0; j < Row; j++) {
		M.push_back(vect);
	}
	Li = vect;
	vect.clear();

	//
	
	while (true)
	{
		for (int j = 0; j < Row; j++) {

			for (int i = 0; i < Col; i++) {

				if (sorted_H._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < Col; k++) {

						if (sorted_H._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}

					temp = log((1 + ta) / (1 - ta));
					if (abs(temp) > LLR_MAX) {
						vect.push_back(temp > 0 ? LLR_MAX : (-1)*LLR_MAX);
					}
					else vect.push_back(temp);
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}

		for (int i = 0; i < Col; i++){
			Li[i] = orig_rx.at(permutation_seq.at(i));
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] += Eji[j][i];
			}
		}

		//
		if (I_count >= decoding_info.SPA_I - 1) {
			for (int i = 0; i < Col; ++i) {
				//decoding_info.rx_signal_seq.at(i) = Li.at(i);
				if (Li.at(i) > 0) codeword_seq.at(i) = 0;
				else codeword_seq.at(i) = 1;
			}
			break;
		}
		// End
		I_count++;
		for (int i = 0; i < Row; i++) {
			for (int j = Row; j < Col; j++) {
				if (sorted_H._matrix[i][j] == 1) {
					M[i][j] = Li.at(j) - Eji[i][j];
				}
			}
		}
		
		/*
		cout << "before permut\n";
		for (int i = 0; i < Col; i++) {
			cout << Li.at(i) << ",";
		}
		cout << endl;
		for (int i = 0; i < Col; i++) {
			cout << permutation_seq.at(i) << ",";
		}
		cout << endl;
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				cout << M[i][j]<<",";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		*/
		//permutation
		//tmp_H = sorted_H;
		//tmp_M = M;
		//ABP_iteration_Permutation(Li, tmp_H, sorted_H, permutation_seq,
			//tmp_M, M);
		//
		/*
		cout << "after permut\n";
		for (int i = 0; i < Col; i++) {
			cout << Li.at(i) << ",";
		}
		cout << endl;
		for (int i = 0; i < Col; i++) {
			cout << permutation_seq.at(i) << ",";
		}
		cout << endl;
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				cout << M[i][j] <<",";
			}
			cout << endl;
		}
		cout << endl;
		cout << endl;
		*/
		vect.clear();
		Eji.clear();
		zi.clear();
	}
	Desort_Function(permutation_seq, codeword_seq, decoding_info.estimated_codeword);
	Desort_Function(permutation_seq, Li, decoding_info.rx_signal_seq);
}

//ABP_permutation
// need to initial orig_permutation_seq
void ABP_Permutation(
	vector<double>			&Rx_signal,
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&orig_permutation_seq)
{
	vector<size_t> permutation_seq(G.Col_number),tmp_permutation_seq(G.Col_number);
	vector<double> tmp_Rx_signal_seq(Rx_signal.size(), 0);
	vector<size_t> Row_Sorter(G.Row_number, 0);
	//
	ABP_Determine_Permutation(Rx_signal, G, sorted_G, permutation_seq, Row_Sorter);
	tmp_Rx_signal_seq = Rx_signal;
	Sort_Function(tmp_Rx_signal_seq, permutation_seq, Rx_signal);
	//
	for (size_t i(0); i < permutation_seq.size(); ++i) {
		tmp_permutation_seq.at(i) = orig_permutation_seq.at(permutation_seq.at(i));
	}
	orig_permutation_seq = tmp_permutation_seq;
}

void ABP_Determine_Permutation(
	vector<double>			&Rx_signal,
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&permutation_seq,
	vector<size_t>			&Row_Sorter)
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
		[&](const double& a, const double& b)
	{return (abs(Rx_signal.at(a)) < abs(Rx_signal.at(b))); });

	// Col. exchange based on permutation sequence
	Sort_Matrix_Col(G, permutation_seq, sorted_G);

	// Obtain systematic form
	ABP_GJ_Elimination(sorted_G, Row_Sorter);

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
}

inline void ABP_iteration_Permutation(
	vector<double>			&Rx_signal,
	MATRIX<__int8>			&G,
	MATRIX<__int8>			&sorted_G,
	vector<size_t>			&orig_permutation_seq,
	vector<vector<double>>		&input_matrix,
	vector<vector<double>>		&output_matrix)
{
	vector<size_t> permutation_seq(G.Col_number), tmp_permutation_seq(G.Col_number);
	vector<double> tmp_Rx_signal_seq(Rx_signal.size(), 0);
	vector<size_t> Row_Sorter(G.Row_number, 0);
	//
	ABP_Determine_Permutation(Rx_signal, G, sorted_G, permutation_seq, Row_Sorter);
	//tmp_Rx_signal_seq = Rx_signal;
	//Sort_Function(tmp_Rx_signal_seq, permutation_seq, Rx_signal);
	//
	Sort_Col(input_matrix, permutation_seq, output_matrix);
	input_matrix = output_matrix;
	Sort_Row(input_matrix, Row_Sorter, output_matrix);
	
	for (size_t i(0); i < permutation_seq.size(); ++i) {
		tmp_permutation_seq.at(i) = orig_permutation_seq.at(permutation_seq.at(i));
	}
	orig_permutation_seq = tmp_permutation_seq;
}

inline void Sort_Col(
	vector<vector<double>>		&input_matrix,
	vector<size_t>				&permutation_seq,
	vector<vector<double>>		&output_matrix)
{
	for (size_t i(0); i < permutation_seq.size(); ++i) {
		for (size_t j(0); j < input_matrix.size(); ++j)
			output_matrix[j][i] = input_matrix[j][permutation_seq.at(i)];
	}
}
inline void Sort_Row(
	vector<vector<double>>		&input_matrix,
	vector<size_t>				&permutation_seq,
	vector<vector<double>>		&output_matrix)
{
	for (size_t i(0); i < permutation_seq.size(); ++i) {
		for (size_t j(0); j < input_matrix.at(0).size(); ++j)
			output_matrix[i][j] = input_matrix[permutation_seq.at(i)][j];
	}
}
//catch record Row permutation
void ABP_GJ_Elimination(MATRIX<__int8> &input_matrix, vector<size_t>& Row_Sorter)
{
	size_t
		temp(0),
		current_position_row(0),
		current_position_col(0);

	//宣告動態向量、矩陣
	/***************************************************************/
	MATRIX<__int8> Register(input_matrix);
	vector<size_t>
		Row_Index(input_matrix.Row_number, 0);
	for (size_t i(0); i < input_matrix.Row_number; ++i)
		Row_Sorter.at(i) = i;
	/***************************************************************/


	while (current_position_row < input_matrix.Row_number) {
		current_position_col = 0;
		/***************************************************************/
		do {
			if (Register._matrix[current_position_row][current_position_col] != 0) break;
			else ++current_position_col;
		} while ((current_position_col < input_matrix.Col_number));
		/***************************************************************/
		Row_Index.at(current_position_row) = current_position_col + 1;

		//使用目前列去消去其他列
		if (current_position_col < input_matrix.Col_number) {
			for (size_t j(0); j < input_matrix.Row_number; ++j) {
				//使用第 i 列去加 j 列 , i 不等於 j		
				if (j != current_position_row)
					if (Register._matrix[j][current_position_col] == 1)
						for (size_t k(0); k < input_matrix.Col_number; ++k)
							if (Register._matrix[current_position_row][k] == 1)
								Register._matrix[j][k] ^= 1;
			};
		};
		/**************************************************************/
		++current_position_row;
	};


	//換列
	/*********************************************************************/
	for (size_t i(0); i < input_matrix.Row_number; ++i) {
		for (size_t j(i); j < input_matrix.Row_number; ++j) {
			if (Row_Index[i] > Row_Index[j]) {
				temp = Row_Index[i];
				Row_Index[i] = Row_Index[j];
				Row_Index[j] = temp;

				temp = Row_Sorter[i];
				Row_Sorter[i] = Row_Sorter[j];
				Row_Sorter[j] = temp;
			};
		};
	};
	
	for (size_t i(0); i < input_matrix.Row_number; ++i)
		for (size_t j(0); j < input_matrix.Col_number; ++j)
			input_matrix._matrix[i][j] = Register._matrix[Row_Sorter.at(i)][j];

}

//Brute Force + ABP(BGA) + A*
//only use (8,4)Hamming code
void BF_ABP_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H, G_, G__, H_inner;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	G_._H = G._H;
	//for BF
	G__.Building_Empty_Matrix(H.Col_number - H.Row_number, H.Col_number);
	G__._matrix = {
	{1,0,0,0,1,1,1,0},
	{0,1,0,0,1,1,0,1},
	{0,0,1,0,1,0,1,1},
	{0,0,0,1,0,1,1,1}
	};

	int i, j, start;
	vector<double>
		subRx(H.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number*(H.Col_number - H.Row_number) / H.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);

	for (i = 0; i < G.Col_number / H.Col_number; i++) {
		//assign sub-Block
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * H.Row_number;
		for (; j < H.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - H.Col_number + H.Row_number);
		}
		//Brute-Force SISO
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (H.Row_number);
		for (; j < H.Col_number; j++) {
			Rx_parity.at(start + j - H.Col_number + H.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;
	for (i = 0; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq.at(i) /= (2 / decoding_info.var);
	}
	//ABP
	ABP(G_, decoding_info);
	//original A*
	A_star_I(G_, decoding_info);
	//set parity LLR estimate
	decoding_info.estimated_codeword.resize(G.Col_number);

	for (i = G_.Col_number; i < G.Col_number; i++) {
		if (Rx_parity.at(i - G_.Col_number) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = tmp_Rx;
}

void ABP_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	ABP(G, decoding_info);
	//can exchange
	A_star_I(G, decoding_info);
	//A_star_3_stack(G, decoding_info);
	/*test
	vector<double> tmp_Rx = decoding_info.rx_signal_seq;
	MATRIX<__int8> H;
	H.Building_Empty_Matrix(G.Col_number - G.Row_number, G.Col_number);
	H._matrix = G._H;
	ASPA_subdecoder(H, decoding_info.rx_signal_seq, decoding_info.var,10);
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.rx_signal_seq.at(i) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = tmp_Rx;
	*/
}

void ASPA_subdecoder(MATRIX<__int8> &H, vector<double> &Rx, double var, size_t SPA_I) {
	vector<size_t> permutation_seq(H.Col_number);
	//initial permutation_seq
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	MATRIX<__int8>  sorted_H, tmp_H;
	tmp_H.Building_Empty_Matrix(H.Row_number, H.Col_number);
	tmp_H = H;
	sorted_H = tmp_H;
	vector<double> orig_rx(H.Col_number);
	for (int i = 0; i < Rx.size(); i++)
		orig_rx.at(i) = Rx.at(i) * (2 / var);

	ABP_Permutation(Rx, tmp_H, sorted_H, permutation_seq);
	//
	// Change into Channel LLR value
	for (int i = 0; i < Rx.size(); i++)
		Rx.at(i) *= (2 / var);
	vector < vector <double> > M, Eji, tmp_M;
	vector <double> Li(tmp_H.Col_number), zi;
	int I_count = 0, Col = tmp_H.Col_number, Row = tmp_H.Row_number;
	double m_temp, temp;
	vector<__int8> codeword_seq(H.Col_number);

	vector<double> vect(Rx.size(), 0);   // for template container
	for (int i = 0; i < Col; i++) {
		vect.at(i) = (Rx.at(i));
	}
	for (int j = 0; j < Row; j++) {
		M.push_back(vect);
	}
	Li = vect;
	vect.clear();

	//

	while (true)
	{
		for (int j = 0; j < Row; j++) {

			for (int i = 0; i < Col; i++) {

				if (sorted_H._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < Col; k++) {

						if (sorted_H._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}

					temp = log((1 + ta) / (1 - ta));
					if (abs(temp) > LLR_MAX) {
						vect.push_back(temp > 0 ? LLR_MAX : (-1)*LLR_MAX);
					}
					else vect.push_back(temp);
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}

		for (int i = 0; i < Col; i++) {
			Li[i] = orig_rx.at(permutation_seq.at(i));
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] += Eji[j][i];
			}
		}

		//
		if (I_count >= SPA_I - 1) {
			break;
		}
		// End
		I_count++;
		for (int i = 0; i < Row; i++) {
			for (int j = Row; j < Col; j++) {
				if (sorted_H._matrix[i][j] == 1) {
					M[i][j] = Li.at(j) - Eji[i][j];
				}
			}
		}
		vect.clear();
		Eji.clear();
		zi.clear();
	}
	Desort_Function(permutation_seq, Li, Rx);
}

void ABP_MRIP_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	vector<double> tmp_Rx = decoding_info.rx_signal_seq;
	ABP(G,decoding_info);
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
		Best_Goal_2(message_length), // 用來接收下一層 A* 回傳的結果
		Child_Node(message_length);

	queue<NODE_PATH> Stack;
	Stack.push(Pointer);

	Best_Goal.metric = FLT_MAX;
	do {
		vector<double> sorting_rx_signal_seq(decoding_info.rx_signal_seq);
		Determine_Permutation(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index);
		Sort_Function(tmp_Rx, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
	} while (0);
	decoding_info.rx_signal_seq = tmp_Rx;
	//test
	static long err_count[128];
	static int count = 0;
	int err_bit = 0;
	count++;
	vector<__int8> Sorted_codeword(codeword_length, 0);
	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	for (int i = 0; i < message_length; i++) {
		if (Sorted_codeword[i] != Hard_RX[i]) {
			err_bit++;
		}
	}
	err_count[err_bit]++;
	if (count >= 5000000) {
		ClearFile("Error_Count.txt");
		for (int i = 0; i < message_length; i++) {
			WriteFile("Error_Count.txt", (double)err_count[i] / (double)count);
		}
	}
	decoding_info.rx_signal_seq = tmp_Rx;
	return;
	//test
	/*
	static long err_loc[128];
	static int count = 0;
	count++;
	vector<__int8> Sorted_codeword(codeword_length, 0);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	Sort_Function(decoding_info.estimated_codeword, Location_Index, codeword_seq);
	for (int i = 0; i < message_length; i++) {
		if (Sorted_codeword[i] != codeword_seq[i]) {
			err_loc[i]++;
		}
	}
	if (count >= 5000000) {
		ClearFile("Error_Loc.txt");
		for (int i = 0; i < message_length; i++) {
			WriteFile("Error_Loc.txt", (double)err_loc[i] / (double)count);
		}
	}
	return;
	*/
	decoding_info.Constraint_i = 4;
	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	do {
		Pointer = Stack.front();
		Stack.pop();
		if (Pointer.metric >= Best_Goal.metric) continue;
		//flipping node and push them into queue
		if (Pointer.D_z < decoding_info.Constraint_i && Pointer.level < message_length) {
			for (int i = Pointer.level; i < message_length; i++) {
				NODE_PATH tmp = Pointer;
				tmp.message_bits.at(i) ^= 1;
				++tmp.D_z;
				tmp.Diff_Index.push_back(i);
				tmp.level = i + 1;
				tmp.metric += Metric_Table._matrix[tmp.message_bits.at(i)][i];
				if (Pointer.metric < Best_Goal.metric)
					Stack.push(tmp);
			}
			decoding_info.STE += 2 * (message_length - Pointer.level);
		}
		//Update best goal
		codeword_seq = MRIP_codeword;
		for (int i = Pointer.level; i < message_length; i++) {
			Pointer.metric += Metric_Table._matrix[codeword_seq.at(i)][i];
		}
		if (Pointer.metric >= Best_Goal.metric) continue;
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
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric)
			Best_Goal = Pointer;
	} while (!Stack.empty());

	//
	Systematic_Linear_Block_Code_Encoder(Sorted_G, Best_Goal.message_bits, codeword_seq);
	Desort_Function(Location_Index, codeword_seq, decoding_info.estimated_codeword);


	//
	decoding_info.STE = decoding_info.STE / (double)message_length;
	decoding_info.COM = decoding_info.COM / (double)message_length;
	decoding_info.Binary_STE = decoding_info.Binary_STE / (double)message_length;

	//
	if (decoding_info.STE > decoding_info.Worst_Case_STE)
		decoding_info.Worst_Case_STE = decoding_info.STE;

	if (decoding_info.COM > decoding_info.Worst_Case_COM)
		decoding_info.Worst_Case_COM = decoding_info.COM;
}