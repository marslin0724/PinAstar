#include "PinAstar/MatrixOperation.h"


void GJ_Elimination(MATRIX<__int8> &input_matrix)
{
	size_t 
		temp(0),
		current_position_row(0),
		current_position_col(0);

	//�ŧi�ʺA�V�q�B�x�}
	/***************************************************************/
	MATRIX<__int8> Register(input_matrix);
	vector<size_t>
		Row_Index(input_matrix.Row_number, 0),
		Row_Sorter(input_matrix.Row_number, 0);
	for (size_t i(0); i < input_matrix.Row_number; ++i) 
		Row_Sorter.at(i) = i;
	/***************************************************************/
	

	while (current_position_row < input_matrix.Row_number){
		current_position_col = 0;
		/***************************************************************/
		do{
			if (Register._matrix[current_position_row][current_position_col] != 0) break;
			else ++current_position_col;
		} while ((current_position_col <  input_matrix.Col_number));
		/***************************************************************/
		Row_Index.at(current_position_row) = current_position_col + 1;
		
		//�ϥΥثe�C�h���h��L�C
		if (current_position_col <  input_matrix.Col_number){
			for (size_t j(0); j < input_matrix.Row_number; ++j){
				//�ϥβ� i �C�h�[ j �C , i ������ j		
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


	//���C
	/*********************************************************************/
	for (size_t i(0); i < input_matrix.Row_number; ++i){
		for (size_t j(i); j < input_matrix.Row_number; ++j){
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
}//end GJ_Elimination

void Sort_Matrix_Col(
	MATRIX<__int8>		&input_matrix,
	vector<size_t>		&permutation_seq, 
	MATRIX<__int8>		&output_matrix)
{
	for (size_t i(0); i < input_matrix.Col_number; ++i){
		for (size_t j(0); j < input_matrix.Row_number; ++j)
			output_matrix._matrix[j][i] = input_matrix._matrix[j][permutation_seq.at(i)];
	}
} //end Sorting_Matrix()

void Sort_Matrix_Row(
	MATRIX<__int16>		&input_matrix,
	vector<size_t>		&permutation_seq,
	MATRIX<__int16>		&output_matrix)
{
	for (size_t i(0); i < input_matrix.Row_number; ++i) {
		for (size_t j(0); j < input_matrix.Col_number; ++j) {
			output_matrix._matrix[i][j] = input_matrix._matrix[permutation_seq.at(i)][j];
		}
	}
}

void H_G_convertor_G_P_I(MATRIX<__int8> & H, vector<int> &permutation_seq, MATRIX<__int8> & G) {
	//�̲ץؼ�: �]�p G = [ P : I ]
	// Permutation matrix set up
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	// Use Col Operation to find independent columns

	MATRIX<__int8> Duplicated_H = H;
	vector<size_t> Searching_Order(H.Col_number);
	size_t j(0), temp(0);
		
	for (size_t i(H.Col_number - H.Row_number); i < H.Col_number; ++i) {
		Searching_Order.at(j) = i;
		++j;
	}
	for (size_t i(0); i < H.Col_number - H.Row_number; ++i) {
		Searching_Order.at(j) = i;
		j++;
	}
		
	for (size_t i(0); i < Duplicated_H.Row_number; ++i) {
		if (Duplicated_H._matrix[i][Searching_Order.at(i)] == 0) {
			size_t searcher = i;
			while (Duplicated_H._matrix[i][Searching_Order.at(++searcher)] == 0) {
			}
			Col_Exchange(Duplicated_H, (int)Searching_Order.at(i), (int)Searching_Order.at(searcher));
			permutation_seq.at(Searching_Order.at(i)) ^= permutation_seq.at(Searching_Order.at(searcher));
			permutation_seq.at(Searching_Order.at(searcher)) ^= permutation_seq.at(Searching_Order.at(i));
			permutation_seq.at(Searching_Order.at(i)) ^= permutation_seq.at(Searching_Order.at(searcher));
		}
		for (size_t j(i + 1); j < Duplicated_H.Col_number; ++j) {
			if (Duplicated_H._matrix[i][Searching_Order.at(j)] == 1) {
				Col_Addition(Duplicated_H, Searching_Order.at(i), Searching_Order.at(j));
			}
		}
	}
	/*
	for (size_t i(0); i < Duplicated_H.Row_number; ++i) {
		if (Duplicated_H._matrix[i][i] == 0) {
			size_t searcher = i;
			while (Duplicated_H._matrix[i][++searcher] == 0) {
			}
			Col_Exchange(Duplicated_H, i, searcher);
			permutation_seq.at(i) ^= permutation_seq.at(searcher);
			permutation_seq.at(searcher) ^= permutation_seq.at(i);
			permutation_seq.at(i) ^= permutation_seq.at(searcher);
		}
		for (size_t j(i + 1); j < Duplicated_H.Col_number; ++j) {
			if (Duplicated_H._matrix[i][j] == 1) {
				Col_Addition(Duplicated_H, i, j);
			}
		}
	}
	*/
	Duplicated_H = H;
	Sort_Matrix_Col_Forward(Duplicated_H, permutation_seq);

	
	
	for (size_t i(0); i < H.Row_number - 1; ++i) {
		for (size_t j(i + 1); j < H.Row_number; ++j) {
			if (Duplicated_H._matrix[j][Searching_Order.at(i)] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
	/*
	for (size_t i(0); i < H.Row_number - 1; ++i) {
		for (size_t j(i + 1); j < H.Row_number; ++j) {
			if (Duplicated_H._matrix[j][i] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
	*/
	
	// Jordan elimination
	for (size_t i(1); i < Duplicated_H.Row_number; ++i) {
		//cout << "i = " << i << " ";
		for (int j(i - 1); j >= 0; --j) {
			//cout << "j = " << j << " ";
			if (Duplicated_H._matrix[j][Searching_Order.at(i)] == 1) {
				//cout << "! ";
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
	/*
	for (size_t i(1); i < Duplicated_H.Row_number; ++i) {
		//cout << "i = " << i << " ";
		for (int j(i - 1); j >= 0; --j) {
			//cout << "j = " << j << " ";
			if (Duplicated_H._matrix[j][i] == 1) {
				//cout << "! ";
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
	*/
	/*
	cout << "Step4 : \n";
	for (size_t i(0); i < Duplicated_H.Row_number; ++i) {
		for (size_t j(0); j < Duplicated_H.Col_number; ++j) {
			cout << (int)Duplicated_H._matrix[i][j];
		}
		cout << "\n";
	}
	system("PAUSE");
	*/
	/*
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < H.Col_number / 2; ++j) {
			cout << (int)Duplicated_H._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < H.Col_number / 2; ++j) {
			cout << (int)Duplicated_H._matrix[i][j + (H.Col_number / 2)];
		}
		cout << endl;
	}

	system("pause");*/

	// Build G: [ I | P ] transfer into [ PT| I ]

	G.Building_Empty_Matrix(Duplicated_H.Col_number - Duplicated_H.Row_number, Duplicated_H.Col_number);
	size_t parity_num = Duplicated_H.Col_number - Duplicated_H.Row_number;//k
	//cout << parity_num << endl;
	for (size_t i(0); i < parity_num; ++i) {
		for (size_t j(0); j < parity_num; ++j) {
			if (j == i) G._matrix[i][j] = 1;
			else G._matrix[i][j] = 0;
		}
	}

	for (size_t i(0); i < parity_num; ++i) {
		for (size_t j(parity_num); j < Duplicated_H.Col_number; ++j) {
			G._matrix[i][j] = Duplicated_H._matrix[j - parity_num][i];
		}
	}
	
	/*
	// Build G: [ I | P ] transfer into [ PT| I ]

	G.Building_Empty_Matrix(Duplicated_H.Col_number - Duplicated_H.Row_number, Duplicated_H.Col_number);
	size_t parity_num = Duplicated_H.Col_number - Duplicated_H.Row_number;
	//cout << parity_num << endl;
	for (size_t i(0); i < parity_num; ++i) {
		for (size_t j(0); j < Duplicated_H.Row_number; ++j) {
			G._matrix[i][j] = Duplicated_H._matrix[j][i + Duplicated_H.Row_number];
		}
	}

	for (size_t i(0); i < parity_num; ++i) {
		for (size_t j(Duplicated_H.Row_number); j < Duplicated_H.Col_number; ++j) {
			if (j - Duplicated_H.Row_number == i) G._matrix[i][j] = 1;
			else G._matrix[i][j] = 0;
		}
	}
	*/

	H = Duplicated_H;

	/*
	cout << endl << endl;
	for (int i = 0; i < G.Row_number; ++i) {
		for (int j = 0; j < 128; ++j) {
			cout << (int)G._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < G.Row_number; ++i) {
		for (int j = 0; j < 192; ++j) {
			cout << (int)G._matrix[i][j + 128];
		}
		cout << endl;
	}

	system("pause");*/

	/*
	G_H_Multiple_Test(G, H);
	system("pause");
	
	cout << endl << "New H:" << endl;
	Duplicated_H.Show_matrix(1);
	system("pause");
	*/
}
void Col_Exchange(MATRIX<__int8> &Matrix, int col1, int col2) {
	for (int i = 0; i < Matrix.Row_number; ++i) {
		Matrix._matrix[i][col1] ^= Matrix._matrix[i][col2];
		Matrix._matrix[i][col2] ^= Matrix._matrix[i][col1];
		Matrix._matrix[i][col1] ^= Matrix._matrix[i][col2];
	}
}

void Row_Exchange(MATRIX<__int8> &Matrix, int row1, int row2) {
	for (int i = 0; i < Matrix.Col_number; ++i) {
		Matrix._matrix[row1][i] ^= Matrix._matrix[row2][i];
		Matrix._matrix[row2][i] ^= Matrix._matrix[row2][i];
		Matrix._matrix[row2][i] ^= Matrix._matrix[row2][i];
	}
}

void Row_Addition(MATRIX<__int8> &Matrix, int row1, int row2) { // value of col1 add on col2 
	//if (row1 == 126 && row2 == 128) {
		//cout << (int)Matrix._matrix[row1][126] << (int)Matrix._matrix[row1][127] << (int)Matrix._matrix[row1][128] << endl;
	//}
	for (int i = 0; i < Matrix.Col_number; ++i) {
		Matrix._matrix[row2][i] ^= Matrix._matrix[row1][i];
	}
}

void Col_Addition(MATRIX<__int8> &Matrix, int col1, int col2) // value of col1 add on col2 
{
	for (int i = 0; i < Matrix.Row_number; ++i) {
		Matrix._matrix[i][col2] ^= Matrix._matrix[i][col1];
	}
}

// Forward  :  H  -> H'
// Backward :  H' -> H

void Sort_Matrix_Col_Forward(MATRIX<__int8> &Matrix, vector<int> &permutation_seq) {
	MATRIX<__int8> M = Matrix;
	//vector<size_t> seq = permutation_seq;
	for (size_t traversal = 0; traversal < Matrix.Col_number; ++traversal) {
		for (size_t row = 0; row < Matrix.Row_number; ++row) {
			Matrix._matrix[row][traversal] = M._matrix[row][permutation_seq.at(traversal)];
		}

		/*
		for (size_t col = 0; col < Matrix.Col_number; ++col) {
			for (size_t row = 0; row < Matrix.Row_number; ++row) {
				Matrix._matrix[row][col] = M._matrix[row][permutation_seq.at(col)];
			}
			//size_t temp = seq.at(col);
			//if (traversal == temp) { // Do exchange
				//if (traversal < 10) cout << traversal << "/" << temp << "/" << col << endl;


				//seq.at(traversal) ^= seq.at(col);
				//seq.at(col) ^= seq.at(traversal);
				//seq.at(traversal) ^= seq.at(col);
				//break;
			//}
		}

		*/
	}
}
void G_H_Multiple_Test(MATRIX<__int8> G, MATRIX<__int8> H) {
	cout << "%%% Examine G and H %%%" << endl;
	cout << "H: row = " << H.Row_number << " , col = " << H.Col_number << endl;
	cout << "G: row = " << G.Row_number << " , col = " << G.Col_number << endl << endl;

	for (size_t i(0); i < G.Row_number; ++i) {
		for (size_t j(0); j < H.Row_number; ++j) {
			size_t temp = 0;
			for (size_t k(0); k < H.Col_number; ++k) {
				temp += (H._matrix[j][k] & G._matrix[i][k]);
			}
			temp %= 2;
			if (temp == 1) cout << "Error: H " << j << "th, G " << i << "th,temp=" << temp << endl;
		}
	}
}

void Transpose_Matrix(MATRIX<__int8> &H, MATRIX<__int8> &Transpose_H) {
	for (size_t i(0); i < H.Row_number; i++) {
		for (size_t j(0); j < H.Col_number; j++) {
			Transpose_H._matrix[j][i] = (int)H._matrix[i][j];
		}
	}
}