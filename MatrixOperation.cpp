#include "PinAstar/MatrixOperation.h"


void GJ_Elimination(MATRIX<__int8> &input_matrix)
{
	size_t 
		temp(0),
		current_position_row(0),
		current_position_col(0);

	//笆篈秖痻皚
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
		
		//ㄏノヘ玡ㄤ
		if (current_position_col <  input_matrix.Col_number){
			for (size_t j(0); j < input_matrix.Row_number; ++j){
				//ㄏノ材 i  j  , i ぃ单 j		
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


	//传
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
	//程沧ヘ夹: 砞璸 G = [ P : I ]
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

inline void Matrix_Mul(MATRIX<__int8> & M1, MATRIX<__int8> & M2, MATRIX<__int8>& res) {
	if (M1.Col_number != M2.Row_number || !(M1.Row_number == res.Row_number)
		|| !(M2.Col_number == res.Col_number)) {
		cout << "Matrix_Mul error";
		return ;
	}
	int tmp;
	for (int i = 0; i < M1.Row_number; i++) {
		for (int j = 0; j < M2.Col_number; j++) {
			tmp = 0;
			for (int k = 0; k < M1.Col_number; k++) {
				tmp ^= M1._matrix[i][k] * M2._matrix[k][j];
			}
			res._matrix[i][j] = tmp;
		}
	}
}

void rref(MATRIX<__int8> &H, MATRIX<__int8> &out_G) {				//G柑
	//郎define
	int C = H.Col_number;
	int R = H.Row_number;
	int** C_swap = new int*[2];						//狦Τcolumn swap 魁
	C_swap[0] = new int[C];
	C_swap[1] = new int[C];
	int** G = new int*[R];
	for (int i = 0; i < R; i++) {
		G[i] = new int[C];
		for (int j = 0; j < C; j++) {
			G[i][j] = H._matrix[i][j];
		}
	}
	int C_index = 0;




	int x = 0, y = 0;
	while (x < R && y < C) {			//ref场だ
		int i = x;
		while (G[i][y] == 0) {			//т獶箂じ传玡
			i = i + 1;
			if (i == R) {				//тぃ1籔传︽
				C_index++;
				for (int h = R; h < C; h++) {		//眖痻皚秨﹍传
					if (G[x][h] == 1) {
						C_swap[0][C_index - 1] = y;
						C_swap[1][C_index - 1] = h;
						for (int k = 0; k < R; k++) {			//ㄢ︽ユ传
							int temp = G[k][y];
							G[k][y] = G[k][h];
							G[k][h] = temp;
						}
						break;
					}
				}
				i = x;
				break;
			}
		}
		if (i == R) y++;				//莱赣ぃ穦祇ネ
		else if (i != x) {				//籔Τ1传︽
			for (int j = y; j < C; j++) {
				int temp = G[i][j];
				G[i][j] = G[x][j];
				G[x][j] = temp;
			}
			for (int k = x + 1; k < R; k++) {
				if (G[k][y] == 1) {
					for (int j = y; j < C; j++) {
						G[k][j] ^= G[x][j];
					}
				}
			}

		}
		else {
			for (int k = x + 1; k < R; k++) {		//┕猭
				if (G[k][y] == 1) {
					for (int j = y; j < C; j++) {
						G[k][j] ^= G[x][j];
					}
				}
			}
		}
		y++;
		x++;
	}

	for (int k = 1, h = 1; k < x; k++, h++) {			//rref场だ

		for (int i = k - 1; i >= 0; i--) {
			if (G[i][h] == 1) {
				for (int j = h; j < C; j++) {
					G[i][j] ^= G[k][j];
				}
			}
		}
	}
	out_G.Building_Empty_Matrix(C - R, C);
	//output
	for (int i = 0; i < R; i++) {
		for (int j = 0; j < C; j++) {
			out_G._matrix[i][j] = G[i][j];
		}
	}
}
