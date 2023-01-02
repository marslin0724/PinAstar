#include "PinAstar/LDPC_OPER.h"

/*

Functions for the LDPC & A* combination project
Including all LDPC functions such as SPA,sequence operation and matrix operation
Note that this part is Inner Codes connecting with the channel

*/


void ReadFile_H_Matrix_ver2(string name, MATRIX<__int8> &H) { // with puncture version
	// File structure:
	// 1st line: Original Row
	// 2nd line: Original Col
	// 3rd line: Puncture Number
	// After 4th line:   (Row,Col)        1
	// After Puncturing: Row = Row - Puncture Number; Col = Col - Puncture Number; 
	fstream fp;
	fp.open(name, ios::out | ios::in);
	if (!fp) cout << "Fail to open file: " << endl;
	char ch;
	size_t i(0), j(0), Row_number(0), Col_number(0), Puncture_number(0);
	fp >> Row_number;
	fp >> Col_number;
	fp >> Puncture_number;
	cout << "Puncture_number: " << Puncture_number << endl;
	H.Building_Empty_Matrix(Row_number, Col_number);
	string line;
	vector<int> VectorTemp;
	getline(fp, line);
	while (getline(fp, line)) {
		int temp = 0;
		for (int i = 0; i != line.size(); ++i) {
			if (isdigit(line[i])) temp = temp * 10 + (line[i] - 48);
			else
			{
				if (temp != 0) VectorTemp.push_back(temp);
				temp = 0;
			}
			if (i == (line.size() - 1)) {
				VectorTemp.push_back(temp);
			}
		}
		H._matrix[(VectorTemp.at(0)) - 1][(VectorTemp.at(1)) - 1] = 1;
		VectorTemp.clear();
	}
	/*
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < H.Col_number / 2; ++j) {
			cout << (int)H._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < H.Col_number / 2; ++j) {
			cout << (int)H._matrix[i][j + (H.Col_number / 2)];
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
	for (int i = 0; i < H.Row_number; ++i) {
		for (int j = 0; j < 64; ++j) {
			cout << (int)H._matrix[i][j];
		}
		cout << endl;
	}

	system("pause");*/
	/*
	if (Puncture_number != 0) {
		H.Row_number -= Puncture_number;
		H.Col_number -= Puncture_number;
		for (int i = 0; i < Puncture_number; ++i) { // Delete Row by Puncturing
			H._matrix.erase(H._matrix.begin());
		}
		for (int i = 0; i < H.Row_number; ++i) {    // Delete Col by Puncturing
			H._matrix.at(i).erase(H._matrix.at(i).begin(), H._matrix.at(i).begin() + (Puncture_number));
			//H._matrix.at(i).erase(H._matrix.at(i).begin() + H.Col_number - H.Row_number - (Puncture_number / 2), H._matrix.at(i).begin() + H.Col_number - H.Row_number);
		}
	}*/
	fp.close();
}

LDPC_FUNC::LDPC_FUNC() {
	
}

void LDPC_FUNC::Col_Exchange(MATRIX<__int8> &Matrix, int col1, int col2) {
	for (int i = 0; i < Matrix.Row_number; ++i) {
		Matrix._matrix[i][col1] ^= Matrix._matrix[i][col2];
		Matrix._matrix[i][col2] ^= Matrix._matrix[i][col1];
		Matrix._matrix[i][col1] ^= Matrix._matrix[i][col2];
	}
}

void LDPC_FUNC::Row_Exchange(MATRIX<__int8> &Matrix, int row1, int row2) {
	for (int i = 0; i < Matrix.Col_number; ++i) {
		Matrix._matrix[row1][i] ^= Matrix._matrix[row2][i];
		Matrix._matrix[row2][i] ^= Matrix._matrix[row2][i];
		Matrix._matrix[row2][i] ^= Matrix._matrix[row2][i];
	}
}

void LDPC_FUNC::Row_Addition(MATRIX<__int8> &Matrix, int row1, int row2) { // value of col1 add on col2 
	//if (row1 == 126 && row2 == 128) {
		//cout << (int)Matrix._matrix[row1][126] << (int)Matrix._matrix[row1][127] << (int)Matrix._matrix[row1][128] << endl;
	//}
	for (int i = 0; i < Matrix.Col_number; ++i) {
		Matrix._matrix[row2][i] ^= Matrix._matrix[row1][i];
	}
}

void LDPC_FUNC::Col_Addition(MATRIX<__int8> &Matrix, int col1, int col2) // value of col1 add on col2 
{
	for (int i = 0; i < Matrix.Row_number; ++i) {
		Matrix._matrix[i][col2] ^= Matrix._matrix[i][col1];
	}
}

// Forward  :  H  -> H'
// Backward :  H' -> H

void LDPC_FUNC::Sort_Matrix_Col_Forward(MATRIX<__int8> &Matrix, vector<int> &permutation_seq){
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

void LDPC_FUNC::Sort_Sequence_Backward(vector<double> &Rx, vector<int> &permutation_seq){
	vector<int> seq = permutation_seq;
	for (size_t col = 0; col < Rx.size(); ++col) {
		for (size_t traversal= 0; traversal < Rx.size(); ++traversal) {
			if (col == seq.at(traversal)) {
				double temp = Rx.at(col);
				Rx.at(col) = Rx.at(traversal);
				Rx.at(traversal) = temp;
				seq.at(col) ^= seq.at(traversal);
				seq.at(traversal) ^= seq.at(col);
				seq.at(col) ^= seq.at(traversal);
				break;
			}
		}
	}
}

void LDPC_FUNC::Sort_Sequence_Backward(vector<__int8> &Rx, vector<int> &permutation_seq) {
	vector<int> seq = permutation_seq;
	for (size_t col = 0; col < Rx.size(); ++col) {
		for (size_t traversal = 0; traversal < Rx.size(); ++traversal) {
			if (col == seq.at(traversal)) {
				double temp = Rx.at(col);
				Rx.at(col) = Rx.at(traversal);
				Rx.at(traversal) = temp;
				seq.at(col) ^= seq.at(traversal);
				seq.at(traversal) ^= seq.at(col);
				seq.at(col) ^= seq.at(traversal);
				break;
			}
		}
	}
}

void LDPC_FUNC::Sort_Sequence_Forward(vector<double> &Rx, vector<int> &permutation_seq) {
	//MATRIX<__int8> M = Matrix;
	vector<double> rtemp = Rx;
	//vector<size_t> seq = permutation_seq;
	for (size_t traversal = 0; traversal < Rx.size(); ++traversal) {
		Rx.at(traversal) = rtemp.at(permutation_seq.at(traversal));
	}
}

void LDPC_FUNC::Sort_Sequence_Forward(vector<__int8> &Rx, vector<int> &permutation_seq) {
	//MATRIX<__int8> M = Matrix;
	vector<__int8> rtemp = Rx;
	//vector<size_t> seq = permutation_seq;
	for (size_t traversal = 0; traversal < Rx.size(); ++traversal) {
		Rx.at(traversal) = rtemp.at(permutation_seq.at(traversal));
	}
}


void LDPC_FUNC::Show_Matrix(MATRIX<__int8> &Matrix,int row_min ,int row_max, int col_min, int col_max)
{
	std::cout << std::endl << "Show Matrix:" << std::endl;
	for (size_t i(row_min); i < row_max; ++i) {
		for (size_t j(col_min); j < col_max; ++j) {
			std::cout << (int)Matrix._matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

vector<double> LDPC_FUNC::Seq_Add(vector<double> cc, vector<double> nn) {
	int num;
	vector<double> result;

	if (cc.size() > nn.size()) num = nn.size();
	else num = cc.size();

	for (int i = 0; i < num; i++) {
		result.push_back(cc[i] + nn[i]);
		//cout << cc[i] + nn[i] << endl;
	}
	return result;
}

vector<size_t> LDPC_FUNC::Message_Generator(size_t length) {

	srand((unsigned long)time(NULL));  //產生亂數
	vector<size_t> bin;

	for (int i = 0; i < length; i++) {
		bin.push_back(rand() * 2 / RAND_MAX);
		//cout << rand() * 2 / RAND_MAX << endl;
	}
	return bin;
}

vector<size_t> LDPC_FUNC::Codeword_Generator(MATRIX<__int8> &Matrix, vector<size_t> &mes) {

	vector <size_t> codes(Matrix.Col_number, 0);
	for (int i = 0; i < Matrix.Col_number; i++) {
		int temp = 0;
		for (int j = 0; j < Matrix.Row_number; j++) {
			temp ^= (mes[j] & Matrix._matrix[j][i]);
		}
		codes.at(i) = temp;
		//codes.push_back(temp);
	}
	return codes;
}


vector<__int8> LDPC_FUNC::Codeword_Generator(MATRIX<__int8> &Matrix, vector<__int8> &mes) {

	vector <__int8> codes(Matrix.Col_number, 0);
	for (int i = 0; i < Matrix.Col_number; i++) {
		int temp = 0;
		for (int j = 0; j < Matrix.Row_number; j++) {
			temp ^= (mes[j] & Matrix._matrix[j][i]);
		}
		codes.at(i) = temp;
		//codes.push_back(temp);
	}
	return codes;
}


vector<size_t> LDPC_FUNC::Estimation(vector <double> LLR) {
	vector<size_t> Est(LLR.size(), 0);
	//cout << "AAA";
	for (size_t i(0); i < LLR.size(); ++i) {
		//cout << LLR.at(i) << " ";
		if (LLR.at(i) == 1) Est.at(i) = 1;
	}
	return Est;
}

int LDPC_FUNC::BER_Counting(vector<size_t> & CodeSeq, vector<size_t> & EstSeq, size_t parity) {
	int error = 0;
	for (size_t i(parity); i < CodeSeq.size(); ++i) {
		if (CodeSeq.at(i) != EstSeq.at(i)) ++error;
	}
	return error;
}

vector<double> LDPC_FUNC::PAM(vector<size_t> codew) {
	vector < double > Pam(codew.size(),0);
	for (int i = 0; i < codew.size(); i++) {
		if (codew.at(i) == 0) Pam.at(i) = 1;
		else  Pam.at(i) = -1;
	}
	return Pam;
}

double LDPC_FUNC::Variance(double SNR_value) {
	//cout << pow(10, SNR_value / 10);
	return (0.5 / pow(10, SNR_value / 10));
}

vector<double> LDPC_FUNC::Noise_Generator(int num, double mean, double var) {

	vector<double> noisee;
	//cout << var << endl;
	unsigned seeed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seeed);
	std::normal_distribution<double> distribution(mean, var);

	for (int i = 0; i < num; i++) {
		noisee.push_back(distribution(generator));
		//cout << noisee.at(i);
	}
	noisee.push_back(mean);
	noisee.push_back(var);

	return noisee;
}







void LDPC_FUNC::Print_Matrix_Col(MATRIX<__int8> Matrix, int col) {
	cout << "The " << col << "th col is:" << endl;
	for (int i = 0; i < Matrix.Row_number; ++i) {
		cout << (int)Matrix._matrix[i][col] << " ";
	}
	cout << endl;
}

void LDPC_FUNC::Print_Matrix_Row(MATRIX<__int8> Matrix, int row) {
	cout << "The " << row << "th row is:" << endl;
	for (int i = 0; i < Matrix.Col_number; ++i) {
		cout << (int)Matrix._matrix[row][i] << " ";
	}
	cout << endl;
}

void LDPC_FUNC::H_Test(MATRIX<__int8> H) {
	cout << "%%% Examine H %%%" << endl;
	for (size_t i(0); i < H.Row_number; ++i) {
		for (size_t j(0); j < H.Row_number; ++j) {
			if ((H._matrix[i][j] == 0) && (i == j)) cout << "?????????????H[" << i << "][" << j << "]=0" << endl;
			else if((H._matrix[i][j] == 1) && (i != j))  cout << "H[" << i << "][" << j << "]=1" << endl;
			else {

			}
		}
	}
}

void LDPC_FUNC::G_H_Multiple_Test(MATRIX<__int8> G, MATRIX<__int8> H) {
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
			if (temp == 1) cout << "Error: H " << j << "th, G " << i << "th,temp=" << temp <<  endl;
		}
	}
}



void LDPC_FUNC::ReadFile_H_Matrix(string name, MATRIX<__int8> &H) {

	fstream fp;
	fp.open(name, ios::out | ios::in);
	if (!fp) cout << "Fail to open file: " << endl;
	char ch;
	size_t i(0), j(0), Row_number(0), Col_number(0);
	fp >> Row_number;
	fp >> Col_number;
	H.Building_Empty_Matrix(Row_number, Col_number);

	fp.get(ch);
	while (fp.get(ch)) {
		if (ch == '\n') {
			++i;
			j = 0;
		}
		else if (ch == ' ') continue;
		else {
			H._matrix[i][j] = ch - '0';
			++j;
		}
	}
	fp.close();

}

void LDPC_FUNC::H_G_convertor_G_I_P(MATRIX<__int8> H, vector<int> &permutation_seq, MATRIX<__int8> & G) {
	//最終目標: 設計 G = [ I : P ]

	// Permutation matrix set up
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	int parity_num = H.Col_number - H.Row_number;
	int redundant_row = 0;

	// Use Col Operation to find independent columns
	MATRIX<__int8> Duplicated_H = H;
	for (int i(Duplicated_H.Row_number-1); i >-1; --i) {
		int I_position = i + parity_num;
		if (Duplicated_H._matrix[i][I_position] == 0) {
			int searcher = I_position;
			while (Duplicated_H._matrix[i][--searcher]==0) {

			}
			
			//cout <<"(" << searcher <<")";
			if (searcher < 0) {
				redundant_row = i + 1;
				//cout << endl << i << endl;
				//system("pause");
				break;
			}
			Col_Exchange(Duplicated_H, I_position, searcher);
			permutation_seq.at(I_position) ^= permutation_seq.at(searcher);
			permutation_seq.at(searcher) ^= permutation_seq.at(I_position);
			permutation_seq.at(I_position) ^= permutation_seq.at(searcher);
		}
		//cout << "A";
		// 確保為真的linear independent
		for (int j(I_position-1); j >= 0; --j) {
			if (Duplicated_H._matrix[i][j] == 1) {
				Col_Addition(Duplicated_H, I_position, j);
			}
		}
		//cout << i << endl;
	}
	//cout << endl;
	//system("pause");

	
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
	system("pause");
	

	Duplicated_H = H;
	Sort_Matrix_Col_Forward(Duplicated_H, permutation_seq);

	for (int j = 0; j < H.Row_number; ++j) {
		cout << (int)Duplicated_H._matrix[j][j+H.Col_number-H.Row_number];
	}
	system("pause");

	// Gaussian Elimination
	for (size_t i(Duplicated_H.Row_number - 1); i > 0; --i) {
		int I_position = i + parity_num;
		for (size_t j(i - 1); j > -1; --j) {
			if (Duplicated_H._matrix[j][I_position] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}

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

	system("pause");

	// Jordan elimination
	for (size_t i(Duplicated_H.Row_number - 2); i > 0; --i){
		//cout << "i = " << i << " ";
		int I_position = i + parity_num;
		for (int j(i+1); j < Duplicated_H.Row_number; ++j) {
			//cout << "j = " << j << " ";
			if (Duplicated_H._matrix[j][I_position] == 1) {
				//cout << "! ";
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
	
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

	system("pause");
	
	
	// Build G: [ I | P ] transfer into [ PT| I ]

	G.Building_Empty_Matrix(Duplicated_H.Col_number - Duplicated_H.Row_number, Duplicated_H.Col_number);
	//cout << parity_num << endl;
	for (size_t i(0); i < parity_num; ++i) {  // col
		for (size_t j(0); j < Duplicated_H.Row_number; ++j) {   // row
			G._matrix[i][j+ parity_num] = Duplicated_H._matrix[j][i];
		}
	}

	for (size_t i(0); i < parity_num; ++i) {
		G._matrix[i][i] = 1;
	}


	
	G_H_Multiple_Test(G, Duplicated_H);
	system("pause");
	cout << endl << "New H:" << endl;
	Duplicated_H.Show_matrix(1);
	system("pause");
	
}

void LDPC_FUNC::H_G_convertor_G_I_P_ver2(MATRIX<__int8> H, vector<int> &permutation_seq, MATRIX<__int8> & G) {
	//最終目標: 設計 G = [ I : P ]
	// Permutation matrix set up
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	int parity_num = H.Col_number - H.Row_number;
	//cout << parity_num << endl;
	int redundant_row = 0, lastrow;

	MATRIX<__int8> Duplicated_H = H;
	//cout << Duplicated_H.Row_number << "," << Duplicated_H.Col_number << endl;
	if (fullrank) lastrow = -1;
	else lastrow = 0;
	// Gaussian Elimination
	for (int i = H.Row_number - 1; i > lastrow; --i) {
		//cout << i << endl;
		int I_position = i + parity_num;
		int k = I_position;
		while (Duplicated_H._matrix[i][k] == 0) {
			if (k > parity_num) k = parity_num - 1;
			else --k;
		}
		if (k != I_position) {
			Col_Exchange(Duplicated_H, I_position, k);
			permutation_seq.at(I_position) ^= permutation_seq.at(k);
			permutation_seq.at(k) ^= permutation_seq.at(I_position);
			permutation_seq.at(I_position) ^= permutation_seq.at(k);
		}

		for (int j(i - 1); j > -1; --j) {
			if (Duplicated_H._matrix[j][I_position] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}
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

	// Jordan elimination
	for (int i(H.Row_number - 2); i > lastrow; --i) {
		int I_position = i + parity_num;
		for (int j(i + 1); j < H.Row_number; ++j) {
			if (Duplicated_H._matrix[j][I_position] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}

	// Build G: [ I | P ] transfer into [ PT| I ]

	G.Building_Empty_Matrix(Duplicated_H.Col_number - Duplicated_H.Row_number, Duplicated_H.Col_number);
	//cout << parity_num << endl;
	for (size_t i(0); i < parity_num; ++i) {  // col
		for (size_t j(0); j < Duplicated_H.Row_number; ++j) {   // row
			G._matrix[i][j + parity_num] = Duplicated_H._matrix[j][i];
		}
	}

	for (size_t i(0); i < parity_num; ++i) {
		G._matrix[i][i] = 1;
	}
	G_H_Multiple_Test(G, Duplicated_H);
	/*
	system("pause");
	cout << endl << "New H:" << endl;
	Duplicated_H.Show_matrix(1);
	system("pause");
	*/
}

void LDPC_FUNC::H_G_convertor_G_P_I(MATRIX<__int8> H, vector<int> &permutation_seq, MATRIX<__int8> & G) {
	//最終目標: 設計 G = [ P : I ]
	// Permutation matrix set up
	for (size_t i(0); i < permutation_seq.size(); ++i)
		permutation_seq.at(i) = i;

	// Use Col Operation to find independent columns

	MATRIX<__int8> Duplicated_H = H;

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

	Duplicated_H = H;
	Sort_Matrix_Col_Forward(Duplicated_H, permutation_seq);

	for (size_t i(0); i < H.Row_number - 1; ++i) {
		for (size_t j(i + 1); j < H.Row_number; ++j) {
			if (Duplicated_H._matrix[j][i] == 1) {
				Row_Addition(Duplicated_H, i, j);
			}
		}
	}

	// Jordan elimination
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
	G_H_Multiple_Test(G, Duplicated_H);
	system("pause");
	
	cout << endl << "New H:" << endl;
	Duplicated_H.Show_matrix(1);
	system("pause");
	*/
}




void SPA(MATRIX<__int8> &H, DECODING_INFO &decoding_info) {
	vector<double> R_temp;
	if (decoding_info.Return_Est_1_or_LLR_0 == 0) R_temp = decoding_info.rx_signal_seq;
	// Change into Channel LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) decoding_info.rx_signal_seq.at(i) *= (2 / decoding_info.var);
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp;
	long double temp;
	vector<__int8> Current_Result(Col, 0), New_Result(Col, 0);

	//cout << "Hi";
	//cout << Col << Row << endl;
	//system("pause");

	vector<double> vect(decoding_info.rx_signal_seq.size(), 0);   // for template container
    
	for (int i = 0; i < Col; i++) {
		vect.at(i) = (decoding_info.rx_signal_seq.at(i));
	}
	for (int j = 0; j < Row; j++) {
		M.push_back(vect);
	}
	vect.clear();

	while (TRUE)
	{

		for (int j = 0; j < Row; j++) {

			for (int i = 0; i < Col; i++) {

				if (H._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < Col; k++) {

						if (H._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}
					temp = log((1 + ta) / (1 - ta));
					if (abs(temp) > LLR_MAX) vect.push_back(temp > 0 ? LLR_MAX : (-1)*LLR_MAX);
					else vect.push_back(temp);
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}

		for (int i = 0; i < Col; i++) {
			Li.push_back(decoding_info.rx_signal_seq[i]);
			//Li.push_back(0);
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
			if (Li.at(i) > 0) New_Result.at(i) = 0;
			else New_Result.at(i) = 1;
		}

		// jump out of the loop (finished)
		
		if (New_Result != Current_Result) {
			if (Diff_Candidate_Test) {

			}
			//cout << I_count << endl;
			bool Breaker = TRUE;
			int temp;
			Current_Result = New_Result;
			for (int i = 0; i < Row; ++i) {
				temp = 0;
				for (int j = 0; j < Col; ++j) {
					temp ^= (Current_Result.at(j)&H._matrix[i][j]);
				}
				if (temp == 1) {
					Breaker = FALSE;
					break;
				}
			}
			if (LDPC_H_Check_Equation && Breaker == TRUE) {
				if (decoding_info.Return_Est_1_or_LLR_0 == 1) {  // Return Hard Rx
					for (int i = 0; i < Col; ++i) {
						if (Li.at(i) > 0) decoding_info.estimated_codeword.at(i) = 0;
						else decoding_info.estimated_codeword.at(i) = 1;
					}
				}
				else {                                           // Return LLR
					for (int i = 0; i < Col; ++i) {
						decoding_info.rx_signal_seq.at(i) = Li.at(i)*(decoding_info.var / 2);
					}
					//cout << "!";
				}
				//cout << "Success!" << endl;
				break;
			}
		}
		
		if (I_count > SPA_Iter) {
			if (decoding_info.Return_Est_1_or_LLR_0 == 1) {      // Return Hard Rx
				for (int i = 0; i < Col; ++i) {
					if (Li.at(i) > 0) decoding_info.estimated_codeword.at(i) = 0;
					else decoding_info.estimated_codeword.at(i) = 1;
				}
			}
			else {                                               // Return LLR
				decoding_info.rx_signal_seq = R_temp;
				
				for (int i = 0; i < Col; ++i) {
					decoding_info.rx_signal_seq.at(i) = Li.at(i)*(decoding_info.var / 2);
				}
				//cout << "!";
			}
			//cout << "Failed!" << endl;
			break;
		}
		// End

		I_count++;
		M.clear();
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				m_temp = decoding_info.rx_signal_seq[j];
				if (H._matrix[i][j] == 1) {
					for (int k = 0; k < Row; k++) {
						if (k != i && H._matrix[k][j] == 1) {
							m_temp = m_temp + Eji[k][j];
						}
					}
					vect.push_back(m_temp);
				}
				else vect.push_back(NULL);

			}
			M.push_back(vect);
			vect.clear();
		}
		Eji.clear();
		Li.clear();
		zi.clear();
	}
	//cout << "HHH";
	// 2018/11/12
}

void SPA_Partial(MATRIX<__int8> &H, DECODING_INFO &decoding_info) {
	vector<double> R_temp;
	LDPC_FUNC LDPC;
	R_temp = decoding_info.rx_signal_seq;
	decoding_info.code_seq_permutation.resize(decoding_info.rx_signal_seq.size());
	decoding_info.code_seq_MRIP.resize(decoding_info.rx_signal_seq.size());

	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		decoding_info.code_seq_permutation.at(i) = i;
		decoding_info.code_seq_MRIP.at(i) = 1;
	}
	for (int k = decoding_info.rx_signal_seq.size() - 1; k > 0; k--)
	{
		int temp;
		double ttemp;
		for (int j = 0; j < k; j++)
		{
			if (abs(R_temp.at(j))> abs(R_temp.at(j + 1))) {
				ttemp = R_temp.at(j + 1);
			    R_temp.at(j + 1) = R_temp.at(j);
				R_temp.at(j) = ttemp;
				temp = decoding_info.code_seq_permutation.at(j);
				decoding_info.code_seq_permutation.at(j) = decoding_info.code_seq_permutation.at(j + 1);
				decoding_info.code_seq_permutation.at(j + 1) = temp;
			}
		}
	}
	reverse(R_temp.begin(), R_temp.end());
	reverse(decoding_info.code_seq_permutation.begin(), decoding_info.code_seq_permutation.end());
	/*
	for (int i = 0; i < decoding_info.code_seq_MRIP.size(); ++i) {
		cout << (int)decoding_info.code_seq_permutation.at(i) << " ";
	}
	cout << endl << endl;
	*/
	//cout << decoding_info.rx_signal_seq.size()*(0.75) << endl;
	int start = decoding_info.rx_signal_seq.size()*(0.5);
	for (int i = start; i < decoding_info.rx_signal_seq.size(); ++i) {
		R_temp.at(i) /= 2;
		decoding_info.code_seq_MRIP.at(i) = 0;
	}

	LDPC.Sort_Sequence_Backward(R_temp, decoding_info.code_seq_permutation);
	LDPC.Sort_Sequence_Backward(decoding_info.code_seq_MRIP, decoding_info.code_seq_permutation);
	/*
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		cout << "(" << R_temp.at(i) << "," << (int)decoding_info.code_seq_MRIP.at(i) << ") ";
	}
	cout << endl << endl;
	system("pause"); 
	/*
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		cout << (int)decoding_info.code_seq_permutation.at(i) << " ";
	}
	cout << endl;
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		cout << (int)decoding_info.code_seq_MRIP.at(i) << " ";
	}

	cout << endl << endl;
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		cout << R_temp.at(i) << " ";
	}
	cout << endl << endl;
	system("pause");*/
	
	// Change into Channel LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) R_temp.at(i) *= (2 / decoding_info.var);
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp;
	vector<__int8> Current_Result(Col, 0), New_Result(Col, 0);

	//cout << Col << Row << endl;
	//system("pause");

	vector<double> vect(decoding_info.rx_signal_seq.size(), 0);   // for template container

	for (int i = 0; i < Col; i++) {
		vect.at(i) = (R_temp.at(i));
	}
	for (int j = 0; j < Row; j++) {
		M.push_back(vect);
	}
	vect.clear();

	while (TRUE)
	{

		for (int j = 0; j < Row; j++) {

			for (int i = 0; i < Col; i++) {

				if (H._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < Col; k++) {

						if (H._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}
					vect.push_back(log((1 + ta) / (1 - ta)));
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}

		for (int i = 0; i < Col; i++) {
			Li.push_back(R_temp.at(i));
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
			if (Li.at(i) > 0) New_Result.at(i) = 0;
			else New_Result.at(i) = 1;
		}

		// jump out of the loop (finished)
		/*
		if (New_Result != Current_Result) {
			//cout << I_count << endl;
			bool Breaker = TRUE;
			int temp;
			Current_Result = New_Result;
			for (int i = 0; i < Row; ++i) {
				temp = 0;
				for (int j = 0; j < Col; ++j) {
					temp ^= (Current_Result.at(j)&H._matrix[i][j]);
				}
				if (temp == 1) {
					Breaker = FALSE;
					break;
				}
			}
			if (Breaker == TRUE) {
				if (decoding_info.Return_Est_1_or_LLR_0 == 1) {  // Return Hard Rx
					for (int i = 0; i < Col; ++i) {
						if (Li.at(i) > 0) decoding_info.estimated_codeword.at(i) = 0;
						else decoding_info.estimated_codeword.at(i) = 1;
					}
				}
				else {                                           // Return LLR
					//for (int i = 0; i < Col; ++i) {
						//decoding_info.rx_signal_seq.at(i) = Li.at(i)*(decoding_info.var / 2);
					//}
					//cout << "!";
				}
				//cout << "Success!" << endl;
				break;
			}
		}
		*/
		if (I_count > SPA_Iter) {
			
			for (int i = 0; i < Col; ++i) {
				Li.at(i) *= (decoding_info.var / 2);
			}
			/*
			for (int i = 0; i < Col; ++i) {
				//cout << 
				if (decoding_info.code_seq_MRIP.at(i) == 1) decoding_info.rx_signal_seq.at(i) = Li.at(i);
			}*/
			//int count = 0;
			for (int i = 0; i < Col; ++i) {
				if (decoding_info.code_seq_MRIP.at(i) == 1) {
					//++count;
					//abs(Li.at(i)) < abs(decoding_info.rx_signal_seq.at(i)) ||
					if (abs(Li.at(i)) < abs(decoding_info.rx_signal_seq.at(i)) || Li.at(i)*decoding_info.rx_signal_seq.at(i) < 0) decoding_info.code_seq_MRIP.at(i) = 2;
				}
			}
			R_temp = decoding_info.rx_signal_seq;
			for (int k = decoding_info.rx_signal_seq.size() - 1; k > 0; k--)
			{
				int temp;
				double ttemp;
				for (int j = 0; j < k; j++)
				{
					if (abs(R_temp.at(j)) > abs(R_temp.at(j + 1))) {
						ttemp = R_temp.at(j + 1);
						R_temp.at(j + 1) = R_temp.at(j);
						R_temp.at(j) = ttemp;
						temp = decoding_info.code_seq_MRIP.at(j);
						decoding_info.code_seq_MRIP.at(j) = decoding_info.code_seq_MRIP.at(j + 1);
						decoding_info.code_seq_MRIP.at(j + 1) = temp;
					}
				}
			}
			reverse(decoding_info.code_seq_MRIP.begin(), decoding_info.code_seq_MRIP.end());
			//cout << count << endl << endl;
			break;
		}
		// End

		I_count++;
		M.clear();
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				m_temp = R_temp.at(j);
				if (H._matrix[i][j] == 1) {
					for (int k = 0; k < Row; k++) {
						if (k != i && H._matrix[k][j] == 1) {
							m_temp = m_temp + Eji[k][j];
						}
					}
					vect.push_back(m_temp);
				}
				else vect.push_back(NULL);

			}
			M.push_back(vect);
			vect.clear();
		}
		Eji.clear();
		Li.clear();
		zi.clear();
	}
	//cout << "HHH";
	// 2018/11/12
}

void SPA(MATRIX<__int8> &H, vector<double> &Rx, double var) {

	// Change into Channel LLR value
	for (int i = 0; i < Rx.size(); i++) Rx.at(i) *= (2 / var);
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp;
	vector<__int8> Current_Result(Col, 0), New_Result(Col, 0);

	//cout << Col << Row << endl;
	//system("pause");

	vector<double> vect(Rx.size(), 0);   // for template container

	for (int i = 0; i < Col; i++) {
		vect.at(i) = (Rx.at(i));
	}
	for (int j = 0; j < Row; j++) {
		M.push_back(vect);
	}
	vect.clear();

	while (TRUE)
	{

		for (int j = 0; j < Row; j++) {

			for (int i = 0; i < Col; i++) {

				if (H._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < Col; k++) {

						if (H._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}
					vect.push_back(log((1 + ta) / (1 - ta)));
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}

		for (int i = 0; i < Col; i++) {
			Li.push_back(Rx[i]);
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
			if (Li.at(i) > 0) New_Result.at(i) = 0;
			else New_Result.at(i) = 1;
		}

		// jump out of the loop (finished)

		if (New_Result != Current_Result) {
			//cout << I_count << endl;
			bool Breaker = TRUE;
			int temp;
			Current_Result = New_Result;
			for (int i = 0; i < Row; ++i) {
				temp = 0;
				for (int j = 0; j < Col; ++j) {
					temp ^= (Current_Result.at(j)&H._matrix[i][j]);
				}
				if (temp == 1) {
					Breaker = FALSE;
					break;
				}
			}
			if (Breaker == TRUE) {
				for (int i = 0; i < Col; ++i) {
					Rx.at(i) = Li.at(i);
				}
				break;
			}
		}
		if (I_count > SPA_Iter) {
			for (int i = 0; i < Col; ++i) {
				Rx.at(i) = Li.at(i);
			}
			break;
		}
		// End

		I_count++;
		M.clear();
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				m_temp = Rx[j];
				if (H._matrix[i][j] == 1) {
					for (int k = 0; k < Row; k++) {
						if (k != i && H._matrix[k][j] == 1) {
							m_temp = m_temp + Eji[k][j];
						}
					}
					vect.push_back(m_temp);
				}
				else vect.push_back(NULL);

			}
			M.push_back(vect);
			vect.clear();
		}
		Eji.clear();
		Li.clear();
		zi.clear();
	}
}

void LDPC_VC_RBP(MATRIX<__int8> &H, DECODING_INFO &decoding_info) {
	
	int Col = H.Col_number,
		Row = H.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		Pass,
		binary_temp,
		Target_Row,
		Target_Col,
		Target_Index_Metric;
	double temp;

	vector<__int8> Hard_Result(Col, 0), Parity_Check_Stack(Row, 0);
	vector<double> Channel_LLR(Col, 0), Final_Channel_LLR(Col, 0);

	MATRIX<double>
		V2C_matrix,
		C2V_matrix,
		Residual_matrix,
		Zero_matrix;

	V2C_matrix.Building_Empty_Matrix(Row, Col);
	C2V_matrix.Building_Empty_Matrix(Row, Col);
	Residual_matrix.Building_Empty_Matrix(Row, Col);
	Zero_matrix.Building_Empty_Matrix(Row, Col);
	
	for (int index = 0; index < Col; ++index) {
		Channel_LLR.at(index) = 2 * decoding_info.rx_signal_seq.at(index) / decoding_info.var;
		//cout << Channel_LLR.at(index) << ", ";
	}
	for (int row = 0; row < Row; row++) {
		for (int col = 0; col < Col; col++) {
			if (H._matrix[row][col]) {
				//cout << "!";
				V2C_matrix._matrix[row][col] = Channel_LLR.at(col);
				Residual_matrix._matrix[row][col] = abs(Channel_LLR.at(col));
			}
		}
	}
	/*
	Target_Index_Metric = 0;
	for (int index = 0; index < Col; ++index) {
		if (abs(Channel_LLR.at(index)) > Target_Index_Metric) {
			Target_Index_Metric = abs(Channel_LLR.at(index));
			Target_Col = index;
		}
	}

	for (int vertical = 0; vertical < Row && H._matrix[vertical][Target_Col]; ++vertical) {
		Residual_matrix._matrix[vertical][Target_Col] = 0;
		// Find the C2V of the other variable nodes except the target in the same row
		for (int vnode = 0; vnode < Col && vnode != Target_Col && H._matrix[vertical][vnode]; vnode++) {
			temp = 1;
			for (int index = 0; index < Col && index!=vnode; index++) {
				temp*= tanh(V2C_matrix._matrix[vertical][index] / 2);
			}
			temp = log((1 + temp) / (1 - temp));
			if (abs(temp) < LLR_MAX) C2V_matrix._matrix[vertical][vnode] = temp;
			else C2V_matrix._matrix[vertical][vnode] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));

			for (int cnode = 0; cnode < Row && cnode != vertical && H._matrix[cnode][vnode]; cnode++) {
				temp = V2C_matrix._matrix[cnode][vnode];
				V2C_matrix._matrix[cnode][vnode] = Channel_LLR.at(vnode);
				for (int index = 0; index < Row && index != cnode; index++) {
					V2C_matrix._matrix[cnode][vnode] += C2V_matrix._matrix[index][vnode];
				}
				Residual_matrix._matrix[cnode][vnode] = abs(V2C_matrix._matrix[cnode][vnode] - temp);
			}
		}
	}*/
	/*
	for (int row = 0; row < Row; row++) {
		for (int col = 0; col < Col/4; col++) {
			cout << Residual_matrix._matrix[row][col] << ", ";
		}
		cout << endl;
	}*/

	while (1) {
		for (int row = 0; row < Row; row++) {
			binary_temp = 0;
			for (int col = 0; col < Col; col++) {
				if (H._matrix[row][col] && (V2C_matrix._matrix[row][col] < 0)) binary_temp ^= 1;
			}
			if (binary_temp == 1) {
				Parity_Check_Stack.at(row) = 1;
			}
			else Parity_Check_Stack.at(row) = 0;
			//cout << endl;
		}
		/*
		for (int row = 0; row < Row; row++) {
			for (int col = 0; col < Col; col++) {
				if(H._matrix[row][col])cout << Residual_matrix._matrix[row][col] << ", ";
			}
			cout << endl;
		}
		cout << endl;
		system("pause");*/

		Target_Index_Metric = 0;
		for (int row = 0; row < Row; row++) {
			//if (Parity_Check_Stack.at(row) == 1) continue;
			for (int col = 0; col < Col; col++) {
				if (H._matrix[row][col] && Residual_matrix._matrix[row][col] > Target_Index_Metric) {
					Target_Col = col;
					Target_Row = row;
					Target_Index_Metric = Residual_matrix._matrix[row][col];
				}
			}
		}
		//cout << "Target: " << Target_Row << "," << Target_Col << " with " << Residual_matrix._matrix[Target_Row][Target_Col] << endl << endl;

		/*
		for (int row = 0; row < Row; row++) {
			for (int col = 0; col < Col && H._matrix[row][col]; col++) {
				Residual_matrix._matrix[row][col] = 0;
			}
		}*/
		Residual_matrix._matrix[Target_Row][Target_Col] = 0;

		for (int vnode = 0; vnode < Col; vnode++) {
			if (vnode != Target_Col && H._matrix[Target_Row][vnode]) {
				temp = 1;
				for (int index = 0; index < Col; index++) {
					if(index != vnode && H._matrix[Target_Row][index]) temp *= tanh(V2C_matrix._matrix[Target_Row][index] / 2);
				}
				temp = log((1 + temp) / (1 - temp));
				if (abs(temp) < LLR_MAX) C2V_matrix._matrix[Target_Row][vnode] = temp;
				else C2V_matrix._matrix[Target_Row][vnode] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));

				for (int cnode = 0; cnode < Row; cnode++) {
					if (cnode != Target_Row && H._matrix[cnode][vnode]) {
						temp = V2C_matrix._matrix[cnode][vnode];
						V2C_matrix._matrix[cnode][vnode] = Channel_LLR.at(vnode);
						for (int index = 0; index < Row; index++) {
							if(index != cnode) V2C_matrix._matrix[cnode][vnode] += C2V_matrix._matrix[index][vnode];
						}
						Residual_matrix._matrix[cnode][vnode] = abs(V2C_matrix._matrix[cnode][vnode] - temp);
						//Residual_matrix._matrix[cnode][vnode] = abs(V2C_matrix._matrix[cnode][vnode]);
						//if (Residual_matrix._matrix[cnode][vnode] == 0) cout << "No";
					}
				}
			}
		}

		if (Iteration_Count++ > 8) break;

		Pass = 1;
		for (int col = 0; col < Col; col++) {
			temp = Channel_LLR.at(col);
			for (int row = 0; row < Row; row++) {
				if(H._matrix[row][col]) temp += C2V_matrix._matrix[row][col];
			}
			if (temp > 0) Hard_Result.at(col) = 0;
			else Hard_Result.at(col) = 1;
		}
		for (int row = 0; row < Row; row++) {
			binary_temp = 0;
			for (int col = 0; col < Col; col++) {
				binary_temp ^= (Hard_Result.at(col)&H._matrix[row][col]);
			}
			if (binary_temp) {
				Pass = 0;
				break;
			}
		}
		if (Pass) break;
	}
	//cout <<"Iteration:"<< Iteration_Count << endl;
	for (int col = 0; col < Col; col++) {
		//Channel_LLR.at(col) = 0;
		for (int row = 0; row < Row; row++) {
			if (H._matrix[row][col]) {
				Channel_LLR.at(col) += C2V_matrix._matrix[row][col];
				//cout << "(" << C2V_matrix._matrix[row][col] << ", " << V2C_matrix._matrix[row][col] << ")";
				//cout << "(" << C2V_matrix._matrix[row][col] << ", ";
			}
		}
		//cout << endl;
		decoding_info.rx_signal_seq.at(col) = Channel_LLR.at(col);
		if (Channel_LLR.at(col) > 0) decoding_info.estimated_codeword.at(col) = 0;
		else decoding_info.estimated_codeword.at(col) = 1;
	}
}

void LDPC_Modified_RBP(MATRIX<__int8> &H, DECODING_INFO &decoding_info) {

	int Col = H.Col_number,
		Row = H.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		Pass,
		binary_temp,
		Target_Row,
		Target_Col,
		Target_Index_Metric;
	double temp;

	vector<__int8> Hard_Result(Col, 0), Parity_Check_Stack(Row, 0);
	vector<double> Channel_LLR(Col, 0), Final_Channel_LLR(Col, 0);

	MATRIX<double>
		V2C_matrix,
		C2V_matrix,
		Residual_matrix,
		Zero_matrix;

	V2C_matrix.Building_Empty_Matrix(Row, Col);
	C2V_matrix.Building_Empty_Matrix(Row, Col);
	Residual_matrix.Building_Empty_Matrix(Row, Col);
	Zero_matrix.Building_Empty_Matrix(Row, Col);

	for (int index = 0; index < Col; ++index) {
		Channel_LLR.at(index) = 2 * decoding_info.rx_signal_seq.at(index) / decoding_info.var;
		//cout << Channel_LLR.at(index) << ", ";
	}
	for (int row = 0; row < Row; row++) {
		for (int col = 0; col < Col; col++) {
			if (H._matrix[row][col]) {
				V2C_matrix._matrix[row][col] = Channel_LLR.at(col);
				//Residual_matrix._matrix[row][col] = abs(Channel_LLR.at(col));
			}
		}
	}
	// Parity_Check_Stack


	while (1) {
		/*
		for (int row = 0; row < Row; row++) {
			for (int col = 0; col < Col; col++) {
				if(H._matrix[row][col])cout << Residual_matrix._matrix[row][col] << ", ";
			}
			cout << endl;
		}
		cout << endl;
		system("pause");*/
		//Parity_Check_Stack
		for (int row = 0; row < Row; row++) {
			binary_temp = 0;
			for (int col = 0; col < Col; col++) {
				if (H._matrix[row][col] && (V2C_matrix._matrix[row][col] < 0)) binary_temp ^= 1;
			}
			if (binary_temp == 1) {
				Parity_Check_Stack.at(row) = 1;
			}
			else Parity_Check_Stack.at(row) = 0;
			cout << endl;
		}

		Target_Index_Metric = 0;
		for (int row = 0; row < Row; row++) {
			//cout << Parity_Check_Stack.at(row) << endl;
			//if (Parity_Check_Stack.at(row) == 0) continue;
			for (int col = 0; col < Col; col++) {
				if (H._matrix[row][col] && Residual_matrix._matrix[row][col] > Target_Index_Metric) {
					Target_Col = col;
					Target_Row = row;
					Target_Index_Metric = Residual_matrix._matrix[row][col];
				}
			}
		}
		//cout << "Target: " << Target_Row << "," << Target_Col << " with " << Residual_matrix._matrix[Target_Row][Target_Col] << endl << endl;

		/*
		for (int row = 0; row < Row; row++) {
			for (int col = 0; col < Col && H._matrix[row][col]; col++) {
				Residual_matrix._matrix[row][col] = 0;
			}
		}*/
		Residual_matrix._matrix[Target_Row][Target_Col] = 0;

		for (int vnode = 0; vnode < Col; vnode++) {
			if (vnode != Target_Col && H._matrix[Target_Row][vnode]) {
				temp = 1;
				for (int index = 0; index < Col; index++) {
					if (index != vnode && H._matrix[Target_Row][index]) temp *= tanh(V2C_matrix._matrix[Target_Row][index] / 2);
				}
				temp = log((1 + temp) / (1 - temp));
				if (abs(temp) < LLR_MAX) C2V_matrix._matrix[Target_Row][vnode] = temp;
				else C2V_matrix._matrix[Target_Row][vnode] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));

				for (int cnode = 0; cnode < Row; cnode++) {
					if (cnode != Target_Row && H._matrix[cnode][vnode]) {
						temp = V2C_matrix._matrix[cnode][vnode];
						V2C_matrix._matrix[cnode][vnode] = Channel_LLR.at(vnode);
						for (int index = 0; index < Row; index++) {
							if (index != cnode) V2C_matrix._matrix[cnode][vnode] += C2V_matrix._matrix[index][vnode];
						}
						Residual_matrix._matrix[cnode][vnode] = abs(V2C_matrix._matrix[cnode][vnode] - temp);
						//if (Residual_matrix._matrix[cnode][vnode] == 0) cout << "No";
					}
				}
			}
		}

		if (Iteration_Count++ > 50) break;

		Pass = 1;
		for (int col = 0; col < Col; col++) {
			//temp = Channel_LLR.at(col);
			temp = 0;
			for (int row = 0; row < Row; row++) {
				if (H._matrix[row][col]) temp += C2V_matrix._matrix[row][col];
			}
			if (temp > 0) Hard_Result.at(col) = 0;
			else Hard_Result.at(col) = 1;
		}
		for (int row = 0; row < Row; row++) {
			binary_temp = 0;
			for (int col = 0; col < Col; col++) {
				binary_temp ^= (Hard_Result.at(col)&H._matrix[row][col]);
			}
			if (binary_temp) {
				Pass = 0;
				break;
			}
		}
		if (Pass) break;
	}
	//cout <<"Iteration:"<< Iteration_Count << endl;
	int Sign_Change;
	for (int col = 0; col < Col; col++) {
		Sign_Change = 0;
		temp = Channel_LLR.at(col);

		Channel_LLR.at(col) = 0;
		for (int row = 0; row < Row; row++) {
			if (H._matrix[row][col]) {
				if (temp*C2V_matrix._matrix[row][col] < 0) Sign_Change++;
				Channel_LLR.at(col) += C2V_matrix._matrix[row][col];
				//cout << "(" << C2V_matrix._matrix[row][col] << ", " << V2C_matrix._matrix[row][col] << ")";
				//cout << "(" << C2V_matrix._matrix[row][col] << ", ";
			}
		}
		//if (Sign_Change == 0) Channel_LLR.at(col) = (Channel_LLR.at(col) > 0 ? 10000 : -10000);
		//cout << endl;
		
		//decoding_info.rx_signal_seq.at(col) = Channel_LLR.at(col);
		decoding_info.rx_signal_seq.at(col) /= (Sign_Change + 1);
		
		if (Channel_LLR.at(col) > 0) decoding_info.estimated_codeword.at(col) = 0;
		else decoding_info.estimated_codeword.at(col) = 1;
	}
}

/*
vector<double> LDPC_FUNC::SPA(double var, vector <double> ri, MATRIX<__int8> h_matrix, int Iter) {

	bool ReturnLLR = true; // Determine whether the function return LLR value (else return binary output)


	//cout << varr << endl;
	//cout << endl << "   " << ri.size() <<endl;
	//double constantnum = 4 / (pow(10, var / 10));
	//cout << "const" << constantnum << endl;

	for (int i = 0; i < ri.size(); i++) ri.at(i) *= (-1)*(2 / var);
	//
	//cout << "Ri is:" << endl;
	//print_sequence_d(ri);

	//for (int i = 0; i < ri.size(); i++) ri.at(i) *= ri.at(i);

	vector < vector <double> > M;
	vector < vector <double> > Eji;
	vector <double> Li;
	vector <double> zi;
	int I_count = 0;
	double m_temp;

	//cout << "pow";
	vector<double> vect(ri.size(), 0);   // for template container
	//cout << vect.size();
	for (int i = 0; i < h_matrix.Col_number; i++) {
		//cout << "hey";
		vect.at(i) = (ri.at(i));
		//cout << i << " ";
	}
	for (int j = 0; j < h_matrix.Row_number; j++) {
		M.push_back(vect);
	}
	//cout << "pow";
	//cout << "M is:" <<endl;
	//print_matrix_d(M);
	vect.clear();
	//cout << "Testing:" << endl;
	//print_matrix(h_matrix);
	bool loop = true;
	while (loop == true)
	{

		for (int j = 0; j < h_matrix.Row_number; j++) {

			for (int i = 0; i < h_matrix.Col_number; i++) {

				if (h_matrix._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < h_matrix.Col_number; k++) {

						if (h_matrix._matrix[j][k] == 1 && k != i) {
							//ta = ta * (exp(M[j][k] / 2) - exp((-1) * M[j][k] / 2)) / (exp(M[j][k] / 2) + exp((-1) * M[j][k] / 2));
							ta *= tanh(M[j][k] / 2);

							//ta = ta * ((1 - exp((-1)*M[j][k])) / (1 + exp((-1)*M[j][k]))); // 2018/11/14

						}
					}
					//cout << ta << " " << log((1 + ta) / (1 - ta)) << endl;
					vect.push_back(log((1 + ta) / (1 - ta)));
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}
		//cout << "E:" << endl;
		//print_matrix_d(Eji);

		for (int i = 0; i < h_matrix.Col_number; i++) {
			Li.push_back(ri[i]);
			for (int j = 0; j < h_matrix.Row_number; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
		}

		/*
		cout << "Li\n";
		for (int j = 0; j <30; j++) {
			cout <<  setprecision(2) << Li[j] <<  " ";
		}
		cout << I_count << endl;
		system("pause");
		*/

		//print_sequence_d(Li);





/*
		for (int i = 0; i < Li.size(); i++) {
			if (Li[i] > 0) zi.push_back(0);
			else zi.push_back(1);
		}
		/*
		for (int i = 0; i < 30; i++) {
			cout << zi.at(i);
		}
		cout << endl;
		*/
		//print_sequence_d(zi);

		/*?????????????????????????????????????????????????????????????????
		bool breaker = true;
		int tempi;
		double tempid;
		vector<double> Hz;
		for (int i = 0; i < h_matrix.Row_number; i++) {
			tempi = 0;
			for (int j = 0; j < zi.size(); j++) {
				tempi = (int)(tempi + (h_matrix._matrix[i][j] * zi[j])) % 2;
			}
			tempid = (double)tempi;
			Hz.push_back(tempid);
		}
		*/
		/*  ??????????????????????????????
		int zero = 0;
		for (int i = 0; i < Hz.size(); i++) zero = zero + Hz[i];
		if (zero == 0) {
			//cout << "Decoding successly!" << endl;
			//cout << "Iteration time:" << I_count << endl;
			//cout << endl << "The result in alorithm:" << endl;
			//print_sequence_d(zi);
			//cout << endl;
			return zi;
		}
		*/
		// jump out of the loop (finished)



		/*
		if (I_count > Iter) {
			if (ReturnLLR == true) {
				for (int kk = 0; kk < Li.size(); kk++) {
					//Li.at(kk) = Li.at(kk)*(-1);
				}
				return Li;
			}
			//cout << "Out of limit of iteration!" << endl << "Iteration time: 10" << endl;
			else {
				for (int i = 0; i < Li.size(); i++) {
					if (Li[i] > 0) zi.push_back(0);
					else zi.push_back(1);
				}
				return zi;
			}
		}

		//cout << "The " << I_count << "th result is:" <<endl;
		//print_sequence_d(zi);
		I_count++;
		M.clear();
		for (int i = 0; i < h_matrix.Row_number; i++) {
			for (int j = 0; j < h_matrix.Col_number; j++) {
				m_temp = ri[j];
				if (h_matrix._matrix[i][j] == 1) {
					for (int k = 0; k < h_matrix.Row_number; k++) {
						if (k != i && h_matrix._matrix[k][j] == 1) {
							m_temp = m_temp + Eji[k][j];
						}
					}
					vect.push_back(m_temp);
				}
				else vect.push_back(NULL);

			}
			M.push_back(vect);
			vect.clear();
		}
		//cout << "M:" << endl;
		//print_matrix_d(M);
		Eji.clear();
		Li.clear();
		zi.clear();
		//Hz.clear();
	}
	// 2018/11/12
}
*/
/*
vector<double> LDPC_FUNC::SPA(double var, vector <double> ri, MATRIX<__int8> h_matrix) {

	for (int i = 0; i < ri.size(); i++) ri.at(i) *= (2/var);

	vector < vector <double> > M;
	vector < vector <double> > Eji;
	vector <double> Li;
	vector <double> zi;
	//zi.resize(ri.size());
	int I_count = 0;
	double m_temp;
	//cout << endl << h_matrix.Col_number << " " << h_matrix.Row_number << " " << ri.size() << endl;

	vector<double> vect(ri.size(),0);   // for template container

	for (int i = 0; i < h_matrix.Col_number; i++) {
		vect.at(i)=(ri.at(i));
	}
	for (int j = 0; j < h_matrix.Row_number; j++) {
		M.push_back(vect);
	}

	vect.clear();

	bool loop = true;
	while (loop == true)
	{

		for (int j = 0; j < h_matrix.Row_number; j++) {

			for (int i = 0; i < h_matrix.Col_number; i++) {

				if (h_matrix._matrix[j][i] == 1) {
					double ta = 1;

					for (int k = 0; k < h_matrix.Col_number; k++) {

						if (h_matrix._matrix[j][k] == 1 && k != i) {
							ta *= tanh(M[j][k] / 2);
						}
					}
					vect.push_back(log((1 + ta) / (1 - ta)));
				}
				else vect.push_back(NULL);
			}
			Eji.push_back(vect);
			vect.clear();
		}
		
		for (int i = 0; i < h_matrix.Col_number; i++) {
			Li.push_back(ri[i]);
			for (int j = 0; j < h_matrix.Row_number; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
		}

		for (int i = 0; i < Li.size(); i++) {
			if (Li[i] > 0) zi.push_back(0);
			else zi.push_back(1);
		}
		//print_sequence_d(zi);
		bool breaker = true;
		int tempi;
		double tempid;
		vector<double> Hz;
		//cout  <<"Zi:"<< zi.size() << endl;
		for (int i = 0; i < h_matrix.Row_number; i++) {
			tempi = 0;
			for (int j = 0; j < zi.size(); j++) {
				tempi = (int)(tempi + (h_matrix._matrix[i][j] * zi.at(j))) % 2;
			}
			//cout <<I_count<< ":" << i << "," << tempi << endl;
			tempid = (double)tempi;
			Hz.push_back(tempid);
		}
		for (int i = 0; i < 20; i++) {
			cout << setprecision(1) << zi.at(i) << " ";
		}
		cout << endl;
		
		if (I_count > 10) {
			return zi;
		}

		I_count++;
		M.clear();
		for (int i = 0; i < h_matrix.Row_number; i++) {
			for (int j = 0; j < h_matrix.Col_number; j++) {
				m_temp = ri[j];
				if (h_matrix._matrix[i][j] == 1) {
					for (int k = 0; k < h_matrix.Row_number; k++) {
						if (k != i && h_matrix._matrix[k][j] == 1) {
							m_temp = m_temp + Eji[k][j];
						}
					}
					vect.push_back(m_temp);
				}
				else vect.push_back(NULL);

			}
			M.push_back(vect);
			vect.clear();
		}

		Eji.clear();
		Li.clear();
		zi.clear();
		Hz.clear();
	}
	// 2018/11/12
}
*/

/*

Several of A* sequence
Used for saving sequence data seperately
Note that this part is outer codes in vicinity to source coding

*/

Astar_OuterCode::Astar_OuterCode(int messagelength, int codelength, int AmountOfCodes, int LDPCRow,int LDPCCol) {
	MessageLength = messagelength;
	CodeLength = codelength;
	CodeAmount = AmountOfCodes;
	LDPC_RowNumber = LDPCRow;
	LDPC_ColNumber = LDPCCol;
	// Message Bits
	for (int i = 0; i < AmountOfCodes; ++i) {
		vector<__int8> Message;
		for (int j = 0; j < messagelength; ++j) {
			Message.push_back(0);
		}
		Message_Seq.push_back(Message);
	}

	// CodeBits
	for (int i = 0; i < AmountOfCodes; ++i) {
		vector<__int8> Code;
		for (int j = 0; j < codelength; ++j) {
			Code.push_back(0);
		}
		CodeWord_Seq.push_back(Code);
	}

	// Error Bits
	for (int i = 0; i < AmountOfCodes; ++i) {
		vector<__int8> Error;
		for (int j = 0; j < codelength; ++j) {
			Error.push_back(0);
		}
		Error_Seq.push_back(Error);
	}

	// Received Bits
	for (int i = 0; i < AmountOfCodes; ++i) {
		vector<double> Rec;
		for (int j = 0; j < codelength; ++j) {
			Rec.push_back(0);
		}
		Received_Seq.push_back(Rec);
	}

	// Resize of Tx, Rx
	Tx_Signal_Seq_series.resize(LDPCCol,0);
	Rx_Signal_Seq_series.resize(LDPCCol,0);

}

void Astar_OuterCode::All_Generate_Message(void(*Func)(size_t , vector<__int8> &)) {
	for (int i = 0; i < CodeAmount; ++i) {
		(*Func)(MessageLength, Message_Seq[i]);
	}
}

void Astar_OuterCode::All_Systematic_Linear_Block_Code_Encoder(void(*Func)(MATRIX<__int8> &, vector<__int8> &, vector<__int8> &), MATRIX<__int8> &G) {
	for (int i = 0; i < CodeAmount; ++i) {
		(*Func)(G, Message_Seq[i], CodeWord_Seq[i]);
	}
}

void Astar_OuterCode::Seq_Combine() {
	for (int i = 0; i < CodeAmount; ++i) {
		LDPC_Message_Seq.insert(LDPC_Message_Seq.end(), CodeWord_Seq[i].begin(), CodeWord_Seq[i].end());
	}
}

void Astar_OuterCode::Seq_Seperate() {
	//int limit = LDPC_ColNumber;
	int id = -1;
	int counter = 0;
	for (int i = LDPC_RowNumber; i < LDPC_ColNumber; ++i) {
		if (i%CodeLength == 0) {
			++id;
			counter = 0;
		}
		Received_Seq[id][counter++] = Rx_Signal_Seq_series[i];
	}
}

void Astar_OuterCode::Astar_Initialization()
{

}

void  Astar_OuterCode::ParityCheck_Encoder() 
{
	int LastRow = CodeAmount - 1;
	for (int i = 0; i < CodeLength; ++i) {
		CodeWord_Seq[LastRow][i] = 0;
		for (int j = 0; j < LastRow; ++j) {
			CodeWord_Seq[LastRow][i]+= CodeWord_Seq[j][i];
		}
		//cout << i << ":" << (int)CodeWord_Seq[LastRow][i];
		CodeWord_Seq[LastRow][i] = CodeWord_Seq[LastRow][i] % 2;
		//cout << " " << (int)CodeWord_Seq[LastRow][i] << endl;
		/*
		for (int j = 0; j < LastRow+1; ++j) {
			cout << (int)CodeWord_Seq[j][i] << " ";
		}
		cout << endl;
		*/
	}
	
}

void Astar_OuterCode::ProductCode_Seq_Combine() {

	for (int i = 0; i < CodeLength; ++i) {
		for (int j = 0; j < CodeAmount; ++j) {
			// Inner_CodeWord_Seq[i*CodeAmount + j] = CodeWord_Seq[j][i];
			Inner_CodeWord_Seq.push_back(CodeWord_Seq[j][i]);
		}
	}
	// Combine方式:  (方便做Trellis)
	// 1 2 3                
	// 4 5 6    ->  1 4 7 2 5 8 3 6 9
	// 7 8 9
}

void Astar_OuterCode::ProductCode_Seq_Seperate() {

	for (int col = 0; col < CodeLength; col++) {
		for (int row = 0; row < CodeAmount; row++) {
			Received_Seq[row][col] = Rx_Signal_Seq_series.at(col*CodeAmount + row);
		}
	}
	// Seperate方式:  
	//                           1 2 3
	// 1 4 7 2 5 8 3 6 9   ->    4 5 6
	//                           7 8 9
}

void Astar_OuterCode::Trellis_Decoding() {
	//int counter = 0;
	int choice = 2;   //  1:Trellis 2:A*star
	for (int Module = 0; Module < CodeLength; Module++) {
		
		vector <int> Hard_Seq;
		// Find Hard Decision Result
		for (int node = 0; node < CodeAmount; ++node) {
			if (Rx_Signal_Seq_series.at(Module*CodeAmount + node) > 0) Hard_Seq.push_back(0);
			else Hard_Seq.push_back(1);
		}
		
		//Trellis Decoding
		vector <double> Rec_Seq;
		Rec_Seq.assign(Rx_Signal_Seq_series.begin() + (Module * CodeAmount), Rx_Signal_Seq_series.begin() + ((Module + 1)* CodeAmount));
		//cout << Hard_Seq.size() << " " << Rec_Seq.size() << endl;

		
		BinaryTree T1(Hard_Seq, Rec_Seq);
		if (choice == 1)	T1.BinaryTree_Construction(var);
		else if (choice == 2) 	T1.ConstructAstar();
		else ;
		//Return value to Recieved seqeunce

		// 將值傳回去的function
		
		for (int node = 0; node < CodeAmount; node++) {
			//if (Rx_Signal_Seq_series.at(node + (Module * CodeAmount)) != T1.ResultRecSeq.at(node)) {
				//++counter;
				//cout << Rx_Signal_Seq_series.at(node + Module * CodeAmount) << " " << T1.ResultRecSeq.at(node);
				Rx_Signal_Seq_series.at(node + (Module * CodeAmount)) = T1.ResultRecSeq.at(node);
			//}
		}
	}
	//cout << counter << endl;
}

void Astar_OuterCode::Trellis_Decoding(vector <int>& count) {
	//int counter = 0;
	int choice = 2;   //  1:Trellis 2:A*star
	for (int Module = 0; Module < CodeLength; Module++) {

		vector <int> Hard_Seq;
		// Find Hard Decision Result
		for (int node = 0; node < CodeAmount; ++node) {
			if (Rx_Signal_Seq_series.at(Module*CodeAmount + node) > 0) Hard_Seq.push_back(0);
			else Hard_Seq.push_back(1);
		}

		//Trellis Decoding
		vector <double> Rec_Seq;
		Rec_Seq.assign(Rx_Signal_Seq_series.begin() + (Module * CodeAmount), Rx_Signal_Seq_series.begin() + ((Module + 1)* CodeAmount));
		//cout << Hard_Seq.size() << " " << Rec_Seq.size() << endl;

		BinaryTree T1(Hard_Seq, Rec_Seq);
		if (choice == 1)	T1.BinaryTree_Construction(var);
		else if (choice == 2) 	T1.ConstructAstar();
		else;

		for (int run = 0; run < 8; ++run) {
			if (T1.ResultRecSeq.at(run) > 0) T1.ResultRecSeq.at(run) = 0;
			else  T1.ResultRecSeq.at(run) = 1;
			if (T1.ResultRecSeq.at(run) != Hard_Seq.at(run)) {
				count.at(run)++;
				cout << "test";
			}
			//DiffCounter
		}
		//Return value to Recieved seqeunce

		// 將值傳回去的function

		/*
		for (int node = 0; node < CodeAmount; node++) {
			//if (Rx_Signal_Seq_series.at(node + (Module * CodeAmount)) != T1.ResultRecSeq.at(node)) {
				//++counter;
				//cout << Rx_Signal_Seq_series.at(node + Module * CodeAmount) << " " << T1.ResultRecSeq.at(node);
				Rx_Signal_Seq_series.at(node + (Module * CodeAmount)) = T1.ResultRecSeq.at(node);
			//}
		}
		*/
	}
	//cout << counter << endl;
}

// 3. Reverse Trellis


void  Astar_OuterCode::ParityCheck_Encoder_Rev()
{
	int LastRow = CodeAmount - 1;
	for (int i = 0; i < CodeLength; ++i) {
		CodeWord_Seq[LastRow][i] = 0;
		for (int j = 0; j < LastRow; ++j) {
			CodeWord_Seq[LastRow][i] += CodeWord_Seq[j][i];
		}
		CodeWord_Seq[LastRow][i] = CodeWord_Seq[LastRow][i] % 2;
		Message_Seq[LastRow][i] = CodeWord_Seq[LastRow][i];
	}
}


void Astar_OuterCode::ProductCode_Seq_Combine_Rev() {

	for (int i = 0; i < CodeAmount; ++i) {
		for (int j = 0; j < CodeLength; ++j) {
			// Inner_CodeWord_Seq[i*CodeAmount + j] = CodeWord_Seq[j][i];
			Inner_CodeWord_Seq.push_back(CodeWord_Seq[i][j]);
		}
	}
	// Combine方式:  (方便做A* Decoding)
	// 1 2 3                
	// 4 5 6    ->  1 2 3 4 5 6 7 8 9
	// 7 8 9
}

void Astar_OuterCode::ProductCode_Seq_Seperate_Rev() {

	for (int col = 0; col < CodeAmount; col++) {
		for (int row = 0; row < CodeLength; row++) {
			Received_Seq[col][row] = Rx_Signal_Seq_series.at(col*CodeLength + row);
			//cout << "A";
		}
		//cout << endl;
	}
	// Seperate方式:  
	//                          1 2 3
	// 1 2 3 4 5 6 7 8 9  ->    4 5 6
	//                          7 8 9
}



// Binary-Tree for codeword candidates
BinaryTree::BinaryTree(vector<int> HardSeq, vector<double> RecSeq) {
	HardSequence.assign(HardSeq.begin(), HardSeq.end());
	RecSequence.assign(RecSeq.begin(), RecSeq.end());
}

void BinaryTree::BinaryTree_Construction(double var) {

	// Construct Metric Table
	for (int i = 0; i < HardSequence.size(); ++i) {
		if (HardSequence.at(i) == 1) {
			MetricTable_0.push_back((RecSequence.at(i) > 0) ? RecSequence.at(i) : RecSequence.at(i)*(-1));
			MetricTable_1.push_back(0);
		}
		else {
			MetricTable_0.push_back(0);
			MetricTable_1.push_back((RecSequence.at(i) > 0) ? RecSequence.at(i) : RecSequence.at(i)*(-1));
		}
	}
	// Construct Tree
	root = new TreeNode();
	ConstructTree(root, 0);
	//cout << NodeCount;
	//system("pause");

	// Bit-Flip
	ResultHardSeq.resize(HardSequence.size());
	ResultRecSeq.resize(HardSequence.size());
	TreeNode* TraversalPtr = BestNodePtr;
	double threshold = 2;
	double scalar = 1;
	int FlipCounter = 0;

	for (int j = (HardSequence.size() - 1); j > -1; --j) {
		ResultHardSeq.at(j) = TraversalPtr->BinaryBit;
		if (ResultHardSeq.at(j) != HardSequence.at(j)) {
			//cout << RecSequence.at(j) << " ";
			//ResultRecSeq.at(j) = RecSequence.at(j)*(-1);
			//ResultRecSeq.at(j) = 0;
			if (ResultHardSeq.at(j) == 0) ResultRecSeq.at(j) = (RecSequence.at(j)+ threshold)*scalar*(2/var);
			else ResultRecSeq.at(j) = (RecSequence.at(j) - threshold)*scalar*(2 / var);

			++FlipCounter;
		}
		else {
			ResultRecSeq.at(j) = RecSequence.at(j)*(2 / var);
		}
		TraversalPtr = TraversalPtr->parent;
	}
	/*
	for (int i = 0; i < ResultRecSeq.size(); ++i) {
		cout << ResultRecSeq.at(i) << " ";
	}
	system("pause");
	*/
	/*
	if (FlipCounter != 0) {
		for (int i = 0; i < HardSequence.size() - 1; ++i) {
			ResultRecSeq.at(i) *= scalar;
		}
	}
	*/
	/*
	for (int i = 0; i < HardSequence.size() - 1; ++i) {

	}
	*/

	DeleteTree(root);
}

void BinaryTree::ConstructTree(TreeNode* ptr, int depth) {
	//++NodeCount;
	if (depth == HardSequence.size() - 1) {   // (leaf-1)
		ptr->BinaryBitTotal %= 2;
		++depth;
		if (ptr->BinaryBitTotal == 0) {
			TreeNode* child_0 = new TreeNode(0, depth);
			ptr->upchild = child_0;
			child_0->parent = ptr;
			child_0->MetricTotal = ptr->MetricTotal + MetricTable_0.at(depth - 1);
			ConstructTree(child_0, depth);
		}
		else {
			TreeNode* child_1 = new TreeNode(1, depth);
			ptr->downchild = child_1;
			child_1->parent = ptr;
			child_1->MetricTotal = ptr->MetricTotal + MetricTable_1.at(depth - 1);
			ConstructTree(child_1, depth);
		}
	}

	else if (depth == HardSequence.size()) {   // leaf
		if (ptr->MetricTotal < BestMetric) {
			BestMetric = ptr->MetricTotal;
			BestNodePtr = ptr;
		}
	}

	else {                                    // not leaf nor (leaf-1)
		++depth;
		TreeNode* child_0 = new TreeNode(0, depth);
		ptr->upchild = child_0;
		child_0->parent = ptr;
		child_0->BinaryBitTotal = ptr->BinaryBitTotal + 0;
		child_0->MetricTotal = ptr->MetricTotal + MetricTable_0.at(depth - 1);
		ConstructTree(child_0, depth);

		TreeNode* child_1 = new TreeNode(1, depth);
		ptr->downchild = child_1;
		child_1->parent = ptr;
		child_1->BinaryBitTotal = ptr->BinaryBitTotal + 1;
		child_1->MetricTotal = ptr->MetricTotal + MetricTable_1.at(depth - 1);
		ConstructTree(child_1, depth);
	}
}

void BinaryTree::DeleteTree(TreeNode* ptr) {
	if (ptr->upchild != NULL) DeleteTree(ptr->upchild);
	if (ptr->downchild != NULL) DeleteTree(ptr->downchild);
	delete ptr;
}

void BinaryTree::ConstructAstar() {
	int AddAll = 0;
	for (int i = 0; i < HardSequence.size(); ++i) {
		AddAll = HardSequence.at(i);
	}
	if ((AddAll % 2) != 0) {
		//cout << "A" << endl;
		vector <double> SoftSaver;
		int MinPosition;
		double Reliability = DBL_MAX;
		for (int traversal = 0; traversal < HardSequence.size(); ++traversal) {
			if (RecSequence.at(traversal) < 0) SoftSaver.push_back(RecSequence.at(traversal)*(-1));
			else SoftSaver.push_back(RecSequence.at(traversal));
			if (SoftSaver.at(traversal) < Reliability) {
				MinPosition = traversal;
				Reliability = SoftSaver.at(traversal);
			}

			//cout << setprecision(1) << RecSequence.at(traversal) << " ";
		}
		//cout << endl << RecSequence.at(MinPosition) << endl << endl;
		RecSequence.at(MinPosition) *= (-1);
	}
	ResultRecSeq.assign(RecSequence.begin(), RecSequence.end());
}

int Cycle4_Check(MATRIX<__int8> &H) {
	int Result = TRUE, row_temp, col_temp, H_row = H.Row_number, H_col = H.Col_number;

	for (int row = 0; row < H_row; row++) {
		cout << "Complete row: " << row << endl;
		for (int col = 0; col < H_col; col++) {
			if (H._matrix[row][col] == 1) { // 4 cycle test start 
				for (int col_check = col + 1; col_check < H_col; col_check++) {
					if (H._matrix[row][col_check] == 1) {
						for (int row_check = row + 1; row_check < H_row; row_check++) {
							if (H._matrix[row_check][col_check] == 1) {
								for (int col_check_final = col_check - 1; col_check_final > -1; col_check_final--) {
									if (H._matrix[row_check][col_check_final] == 1 && col_check_final == col) {
										cout << "4cycle : " << row << "," << col << " / " << row << "," << col_check << " / "
											<< row_check << "," << col_check << " / " << row_check << "," << col_check_final << endl;
										Result = FALSE;
									}
								}
							}
						}
					}

				}
				
			}

		}
	}
	if (Result == TRUE) cout << "No four cycle!" << endl;
	else cout << "Exists four cycle!" << endl;
	system("pause");
	return Result;
}

int Cycle8_Check(MATRIX<__int8> &H) { // 未完成
	int Result = TRUE, row_temp, col_temp, H_row = H.Row_number, H_col = H.Col_number;

	for (int row = 0; row < H_row; row++) {
		cout << "Complete row: " << row << endl;
		for (int col = 0; col < H_col; col++) {
			if (H._matrix[row][col] == 1) { // 4 cycle test start 
				for (int col_check = col + 1; col_check < H_col; col_check++) {
					if (H._matrix[row][col_check] == 1) {
						for (int row_check = row + 1; row_check < H_row; row_check++) {
							if (H._matrix[row_check][col_check] == 1) {
								for (int col_check_final = col_check - 1; col_check_final > -1; col_check_final--) {
									if (H._matrix[row_check][col_check_final] == 1) {
										if (col_check_final == col) {
											cout << "4cycle : " << row << "," << col << " / " << row << "," << col_check << " / "
												<< row_check << "," << col_check << " / " << row_check << "," << col_check_final << endl;
											Result = FALSE;
										}
										else {
											
										}
									}
								}
							}
						}
					}

				}

			}

		}
	}
	if (Result == TRUE) cout << "No four cycle!" << endl;
	else cout << "Exists four cycle!" << endl;
	system("pause");
	return Result;
}