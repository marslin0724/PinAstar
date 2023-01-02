#include "AStarDecode.h"

void SPA_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H,G_;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	
	int i, j,start;
	vector<double> 
		subRx(H.Col_number,0),	//給每個sub block用 
		Rx(G.Col_number*(H.Col_number - H.Row_number)/H.Col_number),
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
		//SPA iteration
		SPA_subDecoder2(H, subRx, decoding_info.var,decoding_info.SPA_I);
		//ASPA_subdecoder(H, subRx, decoding_info.var, decoding_info.SPA_I);
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
	/*
	for (i = 0; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq.at(i) /= (2 / decoding_info.var);
	}*/
	//original A*
	A_star_I(G_,decoding_info);
	//set parity LLR estimate
	decoding_info.estimated_codeword.resize(G.Col_number);
	
	for (i = G_.Col_number; i < G.Col_number; i++) {
		if (Rx_parity.at(i - G_.Col_number) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = tmp_Rx;
}

void SPA_subDecoder(MATRIX<__int8> &H, vector<double> &Rx, double var,size_t SPA_I) {

	// Change into Channel LLR value
	for (int i = 0; i < Rx.size(); i++) Rx.at(i) *= (2 / var);
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp,temp;
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
		if (I_count > SPA_I) {
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

//不丟掉block 的parity
void SPA_A_star_ver2(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H, G_;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();

	int i, j, start;
	vector<double>
		subRx(H.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq;

	for (i = 0; i < G.Col_number / H.Col_number; i++) {
		//assign sub-Block
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * H.Row_number;
		for (; j < H.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j- H.Col_number + H.Row_number);
		}
		//SPA iteration
		SPA_subDecoder2(H, subRx, decoding_info.var,decoding_info.SPA_I);
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * ( H.Row_number);
		for (; j < H.Col_number; j++) {
			Rx.at(G_.Col_number + start + j-H.Col_number + H.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;
	
	A_star_I(G, decoding_info);
	decoding_info.rx_signal_seq = tmp_Rx;
}

void SPA_subDecoder2(MATRIX<__int8> &H, vector<double> &Rx, double var, size_t SPA_I) {
	//做小修改
	/*
	H._matrix.resize(H.Row_number+2);
	H._matrix.at(4).resize(H.Col_number);
	H._matrix.at(5).resize(H.Col_number);
	for (int i = 0; i < H.Col_number; i++) {
		H._matrix.at(4).at(i) =  H._matrix.at(0).at(i) ^ H._matrix.at(1).at(i);
		H._matrix.at(5).at(i) =  H._matrix.at(2).at(i) ^ H._matrix.at(3).at(i);
	}*/
	
	//
	// Change into Channel LLR value
	for (int i = 0; i < Rx.size(); i++) Rx.at(i) *= (2 / var);
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number,Row = H.Row_number;
	double m_temp,temp;
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
					
					temp = log((1 + ta) / (1 - ta));
					if (abs(temp) > LLR_MAX+10) {
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
			Li.push_back(Rx[i]);
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
			if (Li.at(i) > 0) New_Result.at(i) = 0;
			else New_Result.at(i) = 1;
		}

		

		
		if (I_count > SPA_I) {
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

void SPA_LLR(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	static double Avg_LLR[100];
	static int count;
	double temp_LLR = 0;
	count++;
	ClearFile("LLR.txt");
	vector<double> R_temp;
	if (decoding_info.Return_Est_1_or_LLR_0 == 0) R_temp = decoding_info.rx_signal_seq;
	// Change into Channel LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) decoding_info.rx_signal_seq.at(i) *= (2 / decoding_info.var);
	//建立H矩陣
	MATRIX<__int8> H;
	//做小修改 +2
	H.Building_Empty_Matrix(G.Col_number - G.Row_number, G.Col_number);
	H._matrix = G._matrix_outer;
	//做小修改
	/*H._matrix.resize(6);
	H._matrix.at(4).resize(8);
	H._matrix.at(5).resize(8);
	
	for (int i = 0; i < G.Col_number; i++) {
		H._matrix.at(4).at(i) = H._matrix.at(0).at(i) ^ H._matrix.at(1).at(i);
		H._matrix.at(5).at(i) = H._matrix.at(2).at(i) ^ H._matrix.at(3).at(i);
	}*/
	/*
	H._matrix = {
	{1,1,1,0,1,0,0,0},
	{1,0,0,1,1,0,0,1},
	{0,1,0,1,1,0,1,0},
	{1,0,1,0,0,1,0,1},
	{0,1,1,0,0,1,1,0},
	{0,0,0,1,0,1,1,1}
	};*/
	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp;
	long double temp;
	vector<__int8> Current_Result(Col, 0), New_Result(Col, 0);
	//計算初始LLR
	/*
	temp_LLR = 0;
	for (int i = 0; i < G.Row_number; i++) {
		temp_LLR += abs(decoding_info.rx_signal_seq.at(i));
	}
	Avg_LLR[0] += temp_LLR / G.Row_number;
	*/
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
			Li.push_back(decoding_info.rx_signal_seq[i]);
			//Li.push_back(0);
			for (int j = 0; j < Row; j++) {
				if (Eji[j][i] != NULL) Li[i] = Li[i] + Eji[j][i];
			}
			if (Li.at(i) > 0) New_Result.at(i) = 0;
			else New_Result.at(i) = 1;
		}

		// jump out of the loop (finished)
		/*
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
		*/
		if (I_count >= 10/*SPA_Iter*/) {
			if (/*decoding_info.Return_Est_1_or_LLR_0 == */1) {      // Return Hard Rx
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
		//計算平均LLR
		temp_LLR = 0;
			for (int i = 0; i < G.Row_number; i++) {
				temp_LLR += abs(Li.at(i));
			}
			Avg_LLR[I_count] += temp_LLR / G.Row_number;

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
	
	/*
	WriteFile("LLR.txt", "Avg_LLR");
	WriteFile("LLR.txt", "-----------");
	for (int i = 0; i < 51; i++) {
		WriteFile("LLR.txt", Avg_LLR[i]/count);
	}
	WriteFile("Hard_test.txt", "-----------");
	*/
}

void SPA_LLR2(MATRIX<__int8> &G, DECODING_INFO &decoding_info/*,double* message, Matrix m*/) {
	//bool* res = (bool*)malloc(sizeof(bool) * C);
	//printf("r:\n");
	//test
	ClearFile("LLR.txt");
	static int LLR_dis[1001];
	static double Var, Mean;
	static int count;
	count++;
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.code_seq.at(i) == 1) {//1 mapping to -1
			decoding_info.rx_signal_seq.at(i) += 2;
		}
	}
	//
	vector<vector<__int8>> H = G._matrix_outer;
	int C = H.at(0).size();
	int R = H.size();
	double* r = new double[C];
	for (int i = 0; i < C; i++) {								//設定機率矩陣r
		r[i] = decoding_info.rx_signal_seq.at(i) * (2 / decoding_info.var);
		//printf(" %.1f ,", r[i]);
	}
	//printf("\n\n");

	double **M = new double*[R];
	double **E = new double*[R];
	double *L = new double[C];
	for (int i = 0; i < R; i++) {
		M[i] = new double[C];
		E[i] = new double[C];
		for (int j = 0; j < C; j++) {
			E[i][j] = 0;
		}
	}


	int I_max = 0;

	do {


		for (int i = 0; i < R; i++) {	//設立M矩
			int k = 0;
			for (int j = 0; j < C; j++) {
				if (H[i][j] == 1) {
					M[i][j] = r[j];
					//printf(" %.1f ,", M[i][j]);
				}
				else {
					M[i][j] = 0;
					//printf(" %.1f ,", M[i][j]);
				}

			}
			//printf("\n");
		}
		//printf("\n\n");
		//printf("M :\n");
		for (int j = 0; j < C; j++) {			//對M矩陣做迭代
			for (int i = 0; i < R; i++) {
				if (H[i][j] == 1) {
					for (int k = 0; k < R; k++) {
						if (k != i) {
							M[i][j] += E[k][j];
						}
					}
				}
			}
		}

		for (int i = 0; i < R; i++) {
			for (int j = 0; j < C; j++) {
				//printf(" %.1f ,", M[i][j]);
			}
			//printf("\n");
		}
		//printf("\n\n");



		//printf("E :\n");
		for (int i = 0; i < R; i++) {				//設立E矩陣
			for (int j = 0; j < C; j++) {
				if (H[i][j] == 1) {
					double temp = 1;
					for (int k = 0; k < C; k++) {
						if (H[i][k] == 1 && k != j) {
							temp *= tanh(M[i][k] / 2);
						}
					}
					if (fabs(1 - temp) <= pow(10, -307)) {	//防止爆掉
						E[i][j] = 706;
						//printf(" %.1f ,", E[i][j]);
					}
					else if (fabs(1 + temp) <= pow(10, -307)) {
						E[i][j] = -706;
						//printf(" %.1f ,", E[i][j]);
					}
					else {
						E[i][j] = log((1 + temp) / (1 - temp));
						//printf(" %.1f ,", E[i][j]);
					}
				}
				else {
					E[i][j] = 0;
					//printf(" %.1f ,", E[i][j]);
				}
			}
			//printf("\n");
		}
		//printf("\n\n");

		//printf("L :\n");
		for (int j = 0; j < C; j++) {				//處理L陣列
			int i = 0;
			L[j] = r[j];
			while (i < R) {
				L[j] += E[i++][j];
			}
			//printf(" %.1f ,", L[j]);
		}
		//printf("\n\n");

		/*
		for (int i = 0; i < C; i++) {			//求原始訊號
			if (L[i] >= 0) decoding_info.estimated_codeword.at(i) = 0;			//為何? 比照調變時的maping
			else decoding_info.estimated_codeword.at(i) = 1;
		}
		for (int i = 0; i < C; i++) {
			//printf(" %d ,", res[i]);
		}
		//printf("\n\n");
		int *check = new int[R], iteration_check = 0;
		for (int i = 0; i < R; i++) {		//檢查陣列是否正確
			int temp = 0;
			for (int j = 0; j < C; j++) {
				temp ^= decoding_info.estimated_codeword.at(j) * H[i][j];
			}
			check[i] = temp;
			//printf(" %d ,", check[i]);
			iteration_check += check[i];
		}
		//printf("\n\n");
		if (iteration_check == 0) break;	//對M矩陣做迭代
		*/
		I_max++;
		//delete[] check;
	} while (I_max < 10);
	//distribution
	
	int tmp;
	for (int i = 0; i < G.Col_number; i++) {
		if (L[i] > 0) {
			tmp = (int)(L[i] * 10);
			LLR_dis[tmp + 500] += 1;
		}
		else {
			tmp = (int)(L[i] * 10);
			LLR_dis[tmp + 499] += 1;//向下累積
		}
		Mean += L[i]/(1000000*C);
		Var += (L[i] - 13.0218)*
			(L[i] - 13.0218) /
			(1000000 * C - 1);
	}
	if (count >= 1000000) {
		WriteFile("LLR.txt", "LLR_distribution");
		WriteFile("LLR.txt", "-----------");
		for (int i = 0; i < 1001; i++) {
			WriteFile("LLR.txt", (double)LLR_dis[i] / ((double)count*C));
		}
		WriteFile("LLR.txt", "-----------");
		WriteFile("LLR.txt", "Mean");
		WriteFile("LLR.txt", Mean);
		WriteFile("LLR.txt", "Var");
		WriteFile("LLR.txt", Var);
	}

	delete[] r;
	delete[] L;
	delete[] M;
	delete[] E;


	
}

void SPA_LLR1(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	static double Avg_LLR[100];
	static int count;
	double temp_LLR = 0;
	count++;
	ClearFile("LLR.txt");
	vector<double> R_temp;
	//convert all Rx to x = 1 distribution
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.code_seq.at(i) == 1) {//1 mapping to -1
			decoding_info.rx_signal_seq.at(i) += 2;
		}
	}
	if (decoding_info.Return_Est_1_or_LLR_0 == 0) R_temp = decoding_info.rx_signal_seq;
	// Change into Channel LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) decoding_info.rx_signal_seq.at(i) *= (2 / decoding_info.var);
	//建立H矩陣
	MATRIX<__int8> H;
	H.Building_Empty_Matrix(G.Col_number - G.Row_number, G.Col_number);
	H._matrix = G._matrix_outer;

	vector < vector <double> > M, Eji;
	vector <double> Li, zi;
	int I_count = 0, Col = H.Col_number, Row = H.Row_number;
	double m_temp;
	long double temp;
	vector<__int8> Current_Result(Col, 0), New_Result(Col, 0);
	//計算初始LLR
	temp_LLR = 0;
	for (int i = 0; i < G.Row_number; i++) {
		temp_LLR += abs(decoding_info.rx_signal_seq.at(i));
	}
	Avg_LLR[0] += temp_LLR / G.Row_number;

	//LLR distribution ,range in -100 ~ +100
	static int LLR_dis[1001];
	static double Var,Mean;
	//
	/*int tmp;
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.rx_signal_seq.at(i) > 0) {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 500] += 1;
		}
		else {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 499] += 1;//向下累積
		}
		Mean += decoding_info.rx_signal_seq.at(i) / (1000000 * Col);
		Var += (decoding_info.rx_signal_seq.at(i) - 3.989)*
			(decoding_info.rx_signal_seq.at(i) - 3.989) /
			(1000000 * Col - 1);
	}*/
	//
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
		/*
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
		*/
		if (I_count >= SPA_Iter) {
			if (/*decoding_info.Return_Est_1_or_LLR_0 == */0) {      // Return Hard Rx
				for (int i = 0; i < Col; ++i) {
					if (Li.at(i) > 0) decoding_info.estimated_codeword.at(i) = 0;
					else decoding_info.estimated_codeword.at(i) = 1;
				}
			}
			else {                                               // Return LLR
				//decoding_info.rx_signal_seq = R_temp;
				//cout << decoding_info.rx_signal_seq.size() << "," << Col << endl;
				for (int i = 0; i < Col; ++i) {
					decoding_info.rx_signal_seq.at(i) = Li.at(i)/**(decoding_info.var / 2)*/;
				}
				//cout << "!";
			}
			//cout << "Failed!" << endl;
			break;
		}
		// End

		I_count++;
		//計算平均LLR
		temp_LLR = 0;
		for (int i = 0; i < G.Row_number; i++) {
			temp_LLR += abs(Li.at(i));
		}
		Avg_LLR[I_count] += temp_LLR / G.Row_number;

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
	//distribution
	
	int tmp;
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.rx_signal_seq.at(i) > 0) {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 500] += 1;
		}
		else {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 499] += 1;//向下累積
		}
		Mean += decoding_info.rx_signal_seq.at(i)/(1000000*Col);
		Var += (decoding_info.rx_signal_seq.at(i) - 12.8269)*
			(decoding_info.rx_signal_seq.at(i) - 12.8269) /
			(1000000 * Col - 1);
	}
	if (count >= 1000000) {
		WriteFile("LLR.txt", "LLR_distribution");
		WriteFile("LLR.txt", "-----------");
		for (int i = 0; i < 1001; i++) {
			WriteFile("LLR.txt", (double)LLR_dis[i] / ((double)count*Col));
		}
		WriteFile("LLR.txt", "-----------");
		WriteFile("LLR.txt", "Mean");
		WriteFile("LLR.txt", Mean);
		WriteFile("LLR.txt", "Var");
		WriteFile("LLR.txt", Var);
	}
	/*
	WriteFile("LLR.txt", "Avg_LLR");
	WriteFile("LLR.txt", "-----------");
	for (int i = 0; i < 51; i++) {
		WriteFile("LLR.txt", Avg_LLR[i]/count);
	}
	WriteFile("Hard_test.txt", "-----------");
	*/
}

void SPA_Astar(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	MATRIX<__int8> H;
	H._matrix = G._matrix_outer;
	vector<int> permutaion_seq;
	H.Col_number = G.Col_number;
	H.Row_number = G.Col_number - G.Row_number;
	
	MATRIX<__int8> tmp_G;
	rref(H,tmp_G);
	for (int i = 0; i < G.Row_number; i++) {
		for (int j = 0; j < G.Col_number; j++) {
			cout << (int)tmp_G._matrix[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

//Brute Force Soft-in/Soft-out
//only use short code
void BF_SISO(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	//convert all Rx to x = 1 distribution
	/*
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.code_seq.at(i) == 1) {//1 mapping to -1
			decoding_info.rx_signal_seq.at(i) += 2;
		}
	}*/
	//
	int message_length = G.Row_number;
	int codeword_length = G.Col_number;
	int status = pow(2, message_length);
	vector<double> table(status);
	static bool flag = false;
	static vector<vector<double>> codeword(status, vector<double>(codeword_length));
	//initial codeword
	if (!flag) {
		vector<vector<__int8>> tmp_codeword(status, vector<__int8>(codeword_length));
		vector<__int8> message_seq(message_length);
		for (int i = 0; i < status; i++) {
			for (int j = 0; j < message_length; j++) {
				message_seq.at(j) = (i & (1 << j)) > 0 ? 1 : 0;
			}
			Systematic_Linear_Block_Code_Encoder(G, message_seq, tmp_codeword.at(i));
			for (int j = 0; j < codeword_length; j++) {
				codeword.at(i).at(j) = (double)(tmp_codeword.at(i).at(j) == 0 ? 1 : -1);
			}
		}
		flag = true;
	}
	//Compute P(y|c) = exp(-||y - c||^2/(2*sigma^2))
	for (int i = 0; i < status; i++) {
		double tmp = 0;
		double tmp2 = 0;
		for (int j = 0; j < codeword_length; j++) {
			tmp2 = (decoding_info.rx_signal_seq.at(j) - codeword.at(i).at(j) );
			tmp += tmp2 * tmp2;
		}
		tmp = tmp * (-1) / (2 * decoding_info.var);
		table.at(i) = exp(tmp);
	}
	//Compute LLR value
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); i++) {
		double tmp_0 = 0, tmp_1 = 0;
		for (int j = 0; j < status; j++) {
			if (codeword.at(j).at(i) == 1) {
				tmp_0 += table.at(j);
			}
			else {
				tmp_1 += table.at(j);
			}
		}
		decoding_info.rx_signal_seq.at(i) = log(tmp_0 / tmp_1);
	}/*
	for (int i = 0; i < codeword_length; ++i) {
		cout << decoding_info.rx_signal_seq.at(i) << ",";
	}
	cout << endl;*/
	//Hard desicion(optional)
	double temp_LLR = 0;
	for (int i = 0; i < codeword_length; ++i) {
		if (decoding_info.rx_signal_seq.at(i) > 0) decoding_info.estimated_codeword.at(i) = 0;
		else decoding_info.estimated_codeword.at(i) = 1;
	}
	/*
	//LLR distribution ,range in -100 ~ +100
	static int LLR_dis[1001];
	static double Var = 0, Mean = 0;
	static int count = 0;
	count++;
	int tmp,Col = G.Col_number;
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.rx_signal_seq.at(i) > 0) {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 500] += 1;
		}
		else {
			tmp = (int)(decoding_info.rx_signal_seq.at(i) * 10);
			LLR_dis[tmp + 499] += 1;//向下累積
		}
		Mean += decoding_info.rx_signal_seq.at(i) / (1000000 * Col);
		Var += (decoding_info.rx_signal_seq.at(i) - 14.7704)*
			(decoding_info.rx_signal_seq.at(i) - 14.7704) /
			(1000000 * Col - 1);
	}
	if (count >= 1000000) {
		ClearFile("LLR.txt");
		WriteFile("LLR.txt", "LLR_distribution");
		WriteFile("LLR.txt", "-----------");
		for (int i = 0; i < 1001; i++) {
			WriteFile("LLR.txt", (double)LLR_dis[i] / ((double)count*Col));
		}
		WriteFile("LLR.txt", "-----------");
		WriteFile("LLR.txt", "Mean");
		WriteFile("LLR.txt", Mean);
		WriteFile("LLR.txt", "Var");
		WriteFile("LLR.txt", Var);
	}*/
}

//Brute Force + A*
//only use (8,4)Hamming code
void BF_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_inner._matrix = G._H;
	H_inner.Col_number = G_.Col_number;
	H_inner.Row_number = G_.Col_number - G_.Row_number;
	//for BF
	G__._matrix = G.subG;
	G__.Col_number = G.subG.at(0).size();
	G__.Row_number = G.subG.size();

	int i, j, start;
	vector<double>
		subRx(G__.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number*(G__.Row_number) / G__.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);

	for (i = 0; i < G.Col_number / G__.Col_number; i++) {
		//assign sub-Block
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - G__.Row_number);
		}
		//Brute-Force SISO
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			Rx_parity.at(start + j - G__.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;
	//test
	//output of SISO is extrinsic LLR
	/*
	for (int i = G_.Row_number; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq[i] = tmp_Rx[i] * (2 / decoding_info.var);
	}
	*/
	//original A*
	A_star_2_stack_CBC_OSC(G_, decoding_info);
	//RS_Hard_desicion(G_, decoding_info);
	//set parity LLR estimate
	decoding_info.estimated_codeword.resize(G.Col_number);

	for (i = G_.Col_number; i < G.Col_number; i++) {
		if (Rx_parity.at(i - G_.Col_number) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = tmp_Rx;
	/*
	//test
	decoding_info.estimated_codeword.resize(G_.Col_number);
	static long err_loc[128];
	static int count = 0;
	count++;
	vector<__int8> Sorted_codeword(G_.Col_number, 0), Sorted_estimated(G_.Col_number, 0);
	vector <size_t>
		Location_Index(G_.Col_number, 0);
	MATRIX<__int8> Sorted_G(G_);
	MATRIX<double> Metric_Table(2, G.Col_number);
	Pre_Procedure(decoding_info.rx_signal_seq, G_, Sorted_G, Location_Index, Metric_Table);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	Sort_Function(decoding_info.estimated_codeword, Location_Index, Sorted_estimated);
	for (int i = 0; i < 36; i++) {
		if (Sorted_estimated[i] != Sorted_codeword[i]) {
			err_loc[i]++;
		}
	}
	if (count >= 5000000) {
		ClearFile("Error_Loc.txt");
		for (int i = 0; i < 36; i++) {
			WriteFile("Error_Loc.txt", (double)err_loc[i] / (double)count);
		}
	}
	decoding_info.estimated_codeword.resize(G.Col_number);
	return;*/
	//
}

void BF_subDecoder(MATRIX<__int8>& G,vector<double>& subRx,double& var) {
	int message_length = G.Row_number;
	int codeword_length = G.Col_number;
	int status = pow(2, message_length);
	vector<double> table(status);
	static bool flag = false;
	static vector<vector<double>> codeword(status, vector<double>(codeword_length));
	//initial codeword
	if (!flag) {
		vector<vector<__int8>> tmp_codeword(status, vector<__int8>(codeword_length));
		vector<__int8> message_seq(message_length);
		for (int i = 0; i < status; i++) {
			for (int j = 0; j < message_length; j++) {
				message_seq.at(j) = (i & (1 << j)) > 0 ? 1 : 0;
			}
			Systematic_Linear_Block_Code_Encoder(G, message_seq, tmp_codeword.at(i));
			for (int j = 0; j < codeword_length; j++) {
				codeword.at(i).at(j) = (double)(tmp_codeword.at(i).at(j) == 0 ? 1 : -1);
			}
		}
		flag = true;
	}
	//Compute P(y|c) = exp(-||y - c||^2/(2*sigma^2))
	for (int i = 0; i < status; i++) {
		double tmp = 0;
		double tmp2 = 0;
		for (int j = 0; j < codeword_length; j++) {
			tmp2 = (subRx.at(j) - codeword.at(i).at(j));
			tmp += tmp2 * tmp2;
		}
		tmp = tmp * (-1) / (2 * var);
		table.at(i) = exp(tmp);
	}
	//Compute LLR value
	for (int i = 0; i < subRx.size(); i++) {
		double tmp_0 = 0, tmp_1 = 0;
		for (int j = 0; j < status; j++) {
			if (codeword.at(j).at(i) == 1) {
				tmp_0 += table.at(j);
			}
			else {
				tmp_1 += table.at(j);
			}
		}
		subRx.at(i) = log(tmp_0 / tmp_1);
	}/*
	for (int i = 0; i < codeword_length; ++i) {
		cout << subRx.at(i) << ",";
	}
	cout << endl;*/
}

void BF_A_star_ver22(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H, G_,G__;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
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
		Rx(G.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq;

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
		//SPA iteration
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (H.Row_number);
		for (; j < H.Col_number; j++) {
			Rx.at(G_.Col_number + start + j - H.Col_number + H.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;

	A_star_I(G, decoding_info);
	decoding_info.rx_signal_seq = tmp_Rx;
}
//改
void BF_A_star_ver2(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static bool flag = false;
	static MATRIX<__int8> H, G_, G__;
	static vector<vector<double*>>
		subRx(G.Col_number / G._matrix_outer.at(0).size(),vector<double*>(G._matrix_outer.at(0).size()));
	vector<double> tmp_Rx = decoding_info.rx_signal_seq;
	if (!flag) {
		H._matrix = G._matrix_outer;
		H.Col_number = H._matrix.at(0).size();
		H.Row_number = H._matrix.size();
		G_._matrix = G._matrix_inner;
		G_.Col_number = G._matrix_inner.at(0).size();
		G_.Row_number = G._matrix_inner.size();
		//for BF
		G__.Building_Empty_Matrix(H.Col_number - H.Row_number, H.Col_number);
		G__._matrix = {
		{1,0,0,0,1,1,1,0},
		{0,1,0,0,1,1,0,1},
		{0,0,1,0,1,0,1,1},
		{0,0,0,1,0,1,1,1}
		};
		int start;
		for (int i = 0; i < G.Col_number / H.Col_number; i++) {
			//message
			start = i * (H.Col_number - H.Row_number);
			for (int j = 0;j < H.Col_number - H.Row_number;j++) {
				subRx.at(i).at(j) = &(decoding_info.rx_signal_seq.at(start + j));
			}
			//check
			start = start = i * H.Row_number;
			for (int j = 0; j < H.Row_number; j++) {
				subRx.at(i).at(H.Col_number - H.Row_number + j) =
					&(decoding_info.rx_signal_seq.at(G_.Col_number + start + j));
			}
		}
		//finish 
		flag = true;
	}

	// SISO
	for (int i = 0; i < G.Col_number / H.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, subRx.at(i), decoding_info.var);
		
	}
	//Change flipping
	for (int i = 0; i < G.Col_number; i++) {
		if (cmp_sign(tmp_Rx[i], decoding_info.rx_signal_seq[i])) {
			decoding_info.rx_signal_seq[i] += tmp_Rx[i] * (2 / decoding_info.var);
		}
	}
	
	//original A*
	//A_star_I(G, decoding_info);
	A_star_3_stack(G, decoding_info);
	decoding_info.rx_signal_seq = tmp_Rx;
}

void pointer_BF_subDecoder(MATRIX<__int8>& G, vector<double*>& subRx, double& var) {
	int message_length = G.Row_number;
	int codeword_length = G.Col_number;
	int status = pow(2, message_length);
	vector<double> table(status);
	static bool flag = false;
	static vector<vector<double>> codeword(status, vector<double>(codeword_length));
	//initial codeword
	if (!flag) {
		vector<vector<__int8>> tmp_codeword(status, vector<__int8>(codeword_length));
		vector<__int8> message_seq(message_length);
		for (int i = 0; i < status; i++) {
			for (int j = 0; j < message_length; j++) {
				message_seq.at(j) = (i & (1 << j)) > 0 ? 1 : 0;
			}
			Systematic_Linear_Block_Code_Encoder(G, message_seq, tmp_codeword.at(i));
			for (int j = 0; j < codeword_length; j++) {
				codeword.at(i).at(j) = (double)(tmp_codeword.at(i).at(j) == 0 ? 1 : -1);
			}
		}
		flag = true;
	}
	//Compute P(y|c) = exp(-||y - c||^2/(2*sigma^2))
	for (int i = 0; i < status; i++) {
		double tmp = 0;
		double tmp2 = 0;
		for (int j = 0; j < codeword_length; j++) {
			tmp2 = (*(subRx.at(j)) - codeword.at(i).at(j));
			tmp += tmp2 * tmp2;
		}
		tmp = tmp * (-1) / (2 * var);
		table.at(i) = exp(tmp);
	}
	//Compute LLR value
	for (int i = 0; i < subRx.size(); i++) {
		double tmp_0 = 0, tmp_1 = 0;
		for (int j = 0; j < status; j++) {
			if (codeword.at(j).at(i) == 1) {
				tmp_0 += table.at(j);
			}
			else {
				tmp_1 += table.at(j);
			}
		}
		*(subRx.at(i)) = log(tmp_0 / tmp_1);
	}/*
	for (int i = 0; i < codeword_length; ++i) {
		cout << *(subRx.at(i)) << ",";
	}
	cout << endl;*/
}
void pointer_BF_subDecoder2(MATRIX<__int8>& G, vector<double*>& subRx, double& var) {
	int message_length = G.Row_number;
	int codeword_length = G.Col_number;
	int status = pow(2, message_length);
	vector<double> table(status);
	static bool flag = false;
	static vector<vector<double>> codeword(status, vector<double>(codeword_length));
	//initial codeword
	if (!flag) {
		vector<vector<__int8>> tmp_codeword(status, vector<__int8>(codeword_length));
		vector<__int8> message_seq(message_length);
		for (int i = 0; i < status; i++) {
			for (int j = 0; j < message_length; j++) {
				message_seq.at(j) = (i & (1 << j)) > 0 ? 1 : 0;
			}
			Systematic_Linear_Block_Code_Encoder(G, message_seq, tmp_codeword.at(i));
			for (int j = 0; j < codeword_length; j++) {
				codeword.at(i).at(j) = (double)(tmp_codeword.at(i).at(j) == 0 ? 1 : -1);
			}
		}
		flag = true;
	}
	//Compute P(y|c) = exp(-||y - c||^2/(2*sigma^2))
	for (int i = 0; i < status; i++) {
		double tmp = 0;
		double tmp2 = 0;
		for (int j = 0; j < codeword_length; j++) {
			tmp2 = (*(subRx.at(j)) - codeword.at(i).at(j));
			tmp += tmp2 * tmp2;
		}
		tmp = tmp * (-1) / (2 * var);
		table.at(i) = exp(tmp);
	}
	//Compute LLR value
	for (int i = 0; i < subRx.size(); i++) {
		double tmp_0 = 0, tmp_1 = 0;
		for (int j = 0; j < status; j++) {
			if (codeword.at(j).at(i) == 1) {
				tmp_0 += table.at(j);
			}
			else {
				tmp_1 += table.at(j);
			}
		}
		*(subRx.at(i)) = log(tmp_0 / tmp_1);
	}/*
	for (int i = 0; i < codeword_length; ++i) {
		cout << *(subRx.at(i)) << ",";
	}
	cout << endl;*/
}
//Brute Force + A*
//only use (8,4)Hamming code
void BF_A_star_interleaver(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H, G_, G__;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
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
			subRx.at(j) = decoding_info.rx_signal_seq.at(G.interleaver.at(start + j));
		}
		start = i * H.Row_number;
		for (; j < H.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - H.Col_number + H.Row_number);
		}
		//Brute-Force SISO
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			Rx.at(G.interleaver.at(start + j)) = subRx.at(j);
		}
		start = i * (H.Row_number);
		for (; j < H.Col_number; j++) {
			Rx_parity.at(start + j - H.Col_number + H.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;
	//test


	//output of SISO is extrinsic LLR
	/*
	for (int i = 0; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq[i] += tmp_Rx[i] * (2 / decoding_info.var);
	}
	*/
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


//add interleaver
void BF_A_star_interleaver_ver2(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static bool flag = false;
	static MATRIX<__int8> H, G_, G__;
	static vector<vector<double*>>
		subRx(G.Col_number / G._matrix_outer.at(0).size(), vector<double*>(G._matrix_outer.at(0).size()));
	vector<double> tmp_Rx = decoding_info.rx_signal_seq;
	if (!flag) {
		H._matrix = G._matrix_outer;
		H.Col_number = H._matrix.at(0).size();
		H.Row_number = H._matrix.size();
		G_._matrix = G._matrix_inner;
		G_.Col_number = G._matrix_inner.at(0).size();
		G_.Row_number = G._matrix_inner.size();
		//for BF
		G__.Building_Empty_Matrix(H.Col_number - H.Row_number, H.Col_number);
		G__._matrix = {
		{1,0,0,0,1,1,1,0},
		{0,1,0,0,1,1,0,1},
		{0,0,1,0,1,0,1,1},
		{0,0,0,1,0,1,1,1}
		};
		int start;
		for (int i = 0; i < G.Col_number / H.Col_number; i++) {
			//message
			start = i * (H.Col_number - H.Row_number);
			for (int j = 0; j < H.Col_number - H.Row_number; j++) {
				subRx.at(i).at(j) = &(decoding_info.rx_signal_seq.at(G.interleaver.at(start + j)));
			}
			//check
			start = start = i * H.Row_number;
			for (int j = 0; j < H.Row_number; j++) {
				subRx.at(i).at(H.Col_number - H.Row_number + j) =
					&(decoding_info.rx_signal_seq.at(G_.Col_number + start + j));
			}
		}
		//finish 
		flag = true;
	}

	// SISO
	for (int i = 0; i < G.Col_number / H.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, subRx.at(i), decoding_info.var);

	}

	//original A*
	A_star_I(G, decoding_info);

	decoding_info.rx_signal_seq = tmp_Rx;
}

inline bool cmp_sign(double a,double b) {
	if (a > 0 && b < 0) return true;
	if (a < 0 && b > 0) return true;
	return false;
}


//Brute Force + SPA + A*
//only use (8,4)Hamming code
void BF_SPA_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8>  G_, G__,H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_inner._matrix = G._H;
	H_inner.Col_number = G_.Col_number;
	H_inner.Row_number = G_.Col_number -G_.Row_number;
	//for BF
	G__._matrix = G.subG;
	G__.Col_number = G.subG.at(0).size();
	G__.Row_number = G.subG.size();

	int i, j, start;
	vector<double>
		subRx(G__.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number*(G__.Row_number) / G__.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);

	for (i = 0; i < G.Col_number / G__.Col_number; i++) {
		//assign sub-Block
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - G__.Row_number);
		}
		//Brute-Force SISO
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			Rx_parity.at(start + j - G__.Row_number) = subRx.at(j);
		}
	}
	decoding_info.rx_signal_seq = Rx;
	/*for (i = 0; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq.at(i) /= (2 / decoding_info.var);
	}*/
	//SPA_subDecoder(H_inner, decoding_info.rx_signal_seq, decoding_info.var, 50);
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

void RS_Hard_desicion(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	int t = (G.Col_number - G.Row_number) / 8,count_err = 0;
	vector<bool> err(G.Col_number/4, 0);
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.rx_signal_seq.at(i) < 0) {
			decoding_info.estimated_codeword.at(i) = 1;
		}
		else decoding_info.estimated_codeword.at(i) = 0;
	}
	//decoding
	for (int i = 0; i < G.Col_number; i++) {
		if (decoding_info.estimated_codeword.at(i) != decoding_info.code_seq.at(i)) {
			if (!err[i / 4]) {
				err[i / 4] = 1;
				count_err++;
			}
		}
	}
	if (count_err <= t) {
		for (int i = 0; i < G.Col_number; i++) {
			decoding_info.estimated_codeword.at(i) = decoding_info.code_seq.at(i);
		}
	}
}

//Sorted all codeword in MRIP
void MRIP_BF_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_inner._matrix = G._H;
	H_inner.Col_number = G_.Col_number;
	H_inner.Row_number = G_.Col_number - G_.Row_number;
	//for BF
	G__._matrix = G.subG;
	G__.Col_number = G.subG.at(0).size();
	G__.Row_number = G.subG.size();

	int i, j, start;
	vector<double>
		subRx(G__.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number*(G__.Row_number) / G__.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);

	for (i = 0; i < G.Col_number / G__.Col_number; i++) {
		//assign sub-Block
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - G__.Row_number);
		}
		//Brute-Force SISO
		BF_subDecoder(G__, subRx, decoding_info.var);
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			Rx_parity.at(start + j - G__.Row_number) = subRx.at(j);
		}
	}
	//decoding_info.rx_signal_seq = Rx;
	//Sorted codeword
	vector<double> Rx_signal(G.Col_number, 0);
	for (int i = 0; i < G_.Col_number; i++) {
		Rx_signal.at(i) = Rx.at(i);
	}
	for (int i = G_.Col_number; i < G.Col_number; i++) {
		Rx_signal.at(i) = Rx_parity.at(i - G_.Col_number);
	}
	vector<size_t> Permutation(G.Col_number, 0);
	MATRIX<__int8> Sorted_G(G);
	Determine_Permutation(Rx_signal, G, Sorted_G, Permutation);
	Sort_Function(Rx_signal, Permutation, decoding_info.rx_signal_seq);
	//save parity
	/*
	for (int i = G_.Col_number; i < G.Col_number; i++) {
		Rx_parity.at(i) = decoding_info.rx_signal_seq.at(i);
	}*/
	//
	decoding_info.rx_signal_seq.resize(G_.Col_number);
	//deal with matrix
	for (int i = 0; i < G_.Row_number; i++) {
		for (int j = 0; j < G_.Col_number; j++) {
			G_._matrix[i][j] = Sorted_G._matrix[i][j];
		}
	}
	//original A*
	A_star_I(G_, decoding_info);
	//message bits
	vector<__int8> rx(G.Row_number, 0), codeword_seq(G.Col_number,0);
	for (int i = 0; i < G.Row_number; i++) {
		rx.at(i) = decoding_info.estimated_codeword.at(i);
	}
	decoding_info.estimated_codeword.resize(G.Col_number);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, rx, codeword_seq);
	Desort_Function(Permutation, codeword_seq, decoding_info.estimated_codeword);
	decoding_info.rx_signal_seq = tmp_Rx;
}

void Subcode_MRIP_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8>  G_, G__, H_inner,H;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_inner._matrix = G._H;
	H_inner.Col_number = G_.Col_number;
	H_inner.Row_number = G_.Col_number - G_.Row_number;
	//for BF
	G__._matrix = G.subG;
	G__.Col_number = G.subG.at(0).size();
	G__.Row_number = G.subG.size();
	Convert_SystematicG_to_H(H, G__);
	int i, j, start;
	vector<double>
		subRx(G__.Col_number, 0),	//給每個sub block用 
		Rx(G.Col_number*(G__.Row_number) / G__.Col_number),
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);

	for (i = 0; i < G.Col_number / G__.Col_number; i++) {
		//assign sub-Block
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(start + j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			subRx.at(j) = decoding_info.rx_signal_seq.at(G_.Col_number + start + j - G__.Row_number);
		}
		//Brute-Force SISO
		//BF_subDecoder(G__, subRx, decoding_info.var);
		ASPA_subdecoder(H, subRx, decoding_info.var, 5);
		start = i * (G__.Row_number);
		for (j = 0; j < G__.Row_number; j++) {
			Rx.at(start + j) = subRx.at(j);
		}
		start = i * (G__.Col_number - G__.Row_number);
		for (; j < G__.Col_number; j++) {
			Rx_parity.at(start + j - G__.Row_number) = subRx.at(j);
		}
	}
	
	for (i = 0; i < G_.Col_number; i++) {
		if (Rx.at(i) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = Rx;
	//original A*
	//A_star_I(G_, decoding_info);
	//A_star_3_stack(G_,decoding_info);
	//decoding_info.STE *= G_.Row_number;
	//decoding_info.COM *= G_.Row_number;
	//set parity LLR estimate
	
	decoding_info.estimated_codeword.resize(G.Col_number);

	for (i = G_.Col_number; i < G.Col_number; i++) {
		if (Rx_parity.at(i - G_.Col_number) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//decoding_info.rx_signal_seq = tmp_Rx;
	decoding_info.rx_signal_seq.resize(G.Col_number);
	for (i = G_.Col_number; i < G.Col_number; i++) {
		decoding_info.rx_signal_seq.at(i) = Rx_parity.at(i - G_.Col_number);
	}
	/*
	//test
	static long err_count[128];
	static long err_loc[128];
	static int count = 0;
	int err_bit = 0;
	count++;
	int idx1,idx2,idx3,idx4;
	double max1, max2, max3 = FLT_MIN, max4 = FLT_MIN;
	vector<pair<double,int>> MRIP;
	i = 0;
	while (i < G.Col_number / G__.Col_number) {
		int start1 = i * (G__.Row_number);
		int start2 = G_.Col_number + i * (G__.Col_number - G__.Row_number);
		max1 = max2 = FLT_MIN;
		for (int i = 0; i < G__.Row_number; i++) {
			if (abs(decoding_info.rx_signal_seq.at(start1 + i)) > abs(max1)) {
				max2 = max1; idx2 = idx1;
				max1 = decoding_info.rx_signal_seq.at(start1 + i); idx1 = start1+i;
			}
			else if (abs(decoding_info.rx_signal_seq.at(start1+i)) > abs(max2)) {
				max2 = decoding_info.rx_signal_seq.at(start1 + i); idx2 = start1+i;
			}
		}
		for (int i = 0; i < G__.Col_number - G__.Row_number; i++) {
			if (abs(decoding_info.rx_signal_seq.at(start2+i)) > abs(max1)) {
				max2 = max1; idx2 = idx1;
				max1 = decoding_info.rx_signal_seq.at(start2 + i); idx1 = start2+i;
			}
			else if (abs(decoding_info.rx_signal_seq.at(start2 + i)) > abs(max2)) {
				max2 = decoding_info.rx_signal_seq.at(start2 + i); idx2 = start2 + i;
			}
		}
		MRIP.push_back(make_pair(max1, idx1));
		MRIP.push_back(make_pair(max2, idx2));
		i++;
	}
	//放入第一個block的四個bit
	for (int i = 0; i < G__.Row_number; i++) {
		if (MRIP.size() == G.Row_number) break;
		if (i == MRIP[0].second || i == MRIP[1].second) continue;
		MRIP.push_back(make_pair(decoding_info.rx_signal_seq.at(i), i));
	}
	for (int i = G_.Col_number; i <G_.Col_number + G__.Col_number - G__.Row_number; i++) {
		if (MRIP.size() == G.Row_number) break;
		if (i == MRIP[0].second || i == MRIP[1].second) continue;
		MRIP.push_back(make_pair(decoding_info.rx_signal_seq.at(i), i));
	}
	sort(MRIP.begin(),MRIP.end(),
		[](const pair<double, int>& a, const pair<double, int>& b) {
		return abs(a.first) > abs(b.first);
	});
	int err = 0;
	for (int i = 0; i < G.Row_number; i++) {
		if (MRIP[i].first < 0 && decoding_info.code_seq.at(MRIP[i].second) == 0 ||
			MRIP[i].first > 0 && decoding_info.code_seq.at(MRIP[i].second) == 1) {
			err++;
			err_loc[i]++;
		}
	}
	err_count[err]++;
	if (count >= 5000000) {
		ClearFile("Error_Count.txt");
		for (int i = 0; i < G.Row_number; i++) {
			WriteFile("Error_Count.txt", (double)err_count[i] / (double)count);
		}
		ClearFile("Error_Loc.txt");
		for (int i = 0; i < G.Row_number; i++) {
			WriteFile("Error_Loc.txt", (double)err_loc[i] / (double)count);
		}
	}
	return;*/
	
	/*
	//test
	static long err_count[128];
	static long err_loc[128];
	static int count = 0;
	MATRIX<__int8> S_G(G);
	vector<__int8> RX(G.Col_number, 0),CX(G.Col_number, 0),MX(G.Row_number,0);
	int err_bit = 0;
	count++;
	vector <size_t>
		Permutation(G.Col_number, 0);
	vector<pair<int,double>> value_sum(G.Col_number/2);
	for (int i = 0; i < G.Col_number/2; i++) {
		value_sum[i] = make_pair(i, abs(decoding_info.rx_signal_seq[2 * i]) + abs(decoding_info.rx_signal_seq[2 * i + 1]));
	}
	sort(value_sum.begin(), value_sum.end(),
		[](const pair<int, double>& a,const pair<int,double>& b) {
		return a.second > b.second;
	});
	//Systematic form
	for (int i = 0; i < G.Col_number / 2; i++) {
		Permutation[2 * i] = value_sum[i].first;
		Permutation[2 * i+1] = value_sum[i].first+1;
	}
	for (int i = 0; i < G.Col_number; i++) {
		decoding_info.estimated_codeword[i] = decoding_info.rx_signal_seq[i] > 0 ? 0 : 1;
	}
	//Determine_Permutation_Matrix(G,S_G,Permutation);
	Sort_Function(decoding_info.estimated_codeword, Permutation, RX);
	Sort_Function(decoding_info.code_seq, Permutation, CX);
	for (int i = 0; i < G.Row_number ; i++) {
		if (RX[i] != CX[i]) {
			err_bit++;
		}
	}
	err_count[err_bit]++;
	if (count >= 100000) {
		ClearFile("Error_Count.txt");
		for (int i = 0; i < G.Row_number; i++) {
			WriteFile("Error_Count.txt", (double)err_count[i] / (double)count);
		}
		ClearFile("Error_Loc.txt");
		for (int i = 0; i < G.Row_number; i++) {
			WriteFile("Error_Loc.txt", (double)err_loc[i] / (double)count);
		}
	}
	
	MX.assign(RX.begin(), RX.begin() + G.Row_number);
	Systematic_Linear_Block_Code_Encoder(S_G, MX, CX);
	Desort_Function(Permutation, CX, decoding_info.estimated_codeword);
	return;*/
	//multiple stack
	/*
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

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	// use SISO_Astar as inital path
	//for (size_t i(0); i < codeword_length; ++i)
		//if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	decoding_info.Constraint_i = 4;
	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
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
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

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
			//A_star_2_stack(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info);
			A_star_MultipleStack_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
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
	*/

	//using queue without Comparison
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
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(decoding_info.rx_signal_seq);
		Determine_Permutation(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index);
		Sort_Function(tmp_Rx, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
	} while (0);
	decoding_info.rx_signal_seq = tmp_Rx;
	//test
	/*
	double metric_sum = 0;
	for (int i = 0; i < codeword_length; i++)
		metric_sum += abs(tmp_Rx.at(i));
	static int dis_zero[1000] = {};
	static int dis_err[1000] = {};
	static long err_count[128];
	static int count = 0;
	static int ccount = 0;
	int err_bit = 0;
	count++;
	vector<__int8> Sorted_codeword(codeword_length, 0);
	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	//如果在一半的MRIP中錯兩個bits，接下來錯誤的機率
	int index = 0;
	vector<int> err_idx;
	for (; index < message_length; index++) {
		if (Sorted_codeword[index] != Hard_RX[index]) {
			err_bit++;
			err_idx.push_back(index);
		}
		if (err_bit == 2) {
			err_bit = 0;
			ccount++;
			index++;
			for (; index < message_length; index++) {
				if (Sorted_codeword[index] != Hard_RX[index]) {
					err_bit++;
				}
			}
			err_count[err_bit]++;
			double metric = 0;
			vector<__int8> test_code = Hard_RX;
			for (int i = 0; i < 2; i++) {
				test_code[err_idx[i]] ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					test_code.at(j) ^= Sorted_G._matrix[err_idx[i]][j];
				}
			}
			for (int i = 0; i < codeword_length; i++) 
				metric += abs(Metric_Table._matrix[test_code[i]][i]);
			if (err_bit == 0) {
				dis_zero[(int)(metric/metric_sum * 1000)]++;
			}
			else {
				dis_err[(int)(metric/metric_sum * 1000)]++;
			}
			break;
		}
	}
	
	if (count >= 5000000) {
		ClearFile("Error_Count.txt");
		ClearFile("dis_zero.txt");
		ClearFile("dis_err.txt");
		for (int i = 0; i < message_length; i++) {
			WriteFile("Error_Count.txt", (double)err_count[i] / (double)ccount);
		}
		for (int i = 0; i < 1000; i++) {
			WriteFile("dis_zero.txt", dis_zero[i]);
			WriteFile("dis_err.txt", dis_err[i]);
		}
	}
	decoding_info.rx_signal_seq = tmp_Rx;
	return;*/
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
				tmp.level = i+1;
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

void A_star_MultipleStack(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
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

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
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
			decoding_info.CandidateCodeWord++;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

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
			A_star_MultipleStack_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info);
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
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

void A_star_MultipleStack_inner(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	do
	{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (Metric_Table._matrix[codeword_seq[j]][j]) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
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
			if (pc_i < decoding_info.Constraint_i - 1) {
				Best_Goal_2 = Best_Goal;
				A_star_MultipleStack_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 1, decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
			else if (pc_i == decoding_info.Constraint_i - 1) {
				Best_Goal_2 = Best_Goal;
				A_star_MultipleStack_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 2, decoding_info);
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
			else {
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					//for subcode SISO Astar using
					if (Metric_Table._matrix[Hard_RX.at(j)][j])
						Child_Node.metric += Metric_Table._matrix[Hard_RX.at(j)][j];
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這邊的G是Sorted_G
					}
				}
				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (Metric_Table._matrix[codeword_seq[j]][j]) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric >= Best_Goal.metric) break;
					}
				}
				Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				
			}
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
}

void A_star_MultipleStack_SC(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
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

	vector<NODE_PATH> Stack(1, Pointer);

	Best_Goal.metric = FLT_MAX;

	bool Update_Flag = false;
	double OSC_metric_thr(0);
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);

	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table,decoding_info);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);

	do {
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= 1)) {
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
				if (Metric_Table._matrix[codeword_seq.at(j)][j]) {
					Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
					if (Pointer.metric > Best_Goal.metric) break;
				}
				else Pointer.Same_Index.push_back(j);
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//
			Update_Flag = Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Update_Flag == TRUE) {
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) {
					break;
				}
				Temp_Sufficient_Condition_thr = 0;
			}
			
		}
		else if ((Pointer.level < message_length) && (Pointer.D_z < 1)) {

			decoding_info.STE += 2;
			decoding_info.Binary_STE += 2;
			for (__int8 new_bit(0); new_bit < 2; new_bit++) {

				Extend_Node_Procedure(Pointer, Child_Node, Metric_Table, new_bit);
				if (new_bit != Hard_RX.at(Pointer.level)) {
					++Child_Node.D_z;
					Child_Node.Diff_Index.push_back(Pointer.level);
				}
				else Child_Node.Same_Index.push_back(Pointer.level);

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
			if (A_star_MultipleStack_SC_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, (size_t)(1 + 1), decoding_info)) {
				if (Best_Goal.metric >= Best_Goal_2.metric) {
					Best_Goal = Best_Goal_2;
					break;
				}
			}
			Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
		}
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

bool A_star_MultipleStack_SC_inner(MATRIX<__int8> &Sorted_G, MATRIX<double> &Metric_Table, vector<__int8> &Hard_RX, vector<__int8> &MRIP_codeword, NODE_PATH &Node, NODE_PATH &Pre_Best_Goal, size_t pc_i, DECODING_INFO &decoding_info)
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
		Best_Goal_2(Pre_Best_Goal),
		Child_Node(message_length);
	vector<NODE_PATH> Stack(1, Node);

	bool Update_Flag = false;
	double OSC_metric_thr(0);
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);

	do
	{
		Pointer = Stack.at(0);

		if ((Pointer.level == message_length) && (Pointer.D_z <= pc_i)) {
			Stack.erase(Stack.begin());
			//
			//Systematic_Linear_Block_Code_Encoder(Sorted_G, Pointer.message_bits, codeword_seq);

			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}

			for (size_t j(message_length); j < codeword_length; ++j) {
				if (Metric_Table._matrix[codeword_seq[j]][j]) {
					Pointer.metric += Metric_Table._matrix[codeword_seq[j]][j];
					if (Pointer.metric >= Best_Goal.metric) break;
				}
				else Pointer.Same_Index.push_back(j);
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			++decoding_info.CandidateCodeWord;
			//
			Update_Flag = Update_Best_Goal_Procedure(Pointer, Best_Goal, Stack);
			if (Update_Flag == TRUE) {
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) {
					Pre_Best_Goal = Best_Goal;
					return true;
				}
				Temp_Sufficient_Condition_thr = 0;
			}
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
				else Child_Node.Same_Index.push_back(Pointer.level);

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
			if (pc_i < decoding_info.Constraint_i - 1) {
				Best_Goal_2 = Best_Goal;
				if (A_star_MultipleStack_SC_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 1, decoding_info)) {
					Pre_Best_Goal = Best_Goal_2;
					return true;
				}
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
			else if (pc_i == decoding_info.Constraint_i - 1) {
				Best_Goal_2 = Best_Goal;
				if (A_star_MultipleStack_SC_inner(Sorted_G, Metric_Table, Hard_RX, MRIP_codeword, Pointer, Best_Goal_2, pc_i + 2, decoding_info)) {
					Pre_Best_Goal = Best_Goal_2;
					return true;
				}
				Update_Best_Goal_Procedure(Best_Goal_2, Best_Goal, Stack);
			}
			else {
				for (size_t j(Child_Node.level); j < message_length; ++j) {
					//for subcode SISO Astar using
					if (Metric_Table._matrix[Hard_RX.at(j)][j])
						Child_Node.metric += Metric_Table._matrix[Hard_RX.at(j)][j];
					else Child_Node.Same_Index.push_back(j);
					Child_Node.message_bits.at(j) = Hard_RX.at(j);
				}

				codeword_seq = MRIP_codeword;
				for (size_t index(0); index < Child_Node.Diff_Index.size(); ++index) {
					codeword_seq.at(Child_Node.Diff_Index.at(index)) ^= 1;
					for (size_t j(message_length); j < codeword_length; ++j) {
						codeword_seq.at(j) ^= Sorted_G._matrix[Child_Node.Diff_Index.at(index)][j]; //這邊的G是Sorted_G
					}
				}
				decoding_info.STE += (codeword_length - Child_Node.level);
				++decoding_info.CandidateCodeWord;
				for (size_t j(message_length); j < codeword_length; ++j) {
					if (Metric_Table._matrix[codeword_seq[j]][j]) {
						Child_Node.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
						if (Child_Node.metric >= Best_Goal.metric) break;
					}
					else Child_Node.Same_Index.push_back(j);
				}
				Update_Flag = Update_Best_Goal_Procedure(Child_Node, Best_Goal, Stack);
				if (Update_Flag == TRUE) {
					q_i = (dmin*1.2) - Best_Goal.D_z;
					if (q_i < 0) q_i = 0;
					for (size_t j(0); j < q_i; ++j) {
						Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
					}
					if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
					if (Best_Goal.metric < Sufficient_Condition_thr) {
						Pre_Best_Goal = Best_Goal;
						return true;
					}
					Temp_Sufficient_Condition_thr = 0;
				}
			}
		}
	} while (!Stack.empty());

	Pre_Best_Goal = Best_Goal;
	return false;
}



void A_star_constraint_i_noComparison(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
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
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;
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
				if (Metric_Table._matrix[tmp.message_bits.at(i)][i]) {
					tmp.metric += Metric_Table._matrix[tmp.message_bits.at(i)][i];
				}
				if (Pointer.metric < Best_Goal.metric)
					Stack.push(tmp);
			}
			decoding_info.STE += 2 * (message_length - Pointer.level);
		}
		//Update best goal
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
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
	} while (!Stack.empty());

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


void Subcode_A_star_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static bool flag = true;
	static MATRIX<__int8> G_;
	if (flag) {
		G_.Building_Empty_Matrix(G.Row_number, G.Col_number / 2);
		flag = false;
	}
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
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	vector<double> tmp_Rx = decoding_info.rx_signal_seq;
	Sort_Function(tmp_Rx, Location_Index, decoding_info.rx_signal_seq);
	decoding_info.rx_signal_seq.resize(G.Col_number / 2);
	for (int i = 0; i < G_.Row_number; i++) {
		for (int j = 0; j < G_.Col_number; j++) {
			G_._matrix[i][j] = Sorted_G._matrix[i][j];
		}
	}
	/*
	//test
	static long err_count[128];
	static int count = 0;
	int err_bit = 0;
	count++;
	vector<__int8> Sorted_codeword(codeword_length, 0);
	Sort_Function(decoding_info.code_seq, Location_Index, Sorted_codeword);
	for (int i = 0; i < message_length; i++) {
		if (decoding_info.rx_signal_seq[i] < 0) decoding_info.estimated_codeword[i] = 1;
		else decoding_info.estimated_codeword[i] = 0;
	}
	for (int i = 0; i < message_length; i++) {
		if (Sorted_codeword[i] != decoding_info.estimated_codeword[i]) {
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
	*/
	//A_star_I(G_, decoding_info);
	//A_star_2_stack(G_, decoding_info);
	A_star_Segment(G_, decoding_info);
	decoding_info.STE *= G_.Row_number;
	decoding_info.COM *= G_.Row_number;
	Hard_RX = decoding_info.estimated_codeword;
	
	
	

	//decoding_info.Constraint_i = 4;
	//Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	//message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(message_seq.begin(), message_seq.begin() + message_length);

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
	decoding_info.rx_signal_seq = tmp_Rx;

}

//Brute Force + A*
//only use (8,4)Hamming code
void BF_SPA_BF_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	MATRIX<__int8> H, G_, G__,H_;
	H._matrix = G._matrix_outer;
	H.Col_number = H._matrix.at(0).size();
	H.Row_number = H._matrix.size();
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_.Col_number = G._H.at(0).size();
	H_.Row_number = G._H.size();
	H_._matrix = G._H;
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
		Rx,
		tmp_Rx = decoding_info.rx_signal_seq,
		Rx_parity(G.Col_number - G_.Col_number);
	Rx.assign(decoding_info.rx_signal_seq.begin(), decoding_info.rx_signal_seq.begin() + G_.Col_number);

	for (i = 0; i < G.Row_number / (H.Col_number - H.Row_number); i++) {
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
	//SPA
	for (int i = 0; i < G.Row_number; i++) Rx.at(i) /= (2 / decoding_info.var);
	SPA_subDecoder2(H_,Rx,decoding_info.var,decoding_info.SPA_I);
	for (int i = G.Row_number; i < G_.Col_number; i++) Rx.at(i) /= (2 / decoding_info.var);
	//
	for (; i < G.Col_number / H.Col_number; i++) {
		//assign sub-Block
		start = i * (H.Col_number - H.Row_number);
		for (j = 0; j < H.Col_number - H.Row_number; j++) {
			subRx.at(j) = Rx.at(start + j);
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
	//test
	//output of SISO is extrinsic LLR
	/*
	for (int i = G_.Row_number; i < G_.Col_number; i++) {
		decoding_info.rx_signal_seq[i] = tmp_Rx[i] * (2 / decoding_info.var);
	}
	*/
	//original A*
	A_star_I(G_, decoding_info);
	//RS_Hard_desicion(G_, decoding_info);
	//set parity LLR estimate
	decoding_info.estimated_codeword.resize(G.Col_number);

	for (i = G_.Col_number; i < G.Col_number; i++) {
		if (Rx_parity.at(i - G_.Col_number) < 0) decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	decoding_info.rx_signal_seq = tmp_Rx;

}


void SISO_Path_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number/G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i*parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);
	int Min_Dz = INT_MAX;
	//SC Update flag if flag == 1 Update SC
	int SC_flag = 0;
	//
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
		decoding_info.Sorted_R = sorting_rx_signal_seq;
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;

	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	decoding_info.Counter = 0;
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
			if ((Min_Dz >= Best_Goal.D_z) && Best_Goal.D_z >= 2) {
				SC_flag = (SC_flag + 1) % 10;
				if (SC_flag != 1) {
					continue;
				}
				Min_Dz = Best_Goal.D_z;
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for(int i(0); i < codeword_length; i++) {
					if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
						decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
						Best_Goal.Same_Index.push_back(i);
				}
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) break;
				Temp_Sufficient_Condition_thr = 0;
			}
		}
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

void Subcode24_12_MRIP_A_star(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, G___, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	H_inner._matrix = G._H;
	H_inner.Col_number = G_.Col_number;
	H_inner.Row_number = G_.Col_number - G_.Row_number;
	//for BF
	static vector<vector<double*>> Rx;
	for (int i = 0; i < 5; i++) Rx.push_back(vector<double*>(24, 0));
	Rx.push_back(vector<double*>(8, 0));
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__.Building_Empty_Matrix(12, 24);
		G___.Building_Empty_Matrix(4, 8);
		G___._matrix = {
		{1,0,0,0,1,1,1,0},
		{0,1,0,0,1,1,0,1},
		{0,0,1,0,1,0,1,1},
		{0,0,0,1,0,1,1,1}
		};
		G__._matrix = {
		{1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0},
		{0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,1,1,1,0,1},
		{0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0,0},
		{0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,0},
		{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1},
		{0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,1,1,0,0,1,1,0,0},
		{0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,1,1,0,0,1,1,0},
		{0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,1,1},
		{0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,1,1,1,0,0,0,1,1},
		{0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,1,0,1,1},
		{0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1},
		{0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,1}
		};
		//golay code (24,12)
		for (int i = 0; i < 5; i++) {
			//message
			for (int j = 0; j < 12; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i * 12 + j));
				Rx.at(i).at(12 + j) = &(tmp_Rx.at(64 + i * 12 + j));
			}
		}
		//HM(8,4) code
		for (int j = 0; j < 4; j++) {
			Rx.at(5).at(j) = &(tmp_Rx.at(60 + j));
			Rx.at(5).at(4 + j) = &(tmp_Rx.at(124 + j));
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder (if using same sub_decoder has bug)
	for (int i = 0; i < 5; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}
	pointer_BF_subDecoder2(G___, Rx.at(5), decoding_info.var);

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
	} while (0);
	//test
	/*
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
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Desort_Function(Location_Index, MRIP_codeword, decoding_info.estimated_codeword);
	return;
	*/
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
	//decoding_info.Constraint_i = 4;
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

void SISO_Path_A_star_CBC2_OSC(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i * parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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
		CBC1_codeword_seq(codeword_length, 0),
		CBC2_codeword_seq(codeword_length, 0);

	MATRIX<__int8> Sorted_G(G);
	MATRIX<double> Metric_Table(2, codeword_length);

	NODE_PATH
		Pointer(message_length),
		Best_Goal(message_length),
		Best_Goal_2(message_length), // 用來接收下一層 A* 回傳的結果
		Child_Node(message_length),
		Temp_Node(message_length);

	queue<NODE_PATH> Stack;
	Stack.push(Pointer);

	Best_Goal.metric = FLT_MAX;
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;

	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	decoding_info.Counter = 0;
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
		else if ((Pointer.D_z == decoding_info.Constraint_i && Pointer.level < message_length)) {
			codeword_seq = MRIP_codeword;
			for (size_t index(0); index < Pointer.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G._matrix[Pointer.Diff_Index.at(index)][j];
				}
			}
			//1-bit
			for (int i = Pointer.level; i < message_length; i++) {
				CBC1_codeword_seq = codeword_seq;
				CBC1_codeword_seq.at(i) ^= 1;
				int error_counter = Pointer.D_z+1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					CBC1_codeword_seq.at(j) ^= Sorted_G._matrix[i][j];
				}
				//CBC
				decoding_info.STE += (codeword_length - message_length)+1;
				for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
					if (CBC1_codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
				}
				if (error_counter < decoding_info.Constraint_j) {
					Temp_Node = Pointer;
					Temp_Node.message_bits.at(i) ^= 1;
					for (size_t j(Temp_Node.level); j < codeword_length; ++j) {
						Temp_Node.metric += Metric_Table._matrix[CBC1_codeword_seq.at(j)][j];
						if (Temp_Node.metric > Best_Goal.metric) break;
					}
					if ((Temp_Node.auxiliary()) < Best_Goal.metric) {
						Best_Goal = Temp_Node;
					}
				}
				//2-bits
				for (int j = i + 1; j < message_length; j++) {
					CBC2_codeword_seq = CBC1_codeword_seq;
					CBC2_codeword_seq.at(j) ^= 1;
					int error_counter = Pointer.D_z+2;
					for (size_t k(message_length); k < codeword_length; ++k) {
						CBC2_codeword_seq.at(k) ^= Sorted_G._matrix[j][k];
					}
					//CBC
					decoding_info.STE += (codeword_length - message_length)+2;
					for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
						if (CBC2_codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
					}
					if (error_counter < decoding_info.Constraint_j) {
						Temp_Node = Pointer;
						Temp_Node.message_bits.at(i) ^= 1;
						Temp_Node.message_bits.at(j) ^= 1;
						for (size_t j(Temp_Node.level); j < codeword_length; ++j) {
							Temp_Node.metric += Metric_Table._matrix[CBC2_codeword_seq.at(j)][j];
							if (Temp_Node.metric > Best_Goal.metric) break;
						}
						if ((Temp_Node.auxiliary()) < Best_Goal.metric) {
							Best_Goal = Temp_Node;
						}
					}
				}
			}
		}
		//Update best goal
		codeword_seq = MRIP_codeword;
		++decoding_info.Counter;
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
		/*
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}*/
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
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

void A_star_constraint_i_Min_Heap(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
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

	auto comp = [&](const NODE_PATH& a, const NODE_PATH& b) {
		decoding_info.COM++;
		return a.metric > b.metric;
	};
	priority_queue<NODE_PATH, vector<NODE_PATH>,decltype(comp)> Stack(comp);

	Best_Goal.metric = FLT_MAX;
	Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);

	for (size_t i(0); i < codeword_length; ++i)
		if (Metric_Table._matrix[0][i] != 0)  Hard_RX.at(i) = 1;
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Pointer.message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Stack.push(Pointer);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;
	do {
		Pointer = Stack.top();
		Stack.pop();

		if (Pointer.D_z < decoding_info.Constraint_i && Pointer.level < message_length) {
			for (int i = Pointer.level; i < message_length; i++) {
				NODE_PATH tmp = Pointer;
				tmp.message_bits.at(i) ^= 1;
				++tmp.D_z;
				tmp.Diff_Index.push_back(i);
				tmp.level = i + 1;
				if (Metric_Table._matrix[tmp.message_bits.at(i)][i]) {
					tmp.metric += Metric_Table._matrix[tmp.message_bits.at(i)][i];
				}
				if (Pointer.metric < Best_Goal.metric)
					Stack.push(tmp);
			}
			decoding_info.STE += 2 * (message_length - Pointer.level);
		}
		//Update best goal
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
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
	} while (!Stack.empty() && Stack.top().metric < Best_Goal.metric);

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

void SISO_Path_A_star_Min_Heap(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i * parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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

	auto comp = [&](const NODE_PATH& a, const NODE_PATH& b) {
		decoding_info.COM++;
		return a.metric > b.metric;
	};
	priority_queue<NODE_PATH, vector<NODE_PATH>, decltype(comp)> Stack(comp);
	

	Best_Goal.metric = FLT_MAX;
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);
	int Min_Dz = INT_MAX;
	//SC Update flag if flag == 1 Update SC
	int SC_flag = 0;
	//
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
		decoding_info.Sorted_R = sorting_rx_signal_seq;
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;

	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Pointer.message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Stack.push(Pointer);;
	decoding_info.Counter = 0;
	do {
		Pointer = Stack.top();
		Stack.pop();
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
			/*
			if ((Min_Dz >= Best_Goal.D_z) && Best_Goal.D_z >= 2) {
				SC_flag = (SC_flag + 1) % 10;
				if (SC_flag != 1) {
					continue;
				}
				Min_Dz = Best_Goal.D_z;
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for (int i(0); i < codeword_length; i++) {
					if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
						decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
						Best_Goal.Same_Index.push_back(i);
				}
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) break;
				Temp_Sufficient_Condition_thr = 0;
			}*/
		}
	} while (!Stack.empty() && Stack.top().metric < Best_Goal.metric);

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

void SISO_Path_A_star_Adaptive(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i * parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);
	int Min_Dz = INT_MAX;
	bool Break_flag = false;
	//
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
		decoding_info.Sorted_R = sorting_rx_signal_seq;
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;

	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Stack.front().message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	decoding_info.Counter = 0;
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
			/*
			if ((Min_Dz >= Best_Goal.D_z) && Best_Goal.D_z >= 2) {
				Min_Dz = Best_Goal.D_z;
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for (int i(0); i < codeword_length; i++) {
					if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
						decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
						Best_Goal.Same_Index.push_back(i);
				}
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) break;
				Temp_Sufficient_Condition_thr = 0;
			}
			*/
		}
	} while (!Stack.empty() && Stack.front().D_z <= 2);
	//
	if (Best_Goal.metric < 1.4 * OSC_metric_thr) Break_flag = true;
	/*
	Min_Dz = Best_Goal.D_z;
	q_i = (dmin*1.4) - Best_Goal.D_z;
	if (q_i < 0) q_i = 0;
	for (int i(0); i < codeword_length; i++) {
		if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
			decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
			Best_Goal.Same_Index.push_back(i);
	}
	for (size_t j(0); j < q_i; ++j) {
		Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
	}
	if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
	if (Best_Goal.metric < Sufficient_Condition_thr) Break_flag = true;
	Temp_Sufficient_Condition_thr = 0;
	*/
	//
	while(!Stack.empty() && !Break_flag){
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
			/*
			if ((Min_Dz >= Best_Goal.D_z) && Best_Goal.D_z >= 2) {
				Min_Dz = Best_Goal.D_z;
				q_i = (dmin*1.2) - Best_Goal.D_z;
				if (q_i < 0) q_i = 0;
				for (int i(0); i < codeword_length; i++) {
					if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
						decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
						Best_Goal.Same_Index.push_back(i);
				}
				for (size_t j(0); j < q_i; ++j) {
					Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
				}
				if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr)Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
				if (Best_Goal.metric < Sufficient_Condition_thr) break;
				Temp_Sufficient_Condition_thr = 0;
			}
			*/
		}
	} 
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

void SISO_Path_A_star_Adaptive_Min_Heap(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i * parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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

	bool Break_flag = false;

	auto comp = [&](const NODE_PATH& a, const NODE_PATH& b) {
		decoding_info.COM++;
		return a.metric > b.metric;
	};
	priority_queue<NODE_PATH, vector<NODE_PATH>,std::reference_wrapper<decltype(comp)>> Stack(std::ref(comp)),Next_Stack(std::ref(comp));


	Best_Goal.metric = FLT_MAX;
	//
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G, Location_Index);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table);
		decoding_info.Sorted_R = sorting_rx_signal_seq;
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;

	Sort_Function(decoding_info.estimated_codeword, Location_Index, Hard_RX);
	message_seq.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G, message_seq, MRIP_codeword);
	Pointer.message_bits.assign(Hard_RX.begin(), Hard_RX.begin() + message_length);
	Stack.push(Pointer);;
	decoding_info.Counter = 0;
	do {
		Pointer = Stack.top();
		Stack.pop();
		if (Pointer.metric >= Best_Goal.metric) continue;
		if (Pointer.D_z > 2) {
			Next_Stack.push(Pointer);
			continue;
		}
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
	} while (!Stack.empty() );
	//
	Stack = Next_Stack;
	if (Best_Goal.metric < 1.4 * OSC_metric_thr) Break_flag = true;
	//
	while (!Stack.empty() && !Break_flag) {
		Pointer = Stack.top();
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
		++decoding_info.Counter;
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
		//Control level
		int error_counter = Child_Node.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer.metric += Metric_Table._matrix[codeword_seq.at(j)][j];
				Pointer.D_z++;
			}
			if (Pointer.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
	}
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


void SISO_Path_A_star_MultiBase_Adaptive(MATRIX<__int8>& G, DECODING_INFO& decoding_info) {
	static MATRIX<__int8>  G_, G__, H_inner;
	G_._matrix = G._matrix_inner;
	G_.Col_number = G._matrix_inner.at(0).size();
	G_.Row_number = G._matrix_inner.size();
	//for BF
	static vector<vector<double*>> Rx;
	static vector<double> tmp_Rx;
	tmp_Rx = decoding_info.rx_signal_seq;
	static bool flag = false;
	if (!flag) {
		G__._matrix = G.subG;
		G__.Col_number = G.subG.at(0).size();
		G__.Row_number = G.subG.size();
		for (int i = 0; i < G.Col_number / G__.Col_number; i++)
			Rx.push_back(vector<double*>(G__.Col_number, 0));
		//reference each block
		int parity_start = G.Col_number*G__.Row_number / G__.Col_number;
		int parity_len = G__.Col_number - G__.Row_number;
		for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
			//message
			for (int j = 0; j < G__.Row_number; j++) {
				Rx.at(i).at(j) = &(tmp_Rx.at(i*G__.Row_number + j));
			}
			//parity

			for (int j = 0; j < parity_len; j++) {
				Rx.at(i).at(G__.Row_number + j) =
					&(tmp_Rx.at(parity_start + i * parity_len + j));
			}
		}
		flag = true;
	}
	//for each sub-Matrix use a BF_subDecoder 
	for (int i = 0; i < G.Col_number / G__.Col_number; i++) {
		//Brute-Force SISO
		pointer_BF_subDecoder(G__, Rx.at(i), decoding_info.var);
	}

	for (int i = 0; i < G.Col_number; i++) {
		if (tmp_Rx.at(i) < 0)
			decoding_info.estimated_codeword.at(i) = 1;
		else  decoding_info.estimated_codeword.at(i) = 0;
	}
	//using queue without Comparison
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
	queue<NODE_PATH> Stack_Base1;
	queue<NODE_PATH> Stack_Base2;
	Stack_Base1.push(Pointer_Base1);
	Stack_Base2.push(Pointer_Base2);

	Best_Goal.metric = FLT_MAX;
	decoding_info.Counter = 0;

	Best_Goal.metric = FLT_MAX;
	double Sufficient_Condition_thr(0), Temp_Sufficient_Condition_thr(0);
	int q_i(0);
	int Min_Dz = INT_MAX;
	bool Break_flag = false;
	//
	//Pre_Procedure(decoding_info.rx_signal_seq, G, Sorted_G, Location_Index, Metric_Table);
	do {//Base1
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		Determine_Permutation(tmp_Rx, G, Sorted_G_Base1, Location_Index_Base1);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index_Base1, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table_Base1);
		decoding_info.Sorted_R = sorting_rx_signal_seq;
	} while (0);
	do {//Base2
		vector<double> sorting_rx_signal_seq(tmp_Rx);
		//use channel output
		Determine_Permutation(decoding_info.rx_signal_seq, G, Sorted_G_Base2, Location_Index_Base2);
		Sort_Function(decoding_info.rx_signal_seq, Location_Index_Base2, sorting_rx_signal_seq);
		Build_Metric_Table(sorting_rx_signal_seq, Metric_Table_Base2);
		decoding_info.Sorted_R_Base2 = sorting_rx_signal_seq;
	} while (0);
	// threshold
	double OSC_metric_thr(0);
	for (size_t i(0); i < codeword_length; ++i) {
		// for OSC threshold
		OSC_metric_thr += abs(decoding_info.rx_signal_seq.at(i));
	}
	OSC_metric_thr = decoding_info.OSC_Alpha*OSC_metric_thr;
	//Base1
	Sort_Function(decoding_info.estimated_codeword, Location_Index_Base1, Hard_RX_Base1);
	message_seq_Base1.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, message_seq_Base1, MRIP_codeword_Base1);
	Stack_Base1.front().message_bits.assign(Hard_RX_Base1.begin(), Hard_RX_Base1.begin() + message_length);
	//Base2
	Sort_Function(decoding_info.estimated_codeword, Location_Index_Base2, Hard_RX_Base2);
	message_seq_Base2.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, message_seq_Base2, MRIP_codeword_Base2);
	Stack_Base2.front().message_bits.assign(Hard_RX_Base2.begin(), Hard_RX_Base2.begin() + message_length);
	decoding_info.Counter = 0;
	do {
		if(!Stack_Base1.empty() && Stack_Base1.front().D_z <= 2) {
			Pointer_Base1 = Stack_Base1.front();
			Stack_Base1.pop();
			if (Pointer_Base1.metric >= Best_Goal.metric) continue;
			//flipping node and push them into queue
			if (Pointer_Base1.D_z < decoding_info.Constraint_i && Pointer_Base1.level < message_length) {
				for (int i = Pointer_Base1.level; i < message_length; i++) {
					NODE_PATH tmp = Pointer_Base1;
					tmp.message_bits.at(i) ^= 1;
					++tmp.D_z;
					tmp.Diff_Index.push_back(i);
					tmp.level = i + 1;
					tmp.metric += Metric_Table_Base1._matrix[tmp.message_bits.at(i)][i];
					if (Pointer_Base1.metric < Best_Goal.metric)
						Stack_Base1.push(tmp);
				}
				decoding_info.STE += 2 * (message_length - Pointer_Base1.level);
			}
			//Update best goal
			codeword_seq = MRIP_codeword_Base1;
			++decoding_info.Counter;
			for (int i = Pointer_Base1.level; i < message_length; i++) {
				Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(i)][i];
			}
			if (Pointer_Base1.metric >= Best_Goal.metric) continue;

			for (size_t index(0); index < Pointer_Base1.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer_Base1.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Pointer_Base1.Diff_Index.at(index)][j];
				}
			}
			//Control level
			int error_counter = Child_Node_Base1.D_z;
			for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
				if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
			}
			if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
				continue;
			}
			//
			for (size_t j(message_length); j < codeword_length; ++j) {
				if (Metric_Table_Base1._matrix[codeword_seq.at(j)][j] > 0) {
					Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
					Pointer_Base1.D_z++;
				}
				if (Pointer_Base1.metric > Best_Goal.metric) break;
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//Update
			if ((Pointer_Base1.auxiliary()) < Best_Goal.metric) {
				Best_Goal = Pointer_Base1;
				Best_Goal.base = Sorted_Base_1;
				if (Best_Goal.metric < OSC_metric_thr) break;
			}
		}
		//Base2
		if (!Stack_Base2.empty() && Stack_Base2.front().D_z <= 2) {
			Pointer_Base2 = Stack_Base2.front();
			Stack_Base2.pop();
			if (Pointer_Base2.metric >= Best_Goal.metric) continue;
			//flipping node and push them into queue
			if (Pointer_Base2.D_z < decoding_info.Constraint_i && Pointer_Base2.level < message_length) {
				for (int i = Pointer_Base2.level; i < message_length; i++) {
					NODE_PATH tmp = Pointer_Base2;
					tmp.message_bits.at(i) ^= 1;
					++tmp.D_z;
					tmp.Diff_Index.push_back(i);
					tmp.level = i + 1;
					tmp.metric += Metric_Table_Base2._matrix[tmp.message_bits.at(i)][i];
					if (Pointer_Base2.metric < Best_Goal.metric)
						Stack_Base2.push(tmp);
				}
				decoding_info.STE += 2 * (message_length - Pointer_Base2.level);
			}
			//Update best goal
			codeword_seq = MRIP_codeword_Base2;
			++decoding_info.Counter;
			for (int i = Pointer_Base2.level; i < message_length; i++) {
				Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(i)][i];
			}
			if (Pointer_Base2.metric >= Best_Goal.metric) continue;

			for (size_t index(0); index < Pointer_Base2.Diff_Index.size(); ++index) {
				codeword_seq.at(Pointer_Base2.Diff_Index.at(index)) ^= 1;
				for (size_t j(message_length); j < codeword_length; ++j) {
					codeword_seq.at(j) ^= Sorted_G_Base2._matrix[Pointer_Base2.Diff_Index.at(index)][j];
				}
			}
			//Control level
			int error_counter = Child_Node_Base2.D_z;
			for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
				if (codeword_seq.at(j) != Hard_RX_Base2.at(j)) ++error_counter;
			}
			if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
				continue;
			}
			//
			for (size_t j(message_length); j < codeword_length; ++j) {
				if (Metric_Table_Base2._matrix[codeword_seq.at(j)][j] > 0) {
					Pointer_Base2.metric += Metric_Table_Base2._matrix[codeword_seq.at(j)][j];
					Pointer_Base2.D_z++;
				}
				if (Pointer_Base2.metric > Best_Goal.metric) break;
			}
			//
			decoding_info.STE += (codeword_length - message_length);
			decoding_info.CandidateCodeWord++;
			//Update
			if ((Pointer_Base2.auxiliary()) < Best_Goal.metric) {
				Best_Goal = Pointer_Base2;
				Best_Goal.base = Sorted_Base_2;
				if (Best_Goal.metric < OSC_metric_thr) break;
			}
		}
	} while ((!Stack_Base1.empty() && Stack_Base1.front().D_z <= 2) ||
		     (!Stack_Base2.empty() && Stack_Base2.front().D_z <= 2));
	//
	if (Best_Goal.metric < 1.4 * OSC_metric_thr) Break_flag = true;
	/*if ((Min_Dz >= Best_Goal.D_z)) {
		Min_Dz = Best_Goal.D_z;
		q_i = (dmin*1.4) - Best_Goal.D_z;
		if (q_i < 0) q_i = 0;
		for (int i(0); i < codeword_length; i++) {
			if (decoding_info.rx_signal_seq.at(i) < 0 && !codeword_seq.at(i) ||
				decoding_info.rx_signal_seq.at(i) > 0 && codeword_seq.at(i))
				Best_Goal.Same_Index.push_back(i);
		}
		if (Best_Goal.base == Sorted_Base_1) {
			for (size_t j(0); j < q_i; ++j) {
				Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
			}
		}
		else {
			for (size_t j(0); j < q_i; ++j) {
				Temp_Sufficient_Condition_thr += abs(decoding_info.Sorted_R_Base2.at(Best_Goal.Same_Index.at(Best_Goal.Same_Index.size() - j - 1)));
			}
		}
		if (Temp_Sufficient_Condition_thr > Sufficient_Condition_thr) Sufficient_Condition_thr = Temp_Sufficient_Condition_thr;
		if (Best_Goal.metric < Sufficient_Condition_thr) flag = true;
		Temp_Sufficient_Condition_thr = 0;
	}*/
	//
	while (!Stack_Base1.empty() && !Break_flag) {
		Pointer_Base1 = Stack_Base1.front();
		Stack_Base1.pop();
		if (Pointer_Base1.metric >= Best_Goal.metric) continue;
		//flipping node and push them into queue
		if (Pointer_Base1.D_z < decoding_info.Constraint_i && Pointer_Base1.level < message_length) {
			for (int i = Pointer_Base1.level; i < message_length; i++) {
				NODE_PATH tmp = Pointer_Base1;
				tmp.message_bits.at(i) ^= 1;
				++tmp.D_z;
				tmp.Diff_Index.push_back(i);
				tmp.level = i + 1;
				tmp.metric += Metric_Table_Base1._matrix[tmp.message_bits.at(i)][i];
				if (Pointer_Base1.metric < Best_Goal.metric)
					Stack_Base1.push(tmp);
			}
			decoding_info.STE += 2 * (message_length - Pointer_Base1.level);
		}
		//Update best goal
		codeword_seq = MRIP_codeword_Base1;
		++decoding_info.Counter;
		for (int i = Pointer_Base1.level; i < message_length; i++) {
			Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(i)][i];
		}
		if (Pointer_Base1.metric >= Best_Goal.metric) continue;

		for (size_t index(0); index < Pointer_Base1.Diff_Index.size(); ++index) {
			codeword_seq.at(Pointer_Base1.Diff_Index.at(index)) ^= 1;
			for (size_t j(message_length); j < codeword_length; ++j) {
				codeword_seq.at(j) ^= Sorted_G_Base1._matrix[Pointer_Base1.Diff_Index.at(index)][j];
			}
		}
		//Control level
		int error_counter = Child_Node_Base1.D_z;
		for (__int16 j(message_length); j < decoding_info.Control_Level; ++j) {
			if (codeword_seq.at(j) != Hard_RX_Base1.at(j)) ++error_counter;
		}
		if (error_counter > decoding_info.Constraint_j && decoding_info.Counter > 2) {
			continue;
		}
		//
		for (size_t j(message_length); j < codeword_length; ++j) {
			if (Metric_Table_Base1._matrix[codeword_seq.at(j)][j] > 0) {
				Pointer_Base1.metric += Metric_Table_Base1._matrix[codeword_seq.at(j)][j];
				Pointer_Base1.D_z++;
			}
			if (Pointer_Base1.metric > Best_Goal.metric) break;
		}
		//
		decoding_info.STE += (codeword_length - message_length);
		decoding_info.CandidateCodeWord++;
		//Update
		if ((Pointer_Base1.auxiliary()) < Best_Goal.metric) {
			Best_Goal = Pointer_Base1;
			Best_Goal.base = Sorted_Base_1;
			if (Best_Goal.metric < OSC_metric_thr) break;
		}
	}
	//
	if (Best_Goal.base == Sorted_Base_1) {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base1, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base1, codeword_seq, decoding_info.estimated_codeword);
	}
	else {
		Systematic_Linear_Block_Code_Encoder(Sorted_G_Base2, Best_Goal.message_bits, codeword_seq);
		Desort_Function(Location_Index_Base2, codeword_seq, decoding_info.estimated_codeword);
	}

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

