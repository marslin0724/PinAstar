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
		SPA_subDecoder(H, subRx, decoding_info.var,decoding_info.SPA_I);
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
		SPA_subDecoder(H, subRx, decoding_info.var,decoding_info.SPA_I);
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
	/*
	for (i = 0; i < G.Col_number; i++) {
		decoding_info.rx_signal_seq.at(i) /= (2 / decoding_info.var);
	}*/
	//original A*
	A_star_I(G, decoding_info);
	decoding_info.rx_signal_seq = tmp_Rx;
}