#include "PinAstar/TURBO_OPER.h"

void Turbo_G_Generator(MATRIX<__int8> &G, int row, int col) { // g1 = [1 1 1], g2 = [1 0 1], RSC encoder 
	
	G.Building_Empty_Matrix(row + 2, col + 4);
	

}

void Turbo_LDPC_encoder(MATRIX<__int8> &G_, MATRIX<__int8> &G, vector<__int8> &message_seq, vector<__int8> &output_codeword_seq) {
	vector<__int8> mes = message_seq, codeword_2;
	//Systematic_Linear_Block_Code_Encoder(G_, mes, output_codeword_seq);
	int temp;
	output_codeword_seq.resize(G_.Col_number);
	for (int i = 0; i < G_.Col_number; ++i) {
		temp = 0;
		for (int j = 0; j < G_.Row_number; j++) {
			temp ^= (mes.at(j)&G_._matrix[j][i]);
			//temp += Message_Seq.at(j)*LinearBlockCode.G._matrix[j][i];
		}
		output_codeword_seq.at(i) = temp;
	}
	Block_Interleaver(mes);
	//Systematic_Linear_Block_Code_Encoder(G, mes, codeword_2);
	codeword_2.resize(G.Col_number);
	for (int i = 0; i < G.Col_number; ++i) {
		temp = 0;
		for (int j = 0; j < G.Row_number; j++) {
			temp ^= (mes.at(j)&G._matrix[j][i]);
			//temp += Message_Seq.at(j)*LinearBlockCode.G._matrix[j][i];
		}
		codeword_2.at(i) = temp;
	}
	
	output_codeword_seq.insert(output_codeword_seq.end(), codeword_2.begin() + message_seq.size(), codeword_2.begin() + codeword_2.size());
}

void Turbo_LDPC_encoder_Punctured_Version(MATRIX<__int8> &G, vector<__int8> &message_seq, vector<__int8> &output_codeword_seq) {
	vector<__int8> mes = message_seq, codeword_2;
	Systematic_Linear_Block_Code_Encoder(G, mes, output_codeword_seq);
	output_codeword_seq.erase(output_codeword_seq.begin() + G.Col_number*0.75, output_codeword_seq.end());
	output_codeword_seq.resize(G.Col_number*0.75);
	Block_Interleaver(mes);
	Systematic_Linear_Block_Code_Encoder(G, mes, codeword_2);
	output_codeword_seq.insert(output_codeword_seq.end(), codeword_2.begin() + message_seq.size(), codeword_2.begin() + G.Col_number*0.75);
}

void Turbo_LDPC_decoder(MATRIX<__int8> &H_, 
	MATRIX<__int8> &H, MATRIX<__int8> &G, 
	DECODING_INFO &decoding_info, 
	vector<int> &permutation_seq_,
	vector<int> &permutation_seq)
{
	//cout << "I am here" << decoding_info.var << endl << endl;
	// Intialization
	int Col = H.Col_number,
		Row = H.Row_number,
		Col_ = H_.Col_number,
		Row_ = H_.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		negative_amount,
		symdrone,
		codeword_length;
	long double temp;
	double Average_LLR, previous, Best_Hard_Result_Metric;
	Best_Hard_Result_Metric = DBL_MAX;

	vector<__int8> Hard_Result(Col, 0),
		Current_Hard_Result(Col, 2),
		Bad_index(Message_Length, 0), 
		Previous_Hard_Result(Col, 0), 
		Best_Hard_Result(Col, 0);

	vector<double>
		Channel_LLR_1(Col_, 0),
		Channel_LLR_2(Col, 0),
		Channel_LLR_Message_1(Message_Length, 0),
		Channel_LLR_Message_2(Message_Length, 0),
		Channel_LLR_Message_Temp(Message_Length, 0),
		Channel_LLR_CodeWord_Temp(Col, 0),
		Channel_LLR_Parity_1(Row_, 0),
		Channel_LLR_Parity_2(Row, 0),
		Channel_LLR_Parity_Temp(Row, 0),
		Extrinsic_Mes_Passing_1to2(Col_, 0),
		Extrinsic_Mes_Passing_2to1(Col, 0),
		Extrinsic_Mes_Passing_1to2_MesPart(Message_Length, 0),
		Extrinsic_Mes_Passing_2to1_MesPart(Message_Length, 0),
		ZeroSequence(Col, 0),
		ZeroSequence_(Col_, 0),
		Accumulated_LLR(Col, 0);

	MATRIX<long double>
		Message_info_1,
		Extrinsis_Mes_1,
		Message_info_2,
		Extrinsis_Mes_2,
		ZeroMatrix;

	Extrinsis_Mes_1.Building_Empty_Matrix(Row_, Col_);
	Message_info_1.Building_Empty_Matrix(Row_, Col_);
	Extrinsis_Mes_2.Building_Empty_Matrix(Row, Col);
	Message_info_2.Building_Empty_Matrix(Row, Col);
	ZeroMatrix.Building_Empty_Matrix(Row, Col);

	LDPC_FUNC LDPC;
	codeword_length = Col + Row_;

	vector<double> Accumulated_Residual(codeword_length, 0), Previous_Metric(codeword_length, 0);
	
	double Original_LLR = 0;
	vector<__int8> Original_Hard_Result(codeword_length, 0), Now_Hard_Result(codeword_length, 0), Good_Hard_Result(codeword_length, 0);
	for (int i = 0; i < codeword_length; ++i) {
		if (decoding_info.rx_signal_seq.at(i) > 0) Original_Hard_Result.at(i) = 0;
		else  Original_Hard_Result.at(i) = 1;

		Original_LLR += abs(decoding_info.rx_signal_seq.at(i));
	}

	//double 

	//cout << Col << "," << Row << "," << Col_ << "," << Row_ << "," << endl;
	for (int i = 0; i < Col_; ++i) {
		Channel_LLR_1.at(i) = 2 * decoding_info.rx_signal_seq.at(i) / decoding_info.var;
	}
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_Message_1.at(i) = Channel_LLR_1.at(i);
		Channel_LLR_Message_2.at(i) = Channel_LLR_1.at(i);
	}
	Block_Interleaver(Channel_LLR_Message_2);
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_2.at(i) = Channel_LLR_Message_2.at(i);
	}
	for (int i = Message_Length; i < Col; ++i) {
		Channel_LLR_2.at(i) = 2 * decoding_info.rx_signal_seq.at(i + Row_) / decoding_info.var;
		Channel_LLR_Parity_2.at(i - Message_Length) = Channel_LLR_2.at(i);
	}
	for (int i = Message_Length; i < Col_; ++i) {
		Channel_LLR_Parity_1.at(i - Message_Length) = Channel_LLR_1.at(i);
	}
	//cout << "A";

	if (Diff_Candidate_Test) {

		decoding_info.Candidate_Metric.clear();
		
		for (int i = 0; i < Message_Length; i++) {
			decoding_info.rx_signal_seq.at(i) = Channel_LLR_Message_2.at(i)*decoding_info.var / 2;
		}
		for (int i = Message_Length; i < Col; i++) {
			decoding_info.rx_signal_seq.at(i) = decoding_info.rx_signal_seq.at(i + (Message_Length / 2));
		}
		decoding_info.rx_signal_seq.resize(Col);
	}
	/*
	// Error Test
	vector<__int8> CodeWord_Temp(192, 0);
	int Record = 0;
	for (int index = 0; index < Col; index++) {
		if (Channel_LLR_2.at(index) > 0) CodeWord_Temp.at(index) = 0;
		else CodeWord_Temp.at(index) = 1;
	}
	++Record;
	*/
	for (int i = 0; i < Col; ++i) {
		if (Channel_LLR_2.at(i) < 0) Previous_Hard_Result.at(i) = 1;
	}
	
	double AveLLR = 0, Final_LLR = 0, Previous_AveLLR = 0, Latest_AveLLR = 0;
	decoding_info.EarlyTerminate = FALSE;
	//cout << endl;
	if (Renew_One_Bit) {

		// Decoder 1
		LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
		LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);
		for (int PreIteration = 0; PreIteration < 5; PreIteration++) {

			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						temp = Channel_LLR_1.at(col);
						for (int vertical = 0; vertical < Row_; vertical++) {
							if (H_._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_1._matrix[vertical][col];
						}
						if (Modified_LLR_v2c) {       /// Modified LLR Osc version
							if (temp*Message_info_1._matrix[row][col] < 0) Message_info_1._matrix[row][col] += temp;
							else Message_info_1._matrix[row][col] = temp;
						}
						else {                    ///  Conventional version
							Message_info_1._matrix[row][col] = temp;
						}
						Message_info_1._matrix[row][col] *= Modified_LLR_v2c_factor;
					}
				}
			}
			Extrinsic_Mes_Passing_1to2.assign(ZeroSequence_.begin(), ZeroSequence_.end());

			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						if (SumProduct) { /// Sum Product
							temp = 1;
						}
						else {            /// Min Sum
							negative_amount = 0;
							temp = DBL_MAX;
						}
						for (int horizontal = 0; horizontal < Col_; horizontal++) {
							if (H_._matrix[row][horizontal] == 1 && horizontal != col) {

								if (SumProduct) { /// Sum Product
									temp *= tanh(Message_info_1._matrix[row][horizontal] / 2);
								}
								else {            /// Min Sum
									if (Message_info_1._matrix[row][horizontal] < 0) ++negative_amount;
									if (abs(Message_info_1._matrix[row][horizontal]) < temp) temp = abs(Message_info_1._matrix[row][horizontal]);
								}
							}
						}
						if (Modified_LLR_c2v) previous = Extrinsis_Mes_1._matrix[row][col];
						if (SumProduct) { ///Sum Product
							temp = log((1 + temp) / (1 - temp));
							if (abs(temp) < LLR_MAX) Extrinsis_Mes_1._matrix[row][col] = temp;
							else Extrinsis_Mes_1._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
						}
						else {           ///Min Sum
							Extrinsis_Mes_1._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
						}
						if (Modified_LLR_c2v && previous*Extrinsis_Mes_1._matrix[row][col] < 0) {
							Extrinsis_Mes_1._matrix[row][col] += previous;
							Extrinsis_Mes_1._matrix[row][col] *= Modified_LLR_c2v_factor;
						}

						Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]);
					}
				}
			}
			// Residual Test
			double Max_Residual = 0;
			int index = -1;
			for (int i = 0; i < Col; ++i) {
				if (abs(Extrinsic_Mes_Passing_1to2.at(i)) > Max_Residual && Extrinsic_Mes_Passing_1to2.at(i)*Channel_LLR_1.at(i) < 0) {
					Max_Residual = abs(Extrinsic_Mes_Passing_1to2.at(i));
					index = i;
				}
			}
			if (index != -1) {
				Channel_LLR_1.at(index) += Extrinsic_Mes_Passing_1to2.at(index);
				Channel_LLR_2.at(index) += Extrinsic_Mes_Passing_1to2.at(index);
			}
			Extrinsis_Mes_1 = ZeroMatrix;
			/*
			if (index != -1) {
				for (int vertical = 0; vertical < Row_; ++vertical) {
					if (H_._matrix[vertical][index]) {
						for (int j = 0; j < Col; ++j) {
							if (H_._matrix[vertical][j]) {
								Extrinsis_Mes_1._matrix[vertical][j] = 0;
							}
						}
					}
				
				}
			}*/

			// Decoder 2

			for (int row = 0; row < Row; ++row) {
				for (int col = 0; col < Col; ++col) {
					if (H._matrix[row][col] == 1) {
						temp = Channel_LLR_2.at(col);
						for (int vertical = 0; vertical < Row; vertical++) {
							if (H._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_2._matrix[vertical][col];
						}
						if (Modified_LLR_v2c) { /// Modified LLR Osc version
							if (temp*Message_info_2._matrix[row][col] < 0) Message_info_2._matrix[row][col] += temp;
							else Message_info_2._matrix[row][col] = temp;
						}
						else {              ///  Conventional version
							Message_info_2._matrix[row][col] = temp;
						}
						Message_info_2._matrix[row][col] *= Modified_LLR_v2c_factor;
					}
				}
			}

			// check to variable
			Extrinsic_Mes_Passing_2to1.assign(ZeroSequence.begin(), ZeroSequence.end());

			for (int row = 0; row < Row; ++row) {
				for (int col = 0; col < Col; ++col) {
					if (H._matrix[row][col] == 1) {
						if (SumProduct) { /// Sum Product
							temp = 1;
						}
						else {            /// Min Sum
							negative_amount = 0;
							temp = DBL_MAX;
						}
						for (int horizontal = 0; horizontal < Col; horizontal++) {
							if (H._matrix[row][horizontal] == 1 && horizontal != col) {
								if (SumProduct) { /// Sum Product
									temp *= tanh(Message_info_2._matrix[row][horizontal] / 2);
								}
								else {            /// Min Sum
									if (Message_info_2._matrix[row][horizontal] < 0) ++negative_amount;
									if (abs(Message_info_2._matrix[row][horizontal]) < temp) temp = abs(Message_info_2._matrix[row][horizontal]);
								}
							}
						}
						if (Modified_LLR_c2v) previous = Extrinsis_Mes_2._matrix[row][col];
						if (SumProduct) { ///Sum Product
							temp = log((1 + temp) / (1 - temp));
							if (abs(temp) < LLR_MAX) Extrinsis_Mes_2._matrix[row][col] = temp;
							else Extrinsis_Mes_2._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
						}
						else {            ///Min Sum
							Extrinsis_Mes_2._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
						}
						if (Modified_LLR_c2v && Extrinsis_Mes_2._matrix[row][col] * previous < 0) {
							Extrinsis_Mes_2._matrix[row][col] += previous;
							Extrinsis_Mes_2._matrix[row][col] *= Modified_LLR_c2v_factor;
						}

						//if (Extrinsis_Mes_2._matrix[row][col] * Message_info_2._matrix[row][col] < 0) Extrinsis_Mes_2._matrix[row][col] *= (Test_parameter);

						Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]);
						//else Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]*0.7);
						//cout << R << ",";
					}
				}
			}

			// Residual Test
		    Max_Residual = 0;
		    index = -1;
			for (int i = 0; i < Col; ++i) {
				if (abs(Extrinsic_Mes_Passing_2to1.at(i)) > Max_Residual && Extrinsic_Mes_Passing_2to1.at(i)*Channel_LLR_2.at(i) < 0) {
					Max_Residual = abs(Extrinsic_Mes_Passing_2to1.at(i));
					index = i;
				}
			}
			if (index != -1) {
				Channel_LLR_2.at(index) += Extrinsic_Mes_Passing_2to1.at(index);
				Channel_LLR_1.at(index) += Extrinsic_Mes_Passing_2to1.at(index);
			}
			Extrinsis_Mes_2 = ZeroMatrix;
			/*
			if (index != -1) {
				for (int vertical = 0; vertical < Row_; ++vertical) {
					if (H._matrix[vertical][index]) {
						for (int j = 0; j < Col; ++j) {
							if (H._matrix[vertical][j]) {
								Extrinsis_Mes_2._matrix[vertical][j] = 0;
							}
						}
					}

				}
			}*/

		}
		LDPC.Sort_Sequence_Forward(Channel_LLR_1, permutation_seq_);
		LDPC.Sort_Sequence_Forward(Channel_LLR_2, permutation_seq);

		
		for (int i = 0; i < Message_Length; ++i) {
			Channel_LLR_Message_1.at(i) = Channel_LLR_1.at(i);
		}
		for (int i = 0; i < Row_; ++i) {
			Channel_LLR_Parity_1.at(i) = Channel_LLR_1.at(i + Message_Length);
		}
		for (int i = 0; i < Message_Length; ++i) {
			Channel_LLR_Message_2.at(i) = Channel_LLR_2.at(i);
		}
		for (int i = 0; i < Row_; ++i) {
			Channel_LLR_Parity_2.at(i) = Channel_LLR_2.at(i + Message_Length);
		}
	}


	// Iteration Start
	while (1) {
		// Decoder 1
		// variable to check

		LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
		
		//
		for (int local_iter = 0; local_iter < Local_Iteration; local_iter++) {

			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						temp = Channel_LLR_1.at(col);
						for (int vertical = 0; vertical < Row_; vertical++) {
							if (H_._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_1._matrix[vertical][col];
						}
						if (Modified_LLR_v2c) {       /// Modified LLR Osc version
							if (temp*Message_info_1._matrix[row][col] < 0) Message_info_1._matrix[row][col] += temp;
							else Message_info_1._matrix[row][col] = temp;
						}
						else {                    ///  Conventional version
							Message_info_1._matrix[row][col] = temp;
						}
						Message_info_1._matrix[row][col] *= Modified_LLR_v2c_factor;
					}
				}
			}

			// check to variable
			Extrinsic_Mes_Passing_1to2.assign(ZeroSequence_.begin(), ZeroSequence_.end());
			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						if (SumProduct) { /// Sum Product
							temp = 1;
						}
						else {            /// Min Sum
							negative_amount = 0;
							temp = DBL_MAX;
						}
						for (int horizontal = 0; horizontal < Col_; horizontal++) {
							if (H_._matrix[row][horizontal] == 1 && horizontal != col) {

								if (SumProduct) { /// Sum Product
									temp *= tanh(Message_info_1._matrix[row][horizontal] / 2);
								}
								else {            /// Min Sum
									if (Message_info_1._matrix[row][horizontal] < 0) ++negative_amount;
									if (abs(Message_info_1._matrix[row][horizontal]) < temp) temp = abs(Message_info_1._matrix[row][horizontal]);
								}
							}
						}
						if (Modified_LLR_c2v) previous = Extrinsis_Mes_1._matrix[row][col];
						if (SumProduct) { ///Sum Product
						temp = log((1 + temp) / (1 - temp));
						if (abs(temp) < LLR_MAX) Extrinsis_Mes_1._matrix[row][col] = temp;
						else Extrinsis_Mes_1._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
						}
						else {           ///Min Sum
							 Extrinsis_Mes_1._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
						}
						if (Modified_LLR_c2v && previous*Extrinsis_Mes_1._matrix[row][col] < 0) {
							Extrinsis_Mes_1._matrix[row][col] += previous;
							Extrinsis_Mes_1._matrix[row][col] *= Modified_LLR_c2v_factor;
						}
						//
						//if (Extrinsis_Mes_1._matrix[row][col] * Message_info_1._matrix[row][col] < 0) Extrinsis_Mes_1._matrix[row][col] *= (Test_parameter);
						//

						Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]);
						//else Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]*0.7);
					}
				}
			}

		}

		LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_1to2, permutation_seq_);

		if (EXIT_Chart_Record && decoding_info.Turbo_Index < decoding_info.LLR_2_extrinsic.size()) {
			temp = 0;
			for (int index = 0; index < Message_Length; index++) {
				temp += Extrinsic_Mes_Passing_2to1_MesPart.at(index);
			}
			decoding_info.LLR_1_priori.at(decoding_info.Turbo_Index) = temp / Message_Length;
			temp = 0;
			for (int index = 0; index < Message_Length; index++) {
				temp += Extrinsic_Mes_Passing_1to2.at(index);
			}
			decoding_info.LLR_1_extrinsic.at(decoding_info.Turbo_Index) = temp/ Message_Length;
		}
		/*
		if (Iteration_Count == Iteration_Number / 4) {
			Final_LLR = 0;
			for (int row = 0; row < Message_Length; ++row) {
				Final_LLR += abs(Channel_LLR_Message_1.at(row) + Extrinsic_Mes_Passing_1to2.at(row) + Extrinsic_Mes_Passing_2to1_MesPart.at(row));
			}
			Previous_AveLLR = Final_LLR;
		}
		//Previous_AveLLR
		if (Iteration_Count == Iteration_Number / 2) {
			Final_LLR = 0;
			for (int row = 0; row < Message_Length; ++row) {
				Final_LLR += abs(Channel_LLR_Message_1.at(row) + Extrinsic_Mes_Passing_1to2.at(row) + Extrinsic_Mes_Passing_2to1_MesPart.at(row));
			}
			Final_LLR /= Message_Length;
			if (Final_LLR - Previous_AveLLR < 0 && Final_LLR < LLR_MAX) {
				cout << Final_LLR << "Inn" << endl;
				for (int i = 0; i < Message_Length; ++i) {
					//Channel_LLR_1.at(i) = 0;
					Extrinsic_Mes_Passing_1to2.at(i) = 0;
				}
				//Channel_LLR_Parity_Temp = Channel_LLR_Parity_1;
				//Channel_LLR_Parity_1 = Channel_LLR_Parity_2;
				//Channel_LLR_Parity_2 = Channel_LLR_Parity_Temp;
				Message_info_1 = ZeroMatrix;
				Message_info_2 = ZeroMatrix;
				Extrinsis_Mes_1 = ZeroMatrix;
				Extrinsis_Mes_2 = ZeroMatrix;
			}
		}*/

		// Decoder 2

		// Le(uk) = L(uk|y) - L(uk) - Lc(yk)
		for (int col = 0; col < Message_Length; col++) {
			Extrinsic_Mes_Passing_1to2_MesPart.at(col) = Extrinsic_Mes_Passing_1to2.at(col);// -Extrinsic_Mes_Passing_2to1_MesPart.at(col);
			//if (Iteration_Count < 100 && Unreliable_Codeword.at(col)) {
				//Extrinsic_Mes_Passing_1to2_MesPart.at(col) = 0;
			//}
		}

		/*
		// Result 1
		cout << endl << Iteration_Count << ":" << endl;
		cout << "D1 to D2" << endl;
		cout << "L(uk):   ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1_MesPart.at(i));
		}
		cout << endl;
		cout << "L(uk|y): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2.at(i) + Channel_LLR_Message_1.at(i));
		}
		cout << endl;
		cout << "Le1(uk): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2_MesPart.at(i));
		}
		cout << endl;*/
		

		// Extrinsic message Le1(uk|y)
		Block_Interleaver(Extrinsic_Mes_Passing_1to2_MesPart);

		//if (Iteration_Count == 0) {
			// LLR = Extrinsic + Channel LLR
			for (int col = 0; col < Message_Length; col++) {
				Channel_LLR_2.at(col) = Channel_LLR_Message_2.at(col) + Extrinsic_Mes_Passing_1to2_MesPart.at(col);  // %50
				//Channel_LLR_2.at(col) = Channel_LLR_Message_2.at(col) * Extrinsic_Mes_Passing_1to2_MesPart.at(col);
			}
			for (int i = Message_Length; i < Col; ++i) {
				Channel_LLR_2.at(i) = Channel_LLR_Parity_2.at(i - Message_Length);
			}
		//}

		/*
		cout << "LLR 2: ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Channel_LLR_2.at(i));
		}
		cout << endl;
		*/
		// variable to check
		LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);

		for (int local_iter = 0; local_iter < Local_Iteration; local_iter++) {

			for (int row = 0; row < Row; ++row) {
				for (int col = 0; col < Col; ++col) {
					if (H._matrix[row][col] == 1) {
						temp = Channel_LLR_2.at(col);
						for (int vertical = 0; vertical < Row; vertical++) {
							if (H._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_2._matrix[vertical][col];
						}
						if (Modified_LLR_v2c) { /// Modified LLR Osc version
							if (temp*Message_info_2._matrix[row][col] < 0) Message_info_2._matrix[row][col] += temp;
							else Message_info_2._matrix[row][col] = temp;
						}
						else {              ///  Conventional version
							Message_info_2._matrix[row][col] = temp;
						}
						Message_info_2._matrix[row][col] *= Modified_LLR_v2c_factor;
					}
				}
			}

			// check to variable
			Extrinsic_Mes_Passing_2to1.assign(ZeroSequence.begin(), ZeroSequence.end());

			for (int row = 0; row < Row; ++row) {
				for (int col = 0; col < Col; ++col) {
					if (H._matrix[row][col] == 1) {
						if (SumProduct) { /// Sum Product
						temp = 1;
						}
						else {            /// Min Sum
							negative_amount = 0;
							temp = DBL_MAX;
						}
						for (int horizontal = 0; horizontal < Col; horizontal++) {
							if (H._matrix[row][horizontal] == 1 && horizontal != col) {
								if (SumProduct) { /// Sum Product
								temp *= tanh(Message_info_2._matrix[row][horizontal] / 2);
								}
								else {            /// Min Sum
									if (Message_info_2._matrix[row][horizontal] < 0) ++negative_amount;
									if (abs(Message_info_2._matrix[row][horizontal]) < temp) temp = abs(Message_info_2._matrix[row][horizontal]);
								}
							}
						}
						if (Modified_LLR_c2v) previous = Extrinsis_Mes_2._matrix[row][col];
						if (SumProduct) { ///Sum Product
							temp = log((1 + temp) / (1 - temp));
							if (abs(temp) < LLR_MAX) Extrinsis_Mes_2._matrix[row][col] = temp;
							else Extrinsis_Mes_2._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
						}
						else {            ///Min Sum
							Extrinsis_Mes_2._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
						}
						if (Modified_LLR_c2v && Extrinsis_Mes_2._matrix[row][col] * previous < 0) {
							Extrinsis_Mes_2._matrix[row][col] += previous;
							Extrinsis_Mes_2._matrix[row][col] *= Modified_LLR_c2v_factor;
						}

						//if (Extrinsis_Mes_2._matrix[row][col] * Message_info_2._matrix[row][col] < 0) Extrinsis_Mes_2._matrix[row][col] *= (Test_parameter);

						Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]);
						//else Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]*0.7);
						//cout << R << ",";
					}
				}
			}

		}
		
		// symdrone check
		if (Syndrome_Check && !(Iteration_Count%5)) {
			bool Past = FALSE;
			for (int bit = 0; bit < Col; ++bit) {
				if (Channel_LLR_2.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Hard_Result.at(bit) = 0;
				else Hard_Result.at(bit) = 1;
			}
			for (int row = 0; row < Row; ++row) {
				symdrone = 0;
				for (int col = 0; col < Col; ++col) {
					symdrone ^= (Hard_Result.at(col)&H._matrix[row][col]);
				}
				if (symdrone) {
					break;
				}
				if (row == Row - 1) Past = 1; //Iteration_Count = Iteration_Number + 1;
			}
			if (Past) {
				Past = FALSE;
				LDPC.Sort_Sequence_Backward(Extrinsic_Mes_Passing_1to2, permutation_seq_);

				for (int bit = 0; bit < Col_; ++bit) {
					if (Channel_LLR_1.at(bit) + Extrinsic_Mes_Passing_1to2.at(bit) > 0) Hard_Result.at(bit) = 0;
					else Hard_Result.at(bit) = 1;
				}
				for (int row = 0; row < Row_; ++row) {
					symdrone = 0;
					for (int col = 0; col < Col_; ++col) {
						symdrone ^= (Hard_Result.at(col)&H_._matrix[row][col]);
					}
					if (symdrone) {
						//cout << row << endl;
						break;
					}
					if (row == Row_ - 1) {
						//cout << "Early Termination!" << Iteration_Count << endl;
						decoding_info.EarlyTerminate = TRUE;
						Iteration_Count = Iteration_Number + 1;
					}
				}
				LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_1to2, permutation_seq);
			}
		}

		// Bit_Flipping
		if (Bit_Flipping && Iteration_Count > Iteration_Number-5) {
			for (int bit = 0; bit < Message_Length; ++bit) {
				AveLLR += abs(Channel_LLR_2.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit));
			}
			AveLLR /= Message_Length;
			if (AveLLR > 30) {
				LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_2to1, permutation_seq);
				break;
			}
			/*
			for (int bit = 0; bit < Message_Length; ++bit) {
				Channel_LLR_2.at(bit) = Channel_LLR_Message_2.at(bit);
			}
			for (int bit = Message_Length; bit < Col; ++bit) {
				Channel_LLR_2.at(bit) = Channel_LLR_Parity_2.at(bit- Message_Length);
			}*/
			for (int bit = 0; bit < Col; ++bit) {
				Channel_LLR_2.at(bit) += Extrinsic_Mes_Passing_2to1.at(bit);
			}

			vector<int> Reliability_Record(Col, 0), Reliable_Row(Row, 0);


			for (int col = 0; col < Col; col++) {
				if (Channel_LLR_2.at(col) > 0) Hard_Result.at(col) = 0;
				else  Hard_Result.at(col) = 1;
			}
			for (int row = 0; row < Row; row++) {
				int syndrome = 0;
				for (int col = 0; col < Col; col++) {
					if (H._matrix[row][col]) syndrome ^= 1;
				}
				if (!syndrome) {
					for (int col = 0; col < Col; col++) {
						if (H._matrix[row][col]) Reliability_Record.at(col)++;
					}
					Reliable_Row.at(row) = 1;
				}
			}
			//bool breaker = FALSE;
			while (1) {
				//breaker = TRUE;
				for (int row = 0; row < Row; row++) {
					if (!Reliable_Row.at(row)) {
						for (int col = 0; col < Col; col++) {
							if (H._matrix[row][col] && Reliability_Record.at(col) != 2) {
								double running_row = (row + 1) % Row;
								while (!H._matrix[running_row][col]) {
									++running_row;
									if (running_row == Row) running_row = 0;
								}
								if (!Reliable_Row.at(running_row)) {
									Channel_LLR_2.at(col) *= (-1);
									Reliable_Row.at(row) = 1;
									Reliable_Row.at(running_row) = 1;
									break;
								}
								else break;

							}
						}
					}
				}
				bool breaker = TRUE;
				/*
				for (int bit = 0; bit < Row; ++bit) {
					if (!Reliable_Row.at(bit)) cout << bit << " Not Fit!" << " ";
				}
				cout << endl;

				for (int col = 0; col < Col; col++) {
					if (Channel_LLR_2.at(col) > 0) Hard_Result.at(col) = 0;
					else  Hard_Result.at(col) = 1;
				}
				for (int row = 0; row < Row; row++) {
					int syndrome = 0;
					for (int col = 0; col < Col; col++) {
						if (H._matrix[row][col]) syndrome ^= 1;
					}
					if (syndrome) {
						cout << row << " ";
					}
				}
				cout << endl << endl;
				*/

				for (int bit = 0; bit < Row; ++bit) {
					if (!Reliable_Row.at(bit)) {
						breaker = FALSE;
						break;
					}
				}
				if (breaker) break;
			}
			//cout << "HEY" << endl;
			LDPC.Sort_Sequence_Forward(Channel_LLR_2, permutation_seq);
			break;
		}
		//cout << "Out" << endl;
		// BF End

		
		LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_2to1, permutation_seq);



		// Test
		if (ML_Likely_Method && (Iteration_Count > 0))
		{
			//cout << codeword_length << endl;
			vector<double> Rx(codeword_length);
			for (int bit = 0; bit < Message_Length; ++bit) {
				Rx.at(bit) = Channel_LLR_Message_2.at(bit) + Extrinsic_Mes_Passing_1to2_MesPart.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit);
				if (Channel_LLR_Message_2.at(bit) + Extrinsic_Mes_Passing_1to2_MesPart.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Now_Hard_Result.at(bit) = 0;
				else Now_Hard_Result.at(bit) = 1;
				//cout << bit << ",";
			}
			for (int bit = Message_Length; bit < Col_; ++bit) {
				Rx.at(bit) = Channel_LLR_Parity_1.at(bit - Message_Length) + Extrinsic_Mes_Passing_1to2.at(bit);
				if (Channel_LLR_Parity_1.at(bit - Message_Length) + Extrinsic_Mes_Passing_1to2.at(bit) > 0) Now_Hard_Result.at(bit) = 0;
				else Now_Hard_Result.at(bit) = 1;
				//cout << bit << ",";
			}
			for (int bit = Message_Length; bit < Col; ++bit) {
				Rx.at(bit + Col_ - Message_Length) = Channel_LLR_Parity_2.at(bit - Message_Length) + Extrinsic_Mes_Passing_2to1.at(bit);
				if (Channel_LLR_Parity_2.at(bit - Message_Length) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Now_Hard_Result.at(bit + Col_- Message_Length) = 0;
				else Now_Hard_Result.at(bit + Col_ - Message_Length) = 1;
				//cout << bit + Col_ - Message_Length << ",";
			}

			if (Good_Hard_Result != Now_Hard_Result) {
				MATRIX<__int8> Sorted_G(G);
				vector <size_t>
					Location_Index(codeword_length, 0);
				MATRIX<double> Metric_Table(2, codeword_length);

				for (int i = 0; i < codeword_length; ++i) {
					if (Rx.at(i)*decoding_info.rx_signal_seq.at(i) < 0) Rx.at(i) = decoding_info.rx_signal_seq.at(i)*(-1);
					else Rx.at(i) = decoding_info.rx_signal_seq.at(i);
				}
				Pre_Procedure(Rx, G, Sorted_G, Location_Index, Metric_Table, decoding_info);
				// Rx呈上G得到剩下的部分
				vector<__int8> M(Message_Length, 0);
				for (int i = 0; i < Message_Length; ++i) {
					if (decoding_info.Sorted_R.at(i) > 0) M.at(i) = 0;
					else M.at(i) = 1;
				}

				int head = 64;
				for (int flip1 = head -1; flip1 < Message_Length; flip1++) {
					if (flip1 == head) {
						M.at(0) ^= 1;
					}
					else if (flip1 > head) {
						M.at(flip1 - 1) ^= 1;
						M.at(flip1) ^= 1;
					}
					for (int flip2 = head - 1; flip2 < Message_Length; flip2++) {
						if (flip1 != flip2) {
							if (flip2 == head) {
								M.at(0) ^= 1;
							}
							else if (flip2 > head) {
								M.at(flip2 - 1) ^= 1;
								M.at(flip2) ^= 1;
							}
						}


						// Observe the parity check bits
						__int8 TempNumber;
						for (int i = 0; i < codeword_length; ++i) {
							TempNumber = 0;
							for (int j = 0; j < Message_Length; ++j) {
								TempNumber ^= (M.at(j)&Sorted_G._matrix[j][i]);
							}
							decoding_info.estimated_codeword.at(i) = TempNumber;
						}
						Desort_Function(Location_Index, decoding_info.estimated_codeword, Now_Hard_Result);

						double metric = 0;
						for (int bit = 0; bit < codeword_length; ++bit) {
							if (Now_Hard_Result.at(bit) != Original_Hard_Result.at(bit))
							{
								metric += abs(decoding_info.rx_signal_seq.at(bit));
							}
						}
						/*
						cout << Iteration_Count << endl;
						for (int bit = 0; bit < codeword_length; ++bit) {
							cout << (Rx.at(bit) > 0 ? 0 : 1) << (int)decoding_info.code_seq.at(bit) << (int)Now_Hard_Result.at(bit) << ", ";
						}
						cout << " (" << metric << ")" << endl << endl;*/
						//cout << "Iteration: " << Iteration_Count << "(" << Best_Hard_Result_Metric << ") " << endl;

						if (metric < Best_Hard_Result_Metric) {
							//cout << "("<< metric << ")";
							Good_Hard_Result = Now_Hard_Result;
							Best_Hard_Result_Metric = metric;
						}
					}
				}
			}
		}

		if (Find_Reliable_Index) {
			for (int bit = 0; bit < Message_Length; ++bit) {
				if (Channel_LLR_Message_2.at(bit) + Extrinsic_Mes_Passing_1to2_MesPart.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Hard_Result.at(bit) = 0;
				else Hard_Result.at(bit) = 1;
			}
			for (int bit = Message_Length; bit < Col; ++bit) {
				if (Channel_LLR_Parity_2.at(bit - Message_Length) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Hard_Result.at(bit) = 0;
				else Hard_Result.at(bit) = 1;
			}

			if (Previous_Hard_Result != Hard_Result) {
				for (int index = 0; index < Message_Length; index++) {
					if (Previous_Hard_Result.at(index) != Hard_Result.at(index)) Bad_index.at(index)++;
				}
				Previous_Hard_Result = Hard_Result;
			}
		}

		// Accumulated LLR
		if (Accumulated_LLR_Calculate) {
			for (int i = 0; i < Message_Length;++i) {
				Accumulated_LLR.at(i) =
					Accumulated_LLR.at(i)*Acc_LLR_Scaling_Factor +
					Extrinsic_Mes_Passing_2to1.at(i) +
					Extrinsic_Mes_Passing_1to2_MesPart.at(i) +
					Channel_LLR_Message_2.at(i);// +Channel_LLR_2.at(i);
			}
			for (int i = Message_Length; i < Col; ++i) {
				Accumulated_LLR.at(i) =
					Accumulated_LLR.at(i)*Acc_LLR_Scaling_Factor +
					Extrinsic_Mes_Passing_2to1.at(i) +
					Channel_LLR_Parity_2.at(i - Message_Length);// +Channel_LLR_2.at(i);
			}

		}

		if (EXIT_Chart_Record && decoding_info.Turbo_Index < decoding_info.LLR_2_extrinsic.size()) {
			temp = 0;
			for (int index = 0; index < Message_Length; index++) {
				temp += Extrinsic_Mes_Passing_1to2_MesPart.at(index);
			}
			decoding_info.LLR_2_priori.at(decoding_info.Turbo_Index) = temp/ Message_Length;
			temp = 0;
			for (int index = 0; index < Message_Length; index++) {
				temp += Extrinsic_Mes_Passing_2to1.at(index);
			}
			decoding_info.LLR_2_extrinsic.at(decoding_info.Turbo_Index) = temp/ Message_Length;
			++decoding_info.Turbo_Index;
		}

		if (Diff_Candidate_Test && (Iteration_Count > (0)) ) {
			//cout << Iteration_Count << ", ";
			// Obtain Soft Sequence
			for (int index = 0; index < Message_Length; index++) {
				Channel_LLR_Message_Temp.at(index) = 
					Channel_LLR_Message_2.at(index) + 
					Extrinsic_Mes_Passing_2to1.at(index) +
					Extrinsic_Mes_Passing_1to2_MesPart.at(index);
			}
			//Block_DeInterleaver(Channel_LLR_Message_Temp);
			for (int index = 0; index < Message_Length; index++) {
				Channel_LLR_CodeWord_Temp.at(index) = Channel_LLR_Message_Temp.at(index);
			}
			for (int index = Message_Length; index < Col; index++) {
				Channel_LLR_CodeWord_Temp.at(index) = Channel_LLR_2.at(index) + Extrinsic_Mes_Passing_2to1.at(index);
			}
			// Hard_Decision
			for (int index = 0; index < Col; index++) {
				if (Channel_LLR_CodeWord_Temp.at(index) > 0) Hard_Result.at(index) = 0;
				else Hard_Result.at(index) = 1;
			}
			if (Current_Hard_Result != Hard_Result) {
				//cout << "IN";
				Current_Hard_Result = Hard_Result;
				bool Check = TRUE;
				for (int i = 0; i < decoding_info.Candidate_Metric.size(); ++i) {
					if (Hard_Result == decoding_info.Iter_Decode_Candidates_HardDecision._matrix.at(i)) {
						Check = FALSE;
						break;
					}
				}
				if (Check == TRUE) {
					decoding_info.Iter_Decode_Candidates_HardDecision._matrix.at(decoding_info.Candidate_Metric.size()) = Current_Hard_Result;
					A_star_PC_out_CBC_OSC_Verified_Function(G, decoding_info, Channel_LLR_CodeWord_Temp, Hard_Result, temp, Average_LLR);
					decoding_info.Iterative_Decoding_Candidates._matrix.at(decoding_info.Candidate_Metric.size()) = Channel_LLR_CodeWord_Temp;
					decoding_info.Iter_Decode_Candidates_Est._matrix.at(decoding_info.Candidate_Metric.size()) = Hard_Result;
					decoding_info.Candidate_Metric.push_back(Average_LLR);
				}
			}
		}
		if (EarlyTermination) {

		}
		/*
		if ((Iteration_Count + 1) % 10 == 0) {
			vector<__int8> CodeWord_Present(192, 0);
			temp = 0;
			for (int index = 0; index < Message_Length; index++) {
				temp += abs(Extrinsic_Mes_Passing_2to1.at(index) +
					Extrinsic_Mes_Passing_1to2_MesPart.at(index) +
					Channel_LLR_Message_2.at(index));
				if (Extrinsic_Mes_Passing_2to1.at(index) +
					Extrinsic_Mes_Passing_1to2_MesPart.at(index) +
					Channel_LLR_Message_2.at(index) > 0) CodeWord_Present.at(index) = 0;
				else CodeWord_Present.at(index) = 1;
			}
			temp /= Message_Length;
			decoding_info.Amplitude_Flipping_number._matrix[9][Record] = temp;

			temp = 0;
			for (int index = Message_Length; index < Col; index++) {
				temp += abs(Extrinsic_Mes_Passing_2to1.at(index) +
					Channel_LLR_Parity_2.at(index - Message_Length));
				if(Extrinsic_Mes_Passing_2to1.at(index) +
					Channel_LLR_Parity_2.at(index - Message_Length) > 0) CodeWord_Present.at(index) = 0;
				else CodeWord_Present.at(index) = 1;
			}
			temp /= (Col - Message_Length);
			decoding_info.Amplitude_Flipping_number._matrix[10][Record] = temp;

			temp = 0;
			for (int index = 0; index < Col; index++) {
				if (CodeWord_Temp.at(index) != CodeWord_Present.at(index)) temp++;
				CodeWord_Temp.at(index) = CodeWord_Present.at(index);
			}
			decoding_info.Amplitude_Flipping_number._matrix[11][Record] = temp;
			Record++;
		}*/

		//
		/*
		if (Iteration_Count == 0) {
			// vector<double> Accumulated_Residual(Col, 0), Previous_Metric(Col, 0);
			// Channel_LLR_Message_2.at(row) += (Extrinsic_Mes_Passing_2to1.at(row) + Extrinsic_Mes_Passing_1to2_MesPart.at(row));
			for (int i = 0; i < Message_Length; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Message_1.at(i) + Extrinsic_Mes_Passing_1to2.at(i) + Extrinsic_Mes_Passing_2to1.at(i);
			}
			for (int i = Message_Length; i < Col_; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Parity_1.at(i - Message_Length) + Extrinsic_Mes_Passing_1to2.at(i);
			}
			//cout << Col_ << "," << codeword_length << "," << Channel_LLR_Parity_2.size() << endl;
			for (int i = Col_; i < codeword_length; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Parity_2.at(i - Col_) + Extrinsic_Mes_Passing_2to1.at(i - Row_);
			}
		}
		else {
			double temmp;
			for (int i = 0; i < Message_Length; ++i) {
				temmp = Channel_LLR_Message_1.at(i) + Extrinsic_Mes_Passing_1to2.at(i) + Extrinsic_Mes_Passing_2to1.at(i);
				Accumulated_Residual.at(i) += abs(temmp - Previous_Metric.at(i));
			}
			for (int i = Message_Length; i < Col_; ++i) {
				temmp = Channel_LLR_Parity_1.at(i - Message_Length) + Extrinsic_Mes_Passing_1to2.at(i);
				Accumulated_Residual.at(i) += abs(temmp - Previous_Metric.at(i));
			}
			for (int i = Col_; i < codeword_length; ++i) {
				temmp = Channel_LLR_Parity_2.at(i - Col_) + Extrinsic_Mes_Passing_2to1.at(i - Row_);
				Accumulated_Residual.at(i) += abs(temmp - Previous_Metric.at(i));
			}

			for (int i = 0; i < Message_Length; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Message_1.at(i) + Extrinsic_Mes_Passing_1to2.at(i) + Extrinsic_Mes_Passing_2to1.at(i);
			}
			for (int i = Message_Length; i < Col_; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Parity_1.at(i - Message_Length) + Extrinsic_Mes_Passing_1to2.at(i);
			}
			for (int i = Col_; i < codeword_length; ++i) {
				Previous_Metric.at(i) = Channel_LLR_Parity_2.at(i - Col_) + Extrinsic_Mes_Passing_2to1.at(i - Row_);
			}
		}*/

		if (Iteration_Count++ > Iteration_Number) {
			//for (int i = 0; i < codeword_length; ++i) {
				//Accumulated_Residual.at(i) /= Iteration_Number;
			//}
			break;
		}
		//cout << "iter:"<< Iteration_Count << " ,";
		// Le(uk) = L(uk|y) - L(uk) - Lc(yk)
		for (int col = 0; col < Message_Length; col++) {
			Extrinsic_Mes_Passing_2to1_MesPart.at(col) = Extrinsic_Mes_Passing_2to1.at(col); // -Extrinsic_Mes_Passing_1to2_MesPart.at(col);
		}

		/*
		// Result 2
		cout << "D2 to D1" << endl;
		cout << "L(uk):   ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2_MesPart.at(i));
		}
		cout << endl;
		cout << "L(uk|y): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1.at(i) + Channel_LLR_Message_1.at(i));
		}
		cout << endl;
		cout << "Le2(uk): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1_MesPart.at(i));
		}
		cout << endl << endl;
		*/
		// Extrinsic message Le2(uk|y)
		Block_DeInterleaver(Extrinsic_Mes_Passing_2to1_MesPart);

		// Decoder 2 end

		// LLR = Extrinsic + Channel LLR
		for (int col = 0; col < Message_Length; col++) {
			Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) + Extrinsic_Mes_Passing_2to1_MesPart.at(col);  //  %50
			//Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) * Extrinsic_Mes_Passing_2to1_MesPart.at(col);
		}
		for (int i = Message_Length; i < Col_; ++i) {
			Channel_LLR_1.at(i) = Channel_LLR_Parity_1.at(i - Message_Length);
		}

	}
	// Iteration End







	for (int row = 0; row < Message_Length; ++row) {
		Final_LLR += abs(Extrinsic_Mes_Passing_2to1.at(row) + Extrinsic_Mes_Passing_1to2_MesPart.at(row));
	}
	Final_LLR /= Message_Length;
	//cout << Best_Hard_Result_Metric << " / " << Original_LLR << ", " << Best_Hard_Result_Metric/ Original_LLR << endl;
	//cout << "/" << Original_LLR << "/" << endl;
	/*
	bool CorrectOne = TRUE;
	vector<__int8> Mes(CodeWord_Temp);
	Mes.erase(Mes.begin() + Message_Length, Mes.end());
	Block_DeInterleaver(Mes);
	for (int index = 0; index < Message_Length; ++index) {
		if (Mes.at(index) != decoding_info.message_seq.at(index)) {
			CorrectOne = FALSE;
			break;
		}
	}
	if (CorrectOne == FALSE) {
		vector<__int8> Codeword(192, 0);
		Mes.assign(CodeWord_Temp.begin(), CodeWord_Temp.begin() + Message_Length);
		Systematic_Linear_Block_Code_Encoder(G, Mes, Codeword);
		if (Codeword == CodeWord_Temp) {
			// Wrong Codeword
			//cout << "Wrong Codeword" << endl;
			for (int index = 0; index < 13; index++) {
				decoding_info.Amplitude_Flipping_number._matrix[3][index] += decoding_info.Amplitude_Flipping_number._matrix[9][index];
				decoding_info.Amplitude_Flipping_number._matrix[4][index] += decoding_info.Amplitude_Flipping_number._matrix[10][index];
				decoding_info.Amplitude_Flipping_number._matrix[5][index] += decoding_info.Amplitude_Flipping_number._matrix[11][index];
			}
			decoding_info.Amplitude_Flipping_number._matrix[3][13]++;
			decoding_info.Amplitude_Flipping_number._matrix[4][13]++;
			decoding_info.Amplitude_Flipping_number._matrix[5][13]++;
		}
		else {
			// Convegence Failed
			//cout << "Convegence Failed" << endl;
			for (int index = 0; index < 13; index++) {
				decoding_info.Amplitude_Flipping_number._matrix[6][index] += decoding_info.Amplitude_Flipping_number._matrix[9][index];
				decoding_info.Amplitude_Flipping_number._matrix[7][index] += decoding_info.Amplitude_Flipping_number._matrix[10][index];
				decoding_info.Amplitude_Flipping_number._matrix[8][index] += decoding_info.Amplitude_Flipping_number._matrix[11][index];
			}
			decoding_info.Amplitude_Flipping_number._matrix[6][13]++;
			decoding_info.Amplitude_Flipping_number._matrix[7][13]++;
			decoding_info.Amplitude_Flipping_number._matrix[8][13]++;
		}

	}
	else {
		//cout << "Success" << endl;
		for (int index = 0; index < 13; index++) {
			decoding_info.Amplitude_Flipping_number._matrix[0][index] += decoding_info.Amplitude_Flipping_number._matrix[9][index];
			decoding_info.Amplitude_Flipping_number._matrix[1][index] += decoding_info.Amplitude_Flipping_number._matrix[10][index];
			decoding_info.Amplitude_Flipping_number._matrix[2][index] += decoding_info.Amplitude_Flipping_number._matrix[11][index];
		}
		decoding_info.Amplitude_Flipping_number._matrix[0][13]++;
		decoding_info.Amplitude_Flipping_number._matrix[1][13]++;
		decoding_info.Amplitude_Flipping_number._matrix[2][13]++;
	}*/

	if (Accumulated_LLR_Calculate) {
		// Estimate
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			Channel_LLR_Message_2.at(row) = Accumulated_LLR.at(row);
		}
		Block_DeInterleaver(Channel_LLR_Message_2);
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_Message_2.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
		// LLR value
		for (int row = 0; row < Col; ++row) {
			decoding_info.rx_signal_seq.at(row) = Accumulated_LLR.at(row)*decoding_info.var / 2;
		}

	}
	else if (Diff_Candidate_Test) {
		for (int index = 0; index < Message_Length; index++) {
			Channel_LLR_Message_Temp.at(index) =
				Channel_LLR_Message_2.at(index) +
				Extrinsic_Mes_Passing_2to1.at(index) +
				Extrinsic_Mes_Passing_1to2_MesPart.at(index);
		}
		//Block_DeInterleaver(Channel_LLR_Message_Temp);
		for (int index = 0; index < Message_Length; index++) {
			Channel_LLR_CodeWord_Temp.at(index) = Channel_LLR_Message_Temp.at(index);
		}
		for (int index = Message_Length; index < Col; index++) {
			Channel_LLR_CodeWord_Temp.at(index) = Channel_LLR_2.at(index) + Extrinsic_Mes_Passing_2to1.at(index);
		}
		
     	{
			A_star_PC_out_CBC_OSC_Verified_Function(G, decoding_info, Channel_LLR_CodeWord_Temp, Hard_Result, temp, Average_LLR);
			decoding_info.Iterative_Decoding_Candidates._matrix.at(decoding_info.Candidate_Metric.size()) = Channel_LLR_CodeWord_Temp;
			decoding_info.Iter_Decode_Candidates_Est._matrix.at(decoding_info.Candidate_Metric.size()) = Hard_Result;
			decoding_info.Candidate_Metric.push_back(Average_LLR);
		}
	}
	else if (Bit_Flipping && AveLLR < LLR_MAX) {
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			Channel_LLR_Message_2.at(row) = (Channel_LLR_2.at(row));
		}
		Block_DeInterleaver(Channel_LLR_Message_2);
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_Message_2.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
	}
	else if (Find_Reliable_Index) {
		int countt = 0;
		for (int index = 0; index < Message_Length; index++) {
			if (Bad_index.at(index) > 0) {
				decoding_info.rx_signal_seq.at(index) *= 100;
				countt++;
			}
		}
		cout << countt << endl;
	}
	else if (ML_Likely_Method && Final_LLR < 10) {
		cout << "In ";
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			decoding_info.estimated_codeword.at(row) = Good_Hard_Result.at(row);
		}
	}
	else {
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			Channel_LLR_Message_2.at(row) += (Extrinsic_Mes_Passing_2to1.at(row) + Extrinsic_Mes_Passing_1to2_MesPart.at(row));
			//Channel_LLR_Message_2.at(row) = Channel_LLR_Message_2.at(row) * Extrinsic_Mes_Passing_2to1.at(row) + Extrinsic_Mes_Passing_1to2_MesPart.at(row);
		}
		Block_DeInterleaver(Channel_LLR_Message_2);
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_Message_2.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
		//decoding_info.Turbo_LLR /= Channel_LLR_Message_2.size();

		//int Error = 0;
		//decoding_info.Turbo_Small_Value_Counter = 0;
		//for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			//if (decoding_info.estimated_codeword.at(row) != decoding_info.message_seq.at(row)) Error++;
			//if (abs(Channel_LLR_Message_2.at(row)) < 10) decoding_info.Turbo_Small_Value_Counter++;
		//}

		//cout << "AverLLR: "<< decoding_info.Turbo_LLR << " /  ";
		/*
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (abs(Channel_LLR_Message_2.at(row)) < 10) cout << Channel_LLR_Message_2.at(row) << ", ";
		}

		cout << endl;
		*/
		/*
		if (Error != 0) {

			//vector<double> Sorted_Rx=

			cout << endl << "Turbo Decoder error: " << Error << endl;
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				if (decoding_info.estimated_codeword.at(row) != decoding_info.message_seq.at(row)) printf("(%d,%.2f) ", row, Channel_LLR_Message_2.at(row)*decoding_info.var/2);
			}
			cout << endl;
			double LLR = 0;
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				LLR += Channel_LLR_Message_2.at(row);
			}
			LLR /= Channel_LLR_Message_2.size();
			cout << "Error LLR: " << LLR << endl <<endl;

		}*/

		// Return LLR
		//if (decoding_info.Turbo_Small_Value_Counter != 0) 
		if (Combine_Two_CodeWords) {
			//if (AveLLR < LLR_MAX) {
			/*
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
					decoding_info.rx_signal_seq.at(row) += (Channel_LLR_Message_2.at(row)*decoding_info.var / 2);
				}

				for (int row = Channel_LLR_Message_2.size(); row < Extrinsic_Mes_Passing_2to1.size(); ++row) {
					decoding_info.rx_signal_seq.at(row) += (decoding_info.rx_signal_seq.at(row) + Extrinsic_Mes_Passing_1to2.at(row)*decoding_info.var / 2);
					decoding_info.rx_signal_seq.at(row + (Message_Length / 2)) +=
						(decoding_info.rx_signal_seq.at(row + (Message_Length / 2)) + (Extrinsic_Mes_Passing_2to1.at(row)*decoding_info.var / 2));
				}
			*/
				for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
					decoding_info.rx_signal_seq.at(row) += (Channel_LLR_Message_2.at(row)*decoding_info.var*CTC_factor / 2);
				}

				for (int row = Channel_LLR_Message_2.size(); row < Extrinsic_Mes_Passing_2to1.size(); ++row) {
					decoding_info.rx_signal_seq.at(row) += ( (Channel_LLR_Parity_1.at(row - Channel_LLR_Message_2.size()) 
						+ Extrinsic_Mes_Passing_1to2.at(row))*decoding_info.var*CTC_factor / 2);
					decoding_info.rx_signal_seq.at(row + (Message_Length / 2)) +=
						( (Channel_LLR_Parity_2.at(row - Channel_LLR_Message_2.size()) + (Extrinsic_Mes_Passing_2to1.at(row))*decoding_info.var*CTC_factor / 2));
				}
				/*
				for (int row = 0; row < codeword_length; ++row) {
					decoding_info.rx_signal_seq.at(row) /= Accumulated_Residual.at(row);
				}*/

			//}
		}
		// Astar_after_TurboDecoding
		if (Astar_after_TurboDecoding)
		{
			Block_Interleaver(Channel_LLR_Message_2);
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				decoding_info.rx_signal_seq.at(row) = Channel_LLR_Message_2.at(row)*decoding_info.var / 2;
			}

			for (int row = Channel_LLR_Message_2.size(); row < Extrinsic_Mes_Passing_2to1.size(); ++row) {
				decoding_info.rx_signal_seq.at(row) = decoding_info.rx_signal_seq.at(row + Row_);
				decoding_info.rx_signal_seq.at(row) += ((Extrinsic_Mes_Passing_2to1.at(row))*decoding_info.var / 2);
			}
		}

	}
}

void LDPC_TS_VC_RBP(MATRIX<__int8> &H_,
	MATRIX<__int8> &H, MATRIX<__int8> &G,
	DECODING_INFO &decoding_info,
	vector<int> &permutation_seq_,
	vector<int> &permutation_seq) {

	int Col = H.Col_number,
		Row = H.Row_number,
		Col_ = H_.Col_number,
		Row_ = H_.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		negative_amount,
		symdrone,
		codeword_length;
	long double temp;
	double Average_LLR, previous, Best_Hard_Result_Metric;
	Best_Hard_Result_Metric = DBL_MAX;

	vector<__int8> Hard_Result(Col, 0),
		Current_Hard_Result(Col, 2),
		Bad_index(Message_Length, 0),
		Previous_Hard_Result(Col, 0),
		Best_Hard_Result(Col, 0);

	vector<double>
		Channel_LLR_1(Col_, 0),
		Channel_LLR_2(Col, 0),
		Extrinsic_Mes_Passing_1to2(Col_, 0),
		Extrinsic_Mes_Passing_2to1(Col, 0),
		Extrinsic_Mes_Passing_1to2_MesPart(Message_Length, 0),
		Extrinsic_Mes_Passing_2to1_MesPart(Message_Length, 0),
		ZeroSequence(Col, 0),
		ZeroSequence_(Col_, 0),
		Accumulated_LLR(Col, 0);

	MATRIX<long double>
		Message_info_1,
		Extrinsis_Mes_1,
		Message_info_2,
		Extrinsis_Mes_2,
		Residual_V2C_1,
		Residual_V2C_2,
		ZeroMatrix;

	Extrinsis_Mes_1.Building_Empty_Matrix(Row_, Col_);
	Message_info_1.Building_Empty_Matrix(Row_, Col_);
	Extrinsis_Mes_2.Building_Empty_Matrix(Row, Col);
	Message_info_2.Building_Empty_Matrix(Row, Col);
	Residual_V2C_1.Building_Empty_Matrix(Row_, Col_);
	Residual_V2C_2.Building_Empty_Matrix(Row, Col);
	ZeroMatrix.Building_Empty_Matrix(Row, Col);

	LDPC_FUNC LDPC;
	codeword_length = Col + Row_;

	for (int i = 0; i < Col_; ++i) {
		Channel_LLR_1.at(i) = 2 * decoding_info.rx_signal_seq.at(i) / decoding_info.var;
	}
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_2.at(i) = Channel_LLR_1.at(i);
	}
	for (int i = Message_Length; i < Col; ++i) {
		Channel_LLR_2.at(i) = 2 * decoding_info.rx_signal_seq.at(i + Row_) / decoding_info.var;
	}

	LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
	LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);

	for (int row = 0; row < Row_; ++row) {
		for (int col = 0; col < Col_; ++col) {
			if (H_._matrix[row][col]) {
				Message_info_1._matrix[row][col] = Channel_LLR_1.at(col);
			}
		}
	}

	for (int row = 0; row < Row; ++row) {
		for (int col = 0; col < Col; ++col) {
			if (H._matrix[row][col]) {
				Message_info_2._matrix[row][col] = Channel_LLR_2.at(col);
			}
		}
	}

	for (int index = 0; index < Col_; ++index) {

	}



	//Residual_V2C_1


	while (1) {

		// Decoder 1
		for (int local_iter = 0; local_iter < Local_Iteration; local_iter++) {




			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						temp = Channel_LLR_1.at(col);
						for (int vertical = 0; vertical < Row_; vertical++) {
							if (H_._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_1._matrix[vertical][col];
						}
						if (Modified_LLR_v2c) {       /// Modified LLR Osc version
							if (temp*Message_info_1._matrix[row][col] < 0) Message_info_1._matrix[row][col] += temp;
							else Message_info_1._matrix[row][col] = temp;
						}
						else {                    ///  Conventional version
							Message_info_1._matrix[row][col] = temp;
						}
						Message_info_1._matrix[row][col] *= Modified_LLR_v2c_factor;
					}
				}
			}

			// check to variable
			Extrinsic_Mes_Passing_1to2.assign(ZeroSequence_.begin(), ZeroSequence_.end());
			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						if (SumProduct) { /// Sum Product
							temp = 1;
						}
						else {            /// Min Sum
							negative_amount = 0;
							temp = DBL_MAX;
						}
						for (int horizontal = 0; horizontal < Col_; horizontal++) {
							if (H_._matrix[row][horizontal] == 1 && horizontal != col) {

								if (SumProduct) { /// Sum Product
									temp *= tanh(Message_info_1._matrix[row][horizontal] / 2);
								}
								else {            /// Min Sum
									if (Message_info_1._matrix[row][horizontal] < 0) ++negative_amount;
									if (abs(Message_info_1._matrix[row][horizontal]) < temp) temp = abs(Message_info_1._matrix[row][horizontal]);
								}
							}
						}
						if (Modified_LLR_c2v) previous = Extrinsis_Mes_1._matrix[row][col];
						if (SumProduct) { ///Sum Product
							temp = log((1 + temp) / (1 - temp));
							if (abs(temp) < LLR_MAX) Extrinsis_Mes_1._matrix[row][col] = temp;
							else Extrinsis_Mes_1._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
						}
						else {           ///Min Sum
							Extrinsis_Mes_1._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
						}
						if (Modified_LLR_c2v && previous*Extrinsis_Mes_1._matrix[row][col] < 0) {
							Extrinsis_Mes_1._matrix[row][col] += previous;
							Extrinsis_Mes_1._matrix[row][col] *= Modified_LLR_c2v_factor;
						}
						//
						//if (Extrinsis_Mes_1._matrix[row][col] * Message_info_1._matrix[row][col] < 0) Extrinsis_Mes_1._matrix[row][col] *= (Test_parameter);
						//

						Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]);
						//else Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]*0.7);
					}
				}
			}

		}

		//LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_1to2, permutation_seq_);










		if (Iteration_Count > Iteration_Number) {
			break;
		}

	}
	LDPC.Sort_Sequence_Forward(Channel_LLR_1, permutation_seq_);
	LDPC.Sort_Sequence_Forward(Channel_LLR_2, permutation_seq);

}

void Turbo_LDPC_Repitition_Method(MATRIX<__int8> &H_,
	MATRIX<__int8> &H, MATRIX<__int8> &G,
	DECODING_INFO &decoding_info,
	vector<int> &permutation_seq_,
	vector<int> &permutation_seq)
{
	int Col = H.Col_number,
		Row = H.Row_number,
		Col_ = H_.Col_number,
		Row_ = H_.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		negative_amount,
		symdrone,
		codeword_length;
	long double temp;
	double Average_LLR, previous;

	vector<__int8> Hard_Result(Col, 0),
		Current_Hard_Result(Col, 2),
		Bad_index(Message_Length, 0),
		Previous_Hard_Result(Col, 0),
		Best_Hard_Result(Col, 0);

	vector<double>
		Channel_LLR_1(Col_, 0),
		Channel_LLR_2(Col, 0),
		Channel_LLR_Message_1(Message_Length, 0),
		Channel_LLR_Message_2(Message_Length, 0),
		Channel_LLR_Message_Temp(Message_Length, 0),
		Channel_LLR_CodeWord_Temp(Col, 0),
		Channel_LLR_Parity_1(Row_, 0),
		Channel_LLR_Parity_2(Row, 0),
		Extrinsic_Mes_Passing_1to2(Col_, 0),
		Extrinsic_Mes_Passing_2to1(Col, 0),
		Extrinsic_Mes_Passing_1to2_MesPart(Message_Length, 0),
		Extrinsic_Mes_Passing_2to1_MesPart(Message_Length, 0),
		ZeroSequence(Col, 0),
		ZeroSequence_(Col_, 0),
		Accumulated_LLR(Col, 0),
		Decoder_Internodes_to_1(Col_, 0),
		Decoder_Internodes_to_2(Col, 0);

	MATRIX<long double>
		Message_info_1,
		Extrinsis_Mes_1,
		Message_info_2,
		Extrinsis_Mes_2,
		Extrinsic_1_to_Internodes,
		Extrinsic_2_to_Internodes,
		Message_1_to_Internodes,
		Message_2_to_Internodes;

	Extrinsis_Mes_1.Building_Empty_Matrix(Row_, Col_);
	Message_info_1.Building_Empty_Matrix(Row_, Col_);
	Extrinsis_Mes_2.Building_Empty_Matrix(Row, Col);
	Message_info_2.Building_Empty_Matrix(Row, Col);

	Extrinsic_1_to_Internodes.Building_Empty_Matrix(Row_, Col_);
	Extrinsic_2_to_Internodes.Building_Empty_Matrix(Row, Col);
	Message_1_to_Internodes.Building_Empty_Matrix(Row_, Col_);
	Message_2_to_Internodes.Building_Empty_Matrix(Row, Col);

	LDPC_FUNC LDPC;
	
	for (int i = 0; i < Col_; ++i) {
		Channel_LLR_1.at(i) = 2 * decoding_info.rx_signal_seq.at(i) / decoding_info.var;
	}
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_Message_1.at(i) = Channel_LLR_1.at(i);
		Channel_LLR_Message_2.at(i) = Channel_LLR_1.at(i);
	}
	Block_Interleaver(Channel_LLR_Message_2);
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_2.at(i) = Channel_LLR_Message_2.at(i);
	}
	for (int i = Message_Length; i < Col; ++i) {
		Channel_LLR_2.at(i) = 2 * decoding_info.rx_signal_seq.at(i + Row_) / decoding_info.var;
		Channel_LLR_Parity_2.at(i - Message_Length) = Channel_LLR_2.at(i);
	}
	for (int i = Message_Length; i < Col_; ++i) {
		Channel_LLR_Parity_1.at(i - Message_Length) = Channel_LLR_1.at(i);
	}
	
	// Iteration Start
	LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
	LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);

	//int Onemoretime;
	while (1) {
		//if (Iteration_Count==0) Onemoretime = 2;
		//else Onemoretime = 1;
		// Decoder 1
		// variable to check

		//LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
		//for (int initialize = 0; initialize < Onemoretime; initialize++) {
			for (int local_iter = 0; local_iter < Local_Iteration; local_iter++) {

				for (int row = 0; row < Row_; ++row) {
					for (int col = 0; col < Col_; ++col) {
						if (H_._matrix[row][col] == 1) {
							temp = Channel_LLR_1.at(col);
							for (int vertical = 0; vertical < Row_; vertical++) {
								if (H_._matrix[vertical][col] == 1) temp -= Extrinsis_Mes_1._matrix[vertical][col];
							}
							if (Modified_LLR_v2c) {       /// Modified LLR Osc version
								if (temp*Message_info_1._matrix[row][col] < 0) Message_info_1._matrix[row][col] += temp;
								else Message_info_1._matrix[row][col] = temp;
							}
							else {                    ///  Conventional version
								Message_info_1._matrix[row][col] = temp;
							}
							Message_info_1._matrix[row][col] *= Modified_LLR_v2c_factor;
						}
					}
				}

				// check to variable
				Extrinsic_Mes_Passing_1to2.assign(ZeroSequence_.begin(), ZeroSequence_.end());
				for (int row = 0; row < Row_; ++row) {
					for (int col = 0; col < Col_; ++col) {
						if (H_._matrix[row][col] == 1) {
							if (SumProduct) { /// Sum Product
								temp = 1;
							}
							else {            /// Min Sum
								negative_amount = 0;
								temp = DBL_MAX;
							}
							for (int horizontal = 0; horizontal < Col_; horizontal++) {
								if (H_._matrix[row][horizontal] == 1 && horizontal != col) {

									if (SumProduct) { /// Sum Product
										temp *= tanh(Message_info_1._matrix[row][horizontal] / 2);
									}
									else {            /// Min Sum
										if (Message_info_1._matrix[row][horizontal] < 0) ++negative_amount;
										if (abs(Message_info_1._matrix[row][horizontal]) < temp) temp = abs(Message_info_1._matrix[row][horizontal]);
									}
								}
							}
							if (Modified_LLR_c2v) previous = Extrinsis_Mes_1._matrix[row][col];
							if (SumProduct) { ///Sum Product
								temp = log((1 + temp) / (1 - temp));
								if (abs(temp) < LLR_MAX) Extrinsis_Mes_1._matrix[row][col] = temp;
								else Extrinsis_Mes_1._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
							}
							else {           ///Min Sum
								Extrinsis_Mes_1._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
							}
							if (Modified_LLR_c2v && previous*Extrinsis_Mes_1._matrix[row][col] < 0) {
								Extrinsis_Mes_1._matrix[row][col] += previous;
								Extrinsis_Mes_1._matrix[row][col] *= Modified_LLR_c2v_factor;
							}
							//
							//if (Extrinsis_Mes_1._matrix[row][col] * Message_info_1._matrix[row][col] < 0) Extrinsis_Mes_1._matrix[row][col] *= (Test_parameter);
							//

							Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]);
							//else Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]*0.7);
						}
					}
				}

			}

			for (int i = 0; i < Col_; ++i) {
				Channel_LLR_1.at(i) += Extrinsic_Mes_Passing_1to2.at(i);
			}

			// variable to check
			//LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);

			for (int local_iter = 0; local_iter < Local_Iteration; local_iter++) {

				for (int row = 0; row < Row; ++row) {
					for (int col = 0; col < Col; ++col) {
						if (H._matrix[row][col] == 1) {
							temp = Channel_LLR_2.at(col);
							for (int vertical = 0; vertical < Row; vertical++) {
								if (H._matrix[vertical][col] == 1) temp -= Extrinsis_Mes_2._matrix[vertical][col];
							}
							if (Modified_LLR_v2c) { /// Modified LLR Osc version
								if (temp*Message_info_2._matrix[row][col] < 0) Message_info_2._matrix[row][col] += temp;
								else Message_info_2._matrix[row][col] = temp;
							}
							else {              ///  Conventional version
								Message_info_2._matrix[row][col] = temp;
							}
							Message_info_2._matrix[row][col] *= Modified_LLR_v2c_factor;
						}
					}
				}

				// check to variable
				Extrinsic_Mes_Passing_2to1.assign(ZeroSequence.begin(), ZeroSequence.end());

				for (int row = 0; row < Row; ++row) {
					for (int col = 0; col < Col; ++col) {
						if (H._matrix[row][col] == 1) {
							if (SumProduct) { /// Sum Product
								temp = 1;
							}
							else {            /// Min Sum
								negative_amount = 0;
								temp = DBL_MAX;
							}
							for (int horizontal = 0; horizontal < Col; horizontal++) {
								if (H._matrix[row][horizontal] == 1 && horizontal != col) {
									if (SumProduct) { /// Sum Product
										temp *= tanh(Message_info_2._matrix[row][horizontal] / 2);
									}
									else {            /// Min Sum
										if (Message_info_2._matrix[row][horizontal] < 0) ++negative_amount;
										if (abs(Message_info_2._matrix[row][horizontal]) < temp) temp = abs(Message_info_2._matrix[row][horizontal]);
									}
								}
							}
							if (Modified_LLR_c2v) previous = Extrinsis_Mes_2._matrix[row][col];
							if (SumProduct) { ///Sum Product
								temp = log((1 + temp) / (1 - temp));
								if (abs(temp) < LLR_MAX) Extrinsis_Mes_2._matrix[row][col] = temp;
								else Extrinsis_Mes_2._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));
							}
							else {            ///Min Sum
								Extrinsis_Mes_2._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));
							}
							if (Modified_LLR_c2v && Extrinsis_Mes_2._matrix[row][col] * previous < 0) {
								Extrinsis_Mes_2._matrix[row][col] += previous;
								Extrinsis_Mes_2._matrix[row][col] *= Modified_LLR_c2v_factor;
							}

							//if (Extrinsis_Mes_2._matrix[row][col] * Message_info_2._matrix[row][col] < 0) Extrinsis_Mes_2._matrix[row][col] *= (Test_parameter);

							Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]);
							//else Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]*0.7);
							//cout << R << ",";
						}
					}
				}

			}

			//LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_2to1, permutation_seq);

			for (int i = 0; i < Col; ++i) {
				Channel_LLR_2.at(i) += Extrinsic_Mes_Passing_2to1.at(i);
			}

		//}

		// Decoder to InterNodes

		//LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
		// Decoder_1_to_Internodes,     Message_1_to_Internodes Extrinsic_1_to_Internodes
		for (int row = 0; row < Row_; ++row) {
			for (int col = 0; col < Col_; ++col) {
				if (H_._matrix[row][col] == 1) {
					temp = Channel_LLR_1.at(col);
					for (int vertical = 0; vertical < Row_; vertical++) {
						if (H_._matrix[vertical][col] == 1) temp -= Extrinsic_1_to_Internodes._matrix[vertical][col];
					}
					Message_1_to_Internodes._matrix[row][col] = temp;
				}
			}
		}
		Decoder_Internodes_to_2.assign(ZeroSequence_.begin(), ZeroSequence_.end());
		for (int row = 0; row < Row_; ++row) {
			for (int col = 0; col < Col_; ++col) {
				if (H_._matrix[row][col] == 1) {
					temp = 1;
					for (int horizontal = 0; horizontal < Col_; horizontal++) {
						if (H_._matrix[row][horizontal] == 1 && horizontal != col) {
							temp *= tanh(Message_1_to_Internodes._matrix[row][horizontal] / 2);
						}
					}
					temp = log((1 + temp) / (1 - temp));
					if (abs(temp) < LLR_MAX) Extrinsic_2_to_Internodes._matrix[row][col] = temp;
					else Extrinsic_2_to_Internodes._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));

					Decoder_Internodes_to_2.at(col) += (Extrinsic_2_to_Internodes._matrix[row][col]);
				}
			}
		}

		for (int row = 0; row < Row_; ++row) {
			for (int col = 0; col < Col_; ++col) {
				if (H_._matrix[row][col] == 1) {
					temp = Channel_LLR_2.at(col);
					for (int vertical = 0; vertical < Row_; vertical++) {
						if (H_._matrix[vertical][col] == 1) temp -= Extrinsic_2_to_Internodes._matrix[vertical][col];
					}
					Message_2_to_Internodes._matrix[row][col] = temp;
				}
			}
		}
		Decoder_Internodes_to_1.assign(ZeroSequence_.begin(), ZeroSequence_.end());
		for (int row = 0; row < Row_; ++row) {
			for (int col = 0; col < Col_; ++col) {
				if (H_._matrix[row][col] == 1) {
					temp = 1;
					for (int horizontal = 0; horizontal < Col_; horizontal++) {
						if (H_._matrix[row][horizontal] == 1 && horizontal != col) {
							temp *= tanh(Message_2_to_Internodes._matrix[row][horizontal] / 2);
						}
					}
					temp = log((1 + temp) / (1 - temp));
					if (abs(temp) < LLR_MAX) Extrinsic_1_to_Internodes._matrix[row][col] = temp;
					else Extrinsic_1_to_Internodes._matrix[row][col] = (temp > 0 ? LLR_MAX : LLR_MAX * (-1));

					Decoder_Internodes_to_1.at(col) += (Extrinsic_1_to_Internodes._matrix[row][col]);
				}
			}
		}

		for (int i = 0; i < Col; ++i) {
			Channel_LLR_1.at(i) += Decoder_Internodes_to_1.at(i);
		}
		for (int i = 0; i < Col; ++i) {
			Channel_LLR_2.at(i) += Decoder_Internodes_to_2.at(i);
		}
		
		if (Iteration_Count++ > Iteration_Number) break;

		/*
		// Le(uk) = L(uk|y) - L(uk) - Lc(yk)
		for (int col = 0; col < Message_Length; col++) {
			Extrinsic_Mes_Passing_2to1_MesPart.at(col) = Extrinsic_Mes_Passing_2to1.at(col); // -Extrinsic_Mes_Passing_1to2_MesPart.at(col);
		}

		// Extrinsic message Le2(uk|y)
		//Block_DeInterleaver(Extrinsic_Mes_Passing_2to1_MesPart);

		// Decoder 2 end

		// LLR = Extrinsic + Channel LLR
		for (int col = 0; col < Message_Length; col++) {
			Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) + Extrinsic_Mes_Passing_2to1_MesPart.at(col);  //  %50
			//Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) * Extrinsic_Mes_Passing_2to1_MesPart.at(col);
		}
		for (int i = Message_Length; i < Col_; ++i) {
			Channel_LLR_1.at(i) = Channel_LLR_Parity_1.at(i - Message_Length);
		}*/
	}
	LDPC.Sort_Sequence_Forward(Channel_LLR_1, permutation_seq_);
	LDPC.Sort_Sequence_Forward(Channel_LLR_2, permutation_seq_);

	// Iteration End

	{
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_2.at(row) + Channel_LLR_1.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
		
	}
}

void Turbo_LDPC_decoder_punctured(MATRIX<__int8> &H_, MATRIX<__int8> &H, DECODING_INFO &decoding_info, vector<int> &permutation_seq_, vector<int> &permutation_seq) {
	//cout << "I am here" << decoding_info.var << endl << endl;
	// Intialization
	int Col = H.Col_number,
		Row = H.Row_number,
		Col_ = H_.Col_number,
		Row_ = H_.Row_number,
		Message_Length = Col - Row,
		Iteration_Count = 0,
		negative_amount,
		symdrone,
		P1_length,
		P2_length;

	P1_length = 64;
	P2_length = 64;

	long double temp;

	vector<__int8> Hard_Result(Col, 0);

	vector<double>
		Channel_LLR_1(Col_, 0),
		Channel_LLR_2(Col, 0),
		Channel_LLR_Message_1(Message_Length, 0),
		Channel_LLR_Message_2(Message_Length, 0),
		Channel_LLR_Parity_1(P1_length, 0),
		Channel_LLR_Parity_2(P2_length, 0),
		Extrinsic_Mes_Passing_1to2(Col_, 0),
		Extrinsic_Mes_Passing_2to1(Col, 0),
		Extrinsic_Mes_Passing_1to2_MesPart(Message_Length, 0),
		Extrinsic_Mes_Passing_2to1_MesPart(Message_Length, 0),
		ZeroSequence(Col, 0),
		ZeroSequence_(Col_, 0),
		Accumulated_LLR(Col, 0);

	MATRIX<long double>
		Message_info_1,
		Extrinsis_Mes_1,
		Message_info_2,
		Extrinsis_Mes_2;

	Extrinsis_Mes_1.Building_Empty_Matrix(Row_, Col_);
	Message_info_1.Building_Empty_Matrix(Row_, Col_);
	Extrinsis_Mes_2.Building_Empty_Matrix(Row, Col);
	Message_info_2.Building_Empty_Matrix(Row, Col);

	LDPC_FUNC LDPC;
	/*
	for (int i = 0; i < Col; ++i) {
		decoding_info.rx_signal_seq.at(i) /= 1000;
	}*/

	//cout << Col << "," << Row << "," << Col_ << "," << Row_ << "," << endl;
	for (int i = 0; i < Message_Length+ P1_length; ++i) {
		Channel_LLR_1.at(i) = 2 * decoding_info.rx_signal_seq.at(i) / decoding_info.var;
	}
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_Message_1.at(i) = Channel_LLR_1.at(i);
		Channel_LLR_Message_2.at(i) = Channel_LLR_1.at(i);
	}
	Block_Interleaver(Channel_LLR_Message_2);
	for (int i = 0; i < Message_Length; ++i) {
		Channel_LLR_2.at(i) = Channel_LLR_Message_2.at(i);
	}
	for (int i = Message_Length; i < Message_Length+ P2_length; ++i) {
		Channel_LLR_2.at(i) = 2 * decoding_info.rx_signal_seq.at(i + P1_length) / decoding_info.var;
		Channel_LLR_Parity_2.at(i - Message_Length) = Channel_LLR_2.at(i);
	}
	for (int i = Message_Length; i < Message_Length +P1_length; ++i) {
		Channel_LLR_Parity_1.at(i - Message_Length) = Channel_LLR_1.at(i);
	}

	//if(LinearBlockCode.Code_number == Turbo_256_128_with_2LDPC_Punc)
	//for (int i = Message_Length * 0.75; i < Col_; ++i) {
		//Channel_LLR_1.at(i) = 0;
		//Channel_LLR_2.at(i) = 0;
		//Channel_LLR_Parity_1.at(i - Message_Length) = 0;
		//Channel_LLR_Parity_2.at(i - Message_Length) = 0;
	//}

	// Iteration Start
	while (1) {
		// Decoder 1
		// variable to check

		LDPC.Sort_Sequence_Backward(Channel_LLR_1, permutation_seq_);
		//for (int local_iter = 0;  local_iter < Local_Iteration; local_iter++) {

			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						temp = Channel_LLR_1.at(col);
						for (int vertical = 0; vertical < Row_; vertical++) {
							if (H_._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_1._matrix[vertical][col];
						}
						Message_info_1._matrix[row][col] = temp;
					}
				}
			}

			// check to variable
			Extrinsic_Mes_Passing_1to2.assign(ZeroSequence_.begin(), ZeroSequence_.end());

			for (int row = 0; row < Row_; ++row) {
				for (int col = 0; col < Col_; ++col) {
					if (H_._matrix[row][col] == 1) {
						/// Sum Product
						//temp = 1;
						/// Min Sum
						negative_amount = 0;
						temp = DBL_MAX;
						for (int horizontal = 0; horizontal < Col_; horizontal++) {
							if (H_._matrix[row][horizontal] == 1 && horizontal != col) {
								/// Sum Product
								//temp *= tanh(Message_info_1._matrix[row][horizontal] / 2);
								/// Min Sum
								if (Message_info_1._matrix[row][horizontal] < 0) ++negative_amount;
								if (abs(Message_info_1._matrix[row][horizontal]) < temp) temp = abs(Message_info_1._matrix[row][horizontal]);
							}
						}
						///Sum Product
						//temp = log((1 + temp) / (1 - temp));
						//if (abs(temp) < LLR_MAX) Extrinsis_Mes_1._matrix[row][col] = temp;
						//else Extrinsis_Mes_1._matrix[row][col] = (temp > 0 ? LLR_MAX: LLR_MAX*(-1));
						///Min Sum
						Extrinsis_Mes_1._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));

						Extrinsic_Mes_Passing_1to2.at(col) += (Extrinsis_Mes_1._matrix[row][col]);
					}
				}
			}

		//}
		LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_1to2, permutation_seq_);

		// Decoder 2

		// Le(uk) = L(uk|y) - L(uk) - Lc(yk)
		for (int col = 0; col < Message_Length; col++) {
			Extrinsic_Mes_Passing_1to2_MesPart.at(col) = Extrinsic_Mes_Passing_1to2.at(col)*Acc_LLR_Scaling_Factor;// -Extrinsic_Mes_Passing_2to1_MesPart.at(col);
		}

		/*
		// Result 1
		cout << endl << Iteration_Count << ":" << endl;
		cout << "D1 to D2" << endl;
		cout << "L(uk):   ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1_MesPart.at(i));
		}
		cout << endl;
		cout << "L(uk|y): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2.at(i) + Channel_LLR_Message_1.at(i));
		}
		cout << endl;
		cout << "Le1(uk): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2_MesPart.at(i));
		}
		cout << endl;*/


		// Extrinsic message Le1(uk|y)
		Block_Interleaver(Extrinsic_Mes_Passing_1to2_MesPart);

		// LLR = Extrinsic + Channel LLR
		for (int col = 0; col < Message_Length; col++) {
			Channel_LLR_2.at(col) = Channel_LLR_Message_2.at(col) + Extrinsic_Mes_Passing_1to2_MesPart.at(col);  // %50
			//Channel_LLR_2.at(col) = Channel_LLR_Message_2.at(col) * Extrinsic_Mes_Passing_1to2_MesPart.at(col) - 1;
		}
		for (int i = Message_Length; i < Message_Length+P2_length; ++i) {
			Channel_LLR_2.at(i) = Channel_LLR_Parity_2.at(i - Message_Length);
		}

		/*
		cout << "LLR 2: ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Channel_LLR_2.at(i));
		}
		cout << endl;
		*/

		// variable to check
		LDPC.Sort_Sequence_Backward(Channel_LLR_2, permutation_seq);
		for (int row = 0; row < Row; ++row) {
			for (int col = 0; col < Col; ++col) {
				if (H._matrix[row][col] == 1) {
					temp = Channel_LLR_2.at(col);
					for (int vertical = 0; vertical < Row; vertical++) {
						if (H._matrix[vertical][col] == 1 && vertical != row) temp += Extrinsis_Mes_2._matrix[vertical][col];
					}
					Message_info_2._matrix[row][col] = temp;
				}
			}
		}

		// check to variable
		Extrinsic_Mes_Passing_2to1.assign(ZeroSequence.begin(), ZeroSequence.end());

		for (int row = 0; row < Row; ++row) {
			for (int col = 0; col < Col; ++col) {
				if (H._matrix[row][col] == 1) {
					/// Sum Product
					//temp = 1;
					/// Min Sum
					negative_amount = 0;
					temp = DBL_MAX;
					for (int horizontal = 0; horizontal < Col; horizontal++) {
						if (H._matrix[row][horizontal] == 1 && horizontal != col) {
							/// Sum Product
							//temp *= tanh(Message_info_2._matrix[row][horizontal] / 2);
							/// Min Sum
							if (Message_info_2._matrix[row][horizontal] < 0) ++negative_amount;
							if (abs(Message_info_2._matrix[row][horizontal]) < temp) temp = abs(Message_info_2._matrix[row][horizontal]);
						}
					}
					///Sum Product
					//temp = log((1 + temp) / (1 - temp));
					//if (abs(temp) < LLR_MAX) Extrinsis_Mes_2._matrix[row][col] = temp;
					//else Extrinsis_Mes_2._matrix[row][col] = (temp > 0 ? LLR_MAX: LLR_MAX*(-1));
					///Min Sum
					Extrinsis_Mes_2._matrix[row][col] = (negative_amount % 2 == 0 ? temp : temp * (-1));

					Extrinsic_Mes_Passing_2to1.at(col) += (Extrinsis_Mes_2._matrix[row][col]);
				}
			}
		}

		// symdrone check
		if (Syndrome_Check) {
			for (int bit = 0; bit < Col; ++bit) {
				if (Channel_LLR_2.at(bit) + Extrinsic_Mes_Passing_2to1.at(bit) > 0) Hard_Result.at(bit) = 0;
				else Hard_Result.at(bit) = 1;
			}
			for (int row = 0; row < Row; ++row) {
				symdrone = 0;
				for (int col = 0; col < Col; ++col) {
					symdrone ^= (Hard_Result.at(col)&H._matrix[row][col]);
				}
				if (symdrone == 1) {
					break;
				}
				if (row == Row - 1) Iteration_Count = Iteration_Number + 1;
			}
		}

		LDPC.Sort_Sequence_Forward(Extrinsic_Mes_Passing_2to1, permutation_seq);

		// Accumulated LLR
		if (Accumulated_LLR_Calculate) {
			for (int i = 0; i < Message_Length; ++i) {
				Accumulated_LLR.at(i) =
					Accumulated_LLR.at(i)*Acc_LLR_Scaling_Factor +
					Extrinsic_Mes_Passing_2to1.at(i) +
					Extrinsic_Mes_Passing_1to2_MesPart.at(i) +
					Channel_LLR_Message_2.at(i);// +Channel_LLR_2.at(i);
			}
			for (int i = Message_Length; i < Col; ++i) {
				Accumulated_LLR.at(i) =
					Accumulated_LLR.at(i)*Acc_LLR_Scaling_Factor +
					Extrinsic_Mes_Passing_2to1.at(i) +
					Channel_LLR_Parity_2.at(i - Message_Length);// +Channel_LLR_2.at(i);
			}

		}

		if (Iteration_Count++ > Iteration_Number) break;

		// Le(uk) = L(uk|y) - L(uk) - Lc(yk)
		for (int col = 0; col < Message_Length; col++) {
			Extrinsic_Mes_Passing_2to1_MesPart.at(col) = Extrinsic_Mes_Passing_2to1.at(col)*Acc_LLR_Scaling_Factor; // -Extrinsic_Mes_Passing_1to2_MesPart.at(col);
		}

		/*
		// Result 2
		cout << "D2 to D1" << endl;
		cout << "L(uk):   ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_1to2_MesPart.at(i));
		}
		cout << endl;
		cout << "L(uk|y): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1.at(i) + Channel_LLR_Message_1.at(i));
		}
		cout << endl;
		cout << "Le2(uk): ";
		for (int i = 0; i < 5; ++i) {
			printf("%.1f ", Extrinsic_Mes_Passing_2to1_MesPart.at(i));
		}
		cout << endl << endl;
		*/
		// Extrinsic message Le2(uk|y)
		Block_DeInterleaver(Extrinsic_Mes_Passing_2to1_MesPart);

		// Decoder 2 end

		// LLR = Extrinsic + Channel LLR
		for (int col = 0; col < Message_Length; col++) {
			Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) + Extrinsic_Mes_Passing_2to1_MesPart.at(col);  //  %50
			//Channel_LLR_1.at(col) = Channel_LLR_Message_1.at(col) * Extrinsic_Mes_Passing_2to1_MesPart.at(col) - 1;
		}
		for (int i = Message_Length; i < Message_Length+P2_length; ++i) {
			Channel_LLR_1.at(i) = Channel_LLR_Parity_1.at(i - Message_Length);
		}
	}
	// Iteration End

	//cout << Channel_LLR_Message_2.at(0) << "," << Extrinsic_Mes_Passing_2to1.at(0) << "," << Extrinsic_Mes_Passing_1to2_MesPart.at(0) << endl;

	if (Accumulated_LLR_Calculate) {
		// Estimate
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			Channel_LLR_Message_2.at(row) = Accumulated_LLR.at(row);
		}
		Block_DeInterleaver(Channel_LLR_Message_2);
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_Message_2.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
		// LLR value
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			decoding_info.rx_signal_seq.at(row) = Accumulated_LLR.at(row)*decoding_info.var / 2;
		}
	}
	else {
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			Channel_LLR_Message_2.at(row) += (Extrinsic_Mes_Passing_2to1.at(row) + Extrinsic_Mes_Passing_1to2_MesPart.at(row));
			//Channel_LLR_Message_2.at(row) = Channel_LLR_Message_2.at(row)*(Extrinsic_Mes_Passing_2to1.at(row) - 1 + Extrinsic_Mes_Passing_1to2_MesPart.at(row));
		}
		Block_DeInterleaver(Channel_LLR_Message_2);
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (Channel_LLR_Message_2.at(row) > 0) decoding_info.estimated_codeword.at(row) = 0;
			else decoding_info.estimated_codeword.at(row) = 1;
		}
		//decoding_info.Turbo_LLR /= Channel_LLR_Message_2.size();

		int Error = 0;
		decoding_info.Turbo_Small_Value_Counter = 0;
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (decoding_info.estimated_codeword.at(row) != decoding_info.message_seq.at(row)) Error++;
			if (abs(Channel_LLR_Message_2.at(row)) < 10) decoding_info.Turbo_Small_Value_Counter++;

		}
		//cout << "AverLLR: "<< decoding_info.Turbo_LLR << " /  ";
		/*
		for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
			if (abs(Channel_LLR_Message_2.at(row)) < 10) cout << Channel_LLR_Message_2.at(row) << ", ";
		}

		cout << endl;
		*/
		/*
		if (Error != 0) {

			//vector<double> Sorted_Rx=

			cout << endl << "Turbo Decoder error: " << Error << endl;
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				if (decoding_info.estimated_codeword.at(row) != decoding_info.message_seq.at(row)) printf("(%d,%.2f) ", row, Channel_LLR_Message_2.at(row)*decoding_info.var/2);
			}
			cout << endl;
			double LLR = 0;
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				LLR += Channel_LLR_Message_2.at(row);
			}
			LLR /= Channel_LLR_Message_2.size();
			cout << "Error LLR: " << LLR << endl <<endl;

		}*/

		// Return LLR
		//if (decoding_info.Turbo_Small_Value_Counter != 0) 
		if (1)
		{
			Block_Interleaver(Channel_LLR_Message_2);
			for (int row = 0; row < Channel_LLR_Message_2.size(); ++row) {
				decoding_info.rx_signal_seq.at(row) = Channel_LLR_Message_2.at(row)*decoding_info.var / 2;
			}

			for (int row = Channel_LLR_Message_2.size(); row < Channel_LLR_Message_2.size() + P2_length; ++row) {
				decoding_info.rx_signal_seq.at(row) = decoding_info.rx_signal_seq.at(row + P1_length);
				decoding_info.rx_signal_seq.at(row) += ((Extrinsic_Mes_Passing_2to1.at(row + P1_length))*decoding_info.var / 2);
			}
		}
	}
}

void Turbo_Encoder_Punctured(vector<__int8> &message_seq, vector<__int8> &output_codeword_seq, int row, int col) {  // Rate = 1/2
	
	vector<__int8> parity_check_code_1(row, 0), parity_check_code_2(row, 0), Interleaved_message(row, 0);
	// g1 = (1,1,1), g2 = (1,0,1)
	__int8 a_k = 0, a_k1 = 0, a_k2 = 0;

	// RSC Encoder 1
	for (int mes = 0; mes < message_seq.size(); mes++) {
		a_k = message_seq.at(mes) ^ a_k1 ^ a_k2;
		parity_check_code_1.at(mes) = a_k ^ a_k2;
		cout << (int)a_k << "," << (int)a_k1 << "," << (int)a_k2 << endl;
		a_k2 = a_k1;
		a_k1 = a_k;
	}
	// Terminate to all zero state
	a_k = a_k1 ^ a_k2;
	parity_check_code_1.push_back(a_k ^ a_k2);
	a_k = a_k1 ^ a_k2;
	parity_check_code_1.push_back(a_k ^ a_k2);

	// Interleaver
	int Interleaver_Row = ROUND_2_INT(sqrt(row)), Interleaver_Col = row / Interleaver_Row, temp = 0;
	MATRIX<__int8> Interleaver;
	Interleaver.Building_Empty_Matrix(Interleaver_Row, Interleaver_Col);

	for (int inter_row = 0; inter_row < Interleaver_Row; ++inter_row) {
		for (int inter_col = 0; inter_col < Interleaver_Col; ++inter_col) {
			Interleaver._matrix[inter_row][inter_col] = message_seq.at(temp);
			++temp;
		}
	}
	temp = 0;
	for (int inter_col = 0; inter_col < Interleaver_Col; ++inter_col) {
		for (int inter_row = 0; inter_row < Interleaver_Row; ++inter_row) {
			Interleaved_message.at(temp) = Interleaver._matrix[inter_row][inter_col];
			++temp;
		}
	}

	// RSC Encoder 2
	for (int mes = 0; mes < message_seq.size(); mes++) {
		a_k = Interleaved_message.at(mes) ^ a_k1 ^ a_k2;
		parity_check_code_2.at(mes) = a_k ^ a_k2;
		a_k2 = a_k1;
		a_k1 = a_k;
	}
	// Terminate to all zero state
	a_k = a_k1 ^ a_k2;
	parity_check_code_2.push_back(a_k ^ a_k2);
	a_k = a_k1 ^ a_k2;
	parity_check_code_2.push_back(a_k ^ a_k2);

	// Puncture !
	output_codeword_seq = message_seq;
	output_codeword_seq.push_back(0);
	output_codeword_seq.push_back(0);
	Remove_Even_Index(parity_check_code_1);
	Remove_Odd_Index(parity_check_code_2);
	output_codeword_seq.insert(output_codeword_seq.end(), parity_check_code_1.begin(), parity_check_code_1.end());
	output_codeword_seq.insert(output_codeword_seq.end(), parity_check_code_2.begin(), parity_check_code_2.end());
}


void BCJR_Decoder(vector<double> Rx, double Eb_No, double code_rate) {
	int message_length = Rx.size() * code_rate;
	double Lc = 4 * code_rate*Eb_No * 2;  //  Q: *2 ???

	MATRIX<double> alpha, beta;
	vector<double> Rx_message_bits(Rx.begin(), Rx.begin() + message_length - 1),
		Rx_parity_check_1(Rx.begin() + message_length, Rx.begin() + Rx.size() * 0.75 - 1),
		Rx_parity_check_2(Rx.begin() + Rx.size() * 0.75, Rx.begin() + Rx.size() - 1);

	// Intialization
	alpha.Building_Empty_Matrix(4, message_length + 1);
	beta.Building_Empty_Matrix(4, message_length + 1);
	Zero_Padding_Even(Rx_parity_check_1);
	Zero_Padding_Odd(Rx_parity_check_2);
	alpha._matrix[0][0] = 1;
	

	Iteration_Number;

}

// Some Tools

void Remove_Even_Index(vector<__int8> &Vec) {
	int size = Vec.size() / 2;
	for (int i = 0; i < size; ++i) {
		Vec.erase(Vec.begin() + 1 + i);
	}
}

void Remove_Odd_Index(vector<__int8> &Vec) {
	int size = Vec.size() / 2;
	for (int i = 0; i < size; ++i) {
		Vec.erase(Vec.begin() + i);
	}
}

void Zero_Padding_Even(vector<double> &Vec) {
	int size = Vec.size();
	for (int i = 0; i < size; ++i) {
		Vec.insert(Vec.begin() + 2 * i + 1, 0);
	}
}

void Zero_Padding_Odd(vector<double> &Vec) {
	int size = Vec.size();
	for (int i = 0; i < size; ++i) {
		Vec.insert(Vec.begin() + 2 * i, 0);
	}
}

double gamma(int State_Before, int State_After, double Rx1, double Rx2, double Lc ) {
	int x1, x2;
	if ((State_Before == 0 && State_After == 0) || (State_Before == 1 && State_After == 2))
	{
		x1 = -1; x2 = -1;
	}
	else if ((State_Before == 2 && State_After == 3) || (State_Before == 3 && State_After == 1))
	{
		x1 = -1; x2 = 1;
	}
	else if ((State_Before == 0 && State_After == 2) || (State_Before == 1 && State_After == 0))
	{
		x1 = 1; x2 = 1;
	}
	else
	{
		x1 = 1; x2 = -1;
	}

	//r = exp(x1*Lu[i] * 0.5 + Lc * 0.5 * (x1*v1 + x2 * v2));

	//double r=exp()
	return 0;
}

void Block_Interleaver(vector<__int8> &input) {
	if (Interleaver_On) {
		int M = sqrt(input.size());
		int N = input.size() / M;
		MATRIX<__int8> Block;
		Block.Building_Empty_Matrix(M, N);

		for (int row = 0; row < M; ++row) {
			for (int col = 0; col < N; ++col) {
				Block._matrix[row][col] = input.at(row*N + col);
			}
		}
		for (int col = 0; col < N; ++col) {
			for (int row = 0; row < M; ++row) {
				input.at(row + col * M) = Block._matrix[row][col];
			}
		}
	}
}

void Block_Interleaver(vector<double> &input) {
	if (Interleaver_On) {
		int M = sqrt(input.size());
		int N = input.size() / M;
		MATRIX<double> Block;
		Block.Building_Empty_Matrix(M, N);

		for (int row = 0; row < M; ++row) {
			for (int col = 0; col < N; ++col) {
				Block._matrix[row][col] = input.at(row*N + col);
			}
		}
		for (int col = 0; col < N; ++col) {
			for (int row = 0; row < M; ++row) {
				input.at(row + col * M) = Block._matrix[row][col];
			}
		}
	}
}

void Block_DeInterleaver(vector<__int8> &input) {
	if (Interleaver_On) {
		int M = sqrt(input.size());
		int N = input.size() / M;
		MATRIX<__int8> Block;
		Block.Building_Empty_Matrix(M, N);

		for (int col = 0; col < N; ++col) {
			for (int row = 0; row < M; ++row) {
				Block._matrix[row][col] = input.at(row + col * M);
			}
		}
		for (int row = 0; row < M; ++row) {
			for (int col = 0; col < N; ++col) {
				input.at(row*N + col) = Block._matrix[row][col];
			}
		}
	}
}

void Block_DeInterleaver(vector<double> &input) {
	if (Interleaver_On) {
		int M = sqrt(input.size());
		int N = input.size() / M;
		MATRIX<double> Block;
		Block.Building_Empty_Matrix(M, N);

		for (int col = 0; col < N; ++col) {
			for (int row = 0; row < M; ++row) {
				Block._matrix[row][col] = input.at(row + col * M);
			}
		}
		for (int row = 0; row < M; ++row) {
			for (int col = 0; col < N; ++col) {
				input.at(row*N + col) = Block._matrix[row][col];
			}
		}
	}
}
