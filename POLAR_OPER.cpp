#include "PinAstar/POLAR_OPER.h"


void Polar_Code_G_Generator(MATRIX<__int8> &G_, vector<int>&frozen, vector<int>&non_frozen, vector<double> &Channel_parameter) {

	//double BEC_error = BEC_Error;
	int Row = G_.Col_number, Col = G_.Row_number;
	//cout << Row << "," << Col << endl;  // 64,128
	//system("pause");
	vector<double> X(Col, BEC_Error);
	vector<double> Y(Col, BEC_Error);
	//vector<int> G_order(Col, 0);
	non_frozen.resize(Row);
	MATRIX<double> s;
	MATRIX<int> RM;

	s.Building_Empty_Matrix(2, Col);
	RM.Building_Empty_Matrix(Col, Col);

	for (int i = 0; i < Col; i++)   // Record Order
	{
		frozen.push_back(0);
	}
	/*
	for (int i = 0; i < Col; i++)   // Record Order
	{
		G_order[i] = i;
	}*/

	for (int i = 0; i < log(Col) / log(2); i++)
	{

		int GG = pow(2, i);
		for (int g = 0; g < GG; g++)
		{
			int p = pow(2, i + 1);
			for (int j = 0; j < Col / p; j++)
			{
				X.at(g*(Col / GG) + (Col / p) + j) = Y.at(g*(Col / GG) + j) * Y.at(g*(Col / GG) + j);
				//	cout << g * (N / G)+(N / p) + j << " * " << g * (N / G)+ j << " ** ";

				X.at(g*(Col / GG) + j) = 2 * Y.at(g * (Col / GG) + j) - X.at(g * (Col / GG) + Col / p + j);
				//X.at(g*(Col / GG) + j) = 2 * Y.at(g * (Col / GG) + j) - pow(X.at(g * (Col / GG) + Col / p + j), 2);
				//	cout << g * (N / G) +  j << " *** " << g * (N / G) +  j << " **** " << g * (N / G) +  N / p + j << " ***** ";
			}
		}
		for (int k = 0; k < Col; k++)
		{
			Y.at(k) = X.at(k);
		}
	}

	/////////////////////////////////////////////////////sorting
	//double s[2][Col];

	for (int k = 0; k < Col; k++)
	{
		//	cout << X[k] << "\t";
		s._matrix[0][k] = X.at(k);           // s[a][b]: a = 0: 紀錄Bhattacharyya參數      (?)
		s._matrix[1][k] = k;                 //          a = 1: 紀錄位置
		//	cout << s[1][k] << endl;
	}

	double tmp;
	for (int i = 0; i < Col; i++)
	{
		for (int j = 0; j < Col - i - 1; j++)

			if (s._matrix[0][j] > s._matrix[0][j + 1])
			{
				tmp = s._matrix[0][j];
				s._matrix[0][j] = s._matrix[0][j + 1];
				s._matrix[0][j + 1] = tmp;
				tmp = s._matrix[1][j];
				s._matrix[1][j] = s._matrix[1][j + 1];
				s._matrix[1][j + 1] = tmp;
			}
	}

	/////////////////////////////////////////////////reserve data_bits

	//vector<double> Channel_Parameter()
	Channel_parameter.resize(Row);
	for (int k = 0; k < Row; k++)  // 128->192
	{
		non_frozen.at(k) = s._matrix[1][k];    // s越小可靠度越高
		Channel_parameter.at(k) = (1 - s._matrix[0][k]);
		frozen.at(non_frozen.at(k)) = 1;
		//cout << s._matrix[0][k] << " ";
	}
	/*
	for (int k = 0; k < Row; k++)  // 128->192
	{
		if (k < 7)   Channel_parameter.at(k) = 3;
		else if (k < 28) Channel_parameter.at(k) = 2;
		else  if(k>128)Channel_parameter.at(k) = 0.5;
		else Channel_parameter.at(k) = 1;
	}*/
	/*
	for (int k = Row+1; k < Col; k++)  // 128->192
	{
		cout << s._matrix[0][k] << " ";
	}*/

	//system("pause");
	//  從小排到大
	
	for (int k = Row-1; k >0; k--)  
	{
		double temp;
		int inttemp;
		for (int j = 0; j <k; j++)
		{
			if (non_frozen.at(j) > non_frozen.at(j + 1)) {
				non_frozen.at(j) ^= non_frozen.at(j + 1);
				non_frozen.at(j + 1) ^= non_frozen.at(j);
				non_frozen.at(j) ^= non_frozen.at(j + 1);
				temp = Channel_parameter.at(j);
				Channel_parameter.at(j) = Channel_parameter.at(j + 1);
				Channel_parameter.at(j + 1) = temp;
			}
		}
	}
	/*
	for (int k = 0; k < Row; k++)  // 128->192
	{
		cout << Channel_parameter.at(k) << " ";
	}
	system("pause");
	*/
	//sort(non_frozen.begin(), non_frozen.begin() + Row);  // 從小排到大
	/*
	for (int i = 0; i < Col; i++)
		if (frozen[i] == 0)
			cout << i << " ";
	cout << endl;
	*/
	//system("pause");
	//cout << endl;

	/////////////////////////////////////////create general generater matix(RM code)

	int n = log(Col) / log(2);
	int gg;

	RM._matrix[0][0] = 1;

	for (int i = 0; i < n; i++)
	{
		gg = pow(2, i);
		for (int j = 0; j < pow(2, i); j++)
		{
			for (int k = 0; k < pow(2, i); k++)
			{
				RM._matrix[j + gg][k] = RM._matrix[j][k];
				RM._matrix[j + gg][k + gg] = RM._matrix[j][k];
			}
		}
	}
	//cout << "Hey";

	////////////////////////////////////////////////print data matrix


	//ofstream ofile("G_nonforzen.txt", ios::out);


	//if (ofile.is_open())
	//{
		//ofile << Row << endl;
		//ofile << Col << endl;
	for (int i = 0; i < Row; i++)  //128
	{
		for (int j = 0; j < Col; j++)   //256
		{
			//cout << non_frozen[i] << endl;
			//ofile << RM._matrix[non_frozen[i]][j];
			G_._matrix[j][i] = RM._matrix[non_frozen[i]][j];
		}
		//ofile << endl;
	}
	/*
	cout << endl << endl;
	for (int i = 0; i < Col; i++)  //128
	{
		for (int j = 0; j < Col - 30; j++)   //256
		{
			cout << (int)RM._matrix[i][j];
		}
		cout << endl;
		//ofile << endl;
	}
	system("pause");*/

	//}
	//else
		//cout << "寫入失敗";
    //cout << "???";
}

void Hybrid_Polar_Astar_Encoder(
	MATRIX<__int8> &G, 
	MATRIX<__int8> &G_Inner, 
	vector<__int8> &message, 
	vector<__int8> &codeword,
	DECODING_INFO &decoding_info) {
	vector<__int8> Outer_Mes(message.end() - Outer_Astar_Length / 2, message.end());
	vector<__int8> Outer_CodeWord(Outer_Astar_Length, 0);
	Systematic_Linear_Block_Code_Encoder(G_Inner, Outer_Mes, Outer_CodeWord);
	vector<__int8> Inner_Mes(message.begin(), message.end());
	Inner_Mes.insert(Inner_Mes.end(), Outer_CodeWord.begin() + Outer_Astar_Length / 2, Outer_CodeWord.end());
	for (int i = 0; i < G.Col_number; ++i) {
		int temp = 0;
		for (int j = 0; j < G.Row_number; j++) {
			temp ^= (Inner_Mes.at(j)&G._matrix[j][i]);
		}
		codeword.at(i) = temp;
	}
}

void Hybrid_Polar_Astar_Decoder(
	MATRIX<__int8> &G,
	MATRIX<__int8> &G_Inner,
	DECODING_INFO &decoding_info) {

	int Row = G.Row_number, Col = G.Col_number;
	int log_Col = log2(Col);

	vector<int>u(Col, 0);
	vector<int>x(Col, 0);
	MATRIX <double>R;
	R.Building_Empty_Matrix(log_Col + 1, Col);
	MATRIX <double>L;
	L.Building_Empty_Matrix(log_Col + 1, Col);
	// Data 進入Polar Decoder
	for (int i = 0; i < Col; i++) //*******************************************************
	{
		L._matrix[0][i] = 2 * decoding_info.rx_signal_seq[i] / decoding_info.var;
	}

	//Decoding Start
	for (int i = 0; i < Col; i++)
	{
		if (decoding_info.frozen[i] == 0)
		{
			R._matrix[log_Col][i] = DBL_MAX;
		}
	}

	int upper, lower;
	short int reg = 0, triger = 1;
	vector<double> Previous_LLR(Col, 0), FlipRecord(Col, 0);

	for (int k = 0; k < BP_Decode_Iter; k++)
	{
		for (int i = 0; i < log_Col; i++)
		{
			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					L._matrix[i + 1][upper] = min_sum(L._matrix[i][lower] + R._matrix[i + 1][lower], L._matrix[i][upper]);
					L._matrix[i + 1][lower] = L._matrix[i][lower] + min_sum(L._matrix[i][upper], R._matrix[i + 1][upper]);
				}
			}
		}


		if (k == 0) {
			vector<double> Outer_Rx(Outer_Astar_Length, 0), Rx_Temp(decoding_info.rx_signal_seq);

			for (int i = (Outer_Astar_Length - 1); i > -1; i--)
			{  
				Outer_Rx.at(i) =
					(R._matrix[log_Col][decoding_info.non_frozen.at(i - Outer_Astar_Length + Row)] +
						L._matrix[log_Col][decoding_info.non_frozen.at(i - Outer_Astar_Length + Row)]);
			}
			/*
			for (int i = 0; i < Outer_Astar_Length; ++i) {
				Outer_Rx.at(i) *= (1 + ((Outer_Astar_Length / 2 - i) / Outer_Astar_Length));
			}*/
			decoding_info.rx_signal_seq.clear();
			decoding_info.rx_signal_seq = Outer_Rx;

			A_star_PC_out_CBC_OSC_Verified(G_Inner, decoding_info);
			for (int i = Outer_Astar_Length - 1; i > -1; i--)
			{
				if ((decoding_info.estimated_codeword.at(i) - 0.5)*(Outer_Rx.at(i)) > 0) {
					cout << i << ",";
					L._matrix[log_Col][decoding_info.non_frozen.at(i - Outer_Astar_Length + Row)] *= (-1);
					//R._matrix[log_Col][decoding_info.non_frozen.at(i - Outer_Astar_Length + Row)] = 0;
				}
			}
			cout << endl;
			decoding_info.rx_signal_seq.clear();
			decoding_info.rx_signal_seq = Rx_Temp;

		}

		for (int i = log_Col - 1; i >= 0; i--)  // -1 ?
		{
			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					R._matrix[i][upper] = min_sum(R._matrix[i + 1][lower] + L._matrix[i][lower], R._matrix[i + 1][upper]);
					R._matrix[i][lower] = min_sum(R._matrix[i + 1][upper], L._matrix[i][upper]) + R._matrix[i + 1][lower];

				}
			}
		}
	}
	


	for (int i = 0; i < Row - Outer_Astar_Length / 2; i++) //*******************************************************   // N/2是對應coderate為1/2
	{
		if (L._matrix[log_Col][decoding_info.non_frozen[i]] + R._matrix[log_Col][decoding_info.non_frozen[i]] >= 0)
		{
			decoding_info.estimated_codeword[i] = 0; // 結果
		}
		else
		{
			decoding_info.estimated_codeword[i] = 1;
		}
		
	}

}
/*
int Sub_SC_Decoder(int layer, int path) //using tree search to decoding 
{
	if (layer == log(N) / log(2) + 1) //the final node of binary tree
		return 0;

	int index = N / pow(2, layer);//lower layer index
	int index2 = pow(2, logN - layer + 1)*(layer_number[layer] / 2);//up layer index


	if (path == 0 && layer != 0)//up layer
	{
		for (int i = 0; i < index; i++)
		{
			X_sc[layer][i + index2] = step(X_sc[layer - 1][i + index2])*step(X_sc[layer - 1][i + index2 + index])*min(abs(X_sc[layer - 1][i + index2]), abs(X_sc[layer - 1][i + index2 + index]));
			//X_sc[layer][i + index2] = f(X_sc[layer - 1][i + index2],X_sc[layer - 1][i + index2 + index]);
			//cout << X_sc[layer][i + index2] << " ";
		}
	}
	if (path == 1)//lower layer
	{
		for (int i = 0; i < index; i++)
		{
			//	cout << layer_bits[layer][i + index2];
			//		system("pause");
			X_sc[layer][i + index2 + index] = pow(-1, layer_bits[layer][i + index2])*X_sc[layer - 1][i + index2] + X_sc[layer - 1][i + index2 + index];
		}
	}




	if (layer == logN) //create table of "u hat" for lower layer using
	{

		int layer_bits_number = cul(layer_number[logN] + 1);

		int gg, reg;
		RM_layer[0][0] = 1;


		//According to different " u hat " of final layer , create "u hat" of the other layer
		for (int i = 0; i < layer_bits_number; i++)
		{
			gg = pow(2, i);
			for (int j = 0; j < pow(2, i); j++)   //create RM code
			{
				for (int k = 0; k < pow(2, i); k++)
				{
					RM_layer[j + gg][k] = RM_layer[j][k];
					RM_layer[j + gg][k + gg] = RM_layer[j][k];
				}
			}

			for (int k = pow(2, i) - 1; k >= 0; k--) //(RM code matrix) * (final layer "u hat")=other layer "u hat"
			{
				reg = 0;
				for (int j = 0; j < pow(2, i); j++)
				{
					//cout << layer_number[logN] - j;
					reg = step3(X_sc[logN][layer_number[logN] - j], layer_number[logN] - j)* RM_layer[k][j] + reg;
					//reg = step4(X_sc[logN][layer_number[logN] - j][list], layer_number[logN] - j,list)* RM_layer[k][j] + reg;
				}

				layer_bits[logN - i][layer_number[logN] - k] = reg % 2;
			}
		}

	}


	layer_number[layer]++;
	Sub_SC_Decoder(layer + 1, 0);
	Sub_SC_Decoder(layer + 1, 1);
}

void Polar_Code_SC_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info)
{
	int Row = G.Row_number, Col = G.Col_number;
	int log_Col = log2(Col);

	vector<int>u(Col, 0);
	vector<int>x(Col, 0);

	Sub_SC_Decoder(0, 0);
	vector<int> message_output(Row, 0);
	for (int i = 0; i < Row; i++)
	{
		if (X_sc[logN][non_frozen[i]] <= 0)
		{
			error_bits++;
		}
	}
	for (int i = 0; i < logN + 1; i++)
	{
		layer_number[i] = 0;
	}
}*/


void Polar_code_BP_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {

	int Row = G.Row_number, Col = G.Col_number;
	int log_Col = log2(Col);

	vector<int>u(Col, 0);
	vector<int>x(Col, 0);
	/*
	cout << endl;
	for (int i = 0; i < Col; ++i) {
		cout << setw(3) << decoding_info.frozen.at(i);
	}
	cout << endl;
	for (int i = 0; i < Col; ++i) {
		cout << setw(3) << decoding_info.non_frozen.at(i);
	}
	cout << endl;
	system("pause");
	*/
	//vector<int> message_output(Row);
	//cout << Row << "," << Col << endl;
	//cout << decoding_info.frozen.size() << "," << decoding_info.non_frozen.size() << endl;

	MATRIX <double>R;
	R.Building_Empty_Matrix(log_Col + 1, Col);
	MATRIX <double>L;
	L.Building_Empty_Matrix(log_Col + 1, Col);

	/*
	for (int i = 0; i < decoding_info.frozen.size(); ++i) {
		cout << decoding_info.frozen.at(i) << " ";
	}
	cout << endl;
	for (int i = 0; i < decoding_info.non_frozen.size(); ++i) {
		cout << decoding_info.non_frozen.at(i) << " ";
	}
	cout << endl << endl;
	*/
	// Data 進入Polar Decoder
	for (int i = 0; i < Col; i++) //*******************************************************
	{
		//X_sc[0][i] = 2 * codeword_awgn[i] / b;
		L._matrix[0][i] = 2 * decoding_info.rx_signal_seq[i] / decoding_info.var;
		//L[0][i] = 0.25;
		//	X_sc[0][i] = codeword_awgn[i];
		//	cout << X_sc[0][i] << " ";
		//cout << decoding_info.var <<endl;
	}


	//Decoding Start
	for (int i = 0; i < Col; i++)
	{
		if (decoding_info.frozen[i] == 0)
		{
			R._matrix[log_Col][i] = DBL_MAX;
		}
		//else R._matrix[log_Col][i] = 1;
		//	cout << R[logN][i] << endl;
	}

	int upper, lower;
	short int reg = 0, triger = 1;
	vector<double> Previous_LLR(Col, 0), FlipRecord(Col, 0);

	for (int k = 0; k < BP_Decode_Iter; k++)
	{
		//cout << "Iteration:" << k << endl;
		for (int i = 0; i < log_Col; i++)
		{

			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					L._matrix[i + 1][upper] = min_sum(L._matrix[i][lower] + R._matrix[i + 1][lower], L._matrix[i][upper]);
					L._matrix[i + 1][lower] = L._matrix[i][lower] + min_sum(L._matrix[i][upper], R._matrix[i + 1][upper]);
					/*L[i + 1][upper] = f(L[i][lower] + R[i + 1][lower], L[i][upper]);
					L[i + 1][lower] = L[i][lower] + f(L[i][upper], R[i + 1][upper]);
					R[i][upper] = f(R[i + 1][lower] + L[i][lower], R[i + 1][upper]);
					R[i][lower] = f(R[i + 1][upper], L[i][upper]) + R[i + 1][lower];*/
				}
			}
		}

		//for (int i = 0; i < log_Col; i++)
		//{
		for (int i = log_Col; i >= 0; i--)
		{
			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					R._matrix[i][upper] = min_sum(R._matrix[i + 1][lower] + L._matrix[i][lower], R._matrix[i + 1][upper]);
					R._matrix[i][lower] = min_sum(R._matrix[i + 1][upper], L._matrix[i][upper]) + R._matrix[i + 1][lower];

					/*L[i + 1][upper] = f(L[i][lower] + R[i + 1][lower], L[i][upper]);
					L[i + 1][lower] = L[i][lower] + f(L[i][upper], R[i + 1][upper]);
					R[i][upper] = f(R[i + 1][lower] + L[i][lower], R[i + 1][upper]);
					R[i][lower] = f(R[i + 1][upper], L[i][upper]) + R[i + 1][lower];*/
				}
			}
		}
		/*
		for (int index = 0; index < Col; index++) {
			if (Previous_LLR.at(index) *(R._matrix[0][index] + L._matrix[0][index]) < 0) FlipRecord.at(index)++;
		}

		for (int index = 0; index < Col; index++) {
			Previous_LLR.at(index) = R._matrix[0][index] + L._matrix[0][index];
		}*/

		/*
		for (int index = 0; index < Col; index++) {
			if (Previous_LLR.at(index)*(R._matrix[0][index]) < 0) {
				if (index >= Col/2) {
					R._matrix[1][index - Col / 2] = min_sum(R._matrix[0][index] + Previous_LLR.at(index) + R._matrix[1][index], R._matrix[0][index - Col / 2]);
					R._matrix[1][index] = R._matrix[0][index] + Previous_LLR.at(index) + min_sum(R._matrix[0][index - Col / 2], R._matrix[1][index - Col / 2]);
				}
				else {
					R._matrix[1][index] = min_sum(R._matrix[0][index + Col / 2] + R._matrix[1][index + Col / 2], Previous_LLR.at(index) + R._matrix[0][index]);
					R._matrix[1][index + Col / 2] = R._matrix[0][index + Col / 2] + min_sum(R._matrix[0][index] + Previous_LLR.at(index), R._matrix[1][index]);
				}
			}
		}

		for (int index = 0; index < Col; index++) {
			Previous_LLR.at(index) = R._matrix[0][index];
		}*/ //FlipRecord


		if ((k == BP_Decode_Iter - 2) && Polar_with_Astar) {
			for (int i = 0; i < Col; i++)
			{
				decoding_info.rx_signal_seq.at(i) = (R._matrix[0][i] + L._matrix[0][i]);
				//if (FlipRecord.at(i) != 0) decoding_info.rx_signal_seq.at(i) /= FlipRecord.at(i);
				//cout << decoding_info.rx_signal_seq.at(i) << ", ";
			}
			//cout << endl << endl;
			A_star_PC_out_CBC_OSC_Verified(G, decoding_info);
			/*
			for (int i = 0; i < Col; i++)
			{
				if ((decoding_info.estimated_codeword.at(i) - 0.5)*L._matrix[0][i] > 0) {
					L._matrix[0][i] *= (-1);
				}
			}*/
			for (int i = 0; i < Col; i++) //*******************************************************   // N/2是對應coderate為1/2
			{
				/*
				if ((decoding_info.estimated_codeword.at(i) - 0.5) * decoding_info.rx_signal_seq.at(i) > 0) {
					R._matrix[0][i] *= (-1);
					L._matrix[0][i] *= (-1);
					//cout << "!" << endl;
				}*/

				if (decoding_info.estimated_codeword.at(i)) {
					//R._matrix[0 + 1][i] = -1000;
					L._matrix[0][i] = -1000;
				}
				else {
					//R._matrix[0 + 1][i] = 1000;
					L._matrix[0][i] = 1000;
				}

			}
		}
		//}
	}
	//cout << "!" << endl;
	//int message_output[Row];  //   ->192


	for (int i = 0; i < Row; i++) //*******************************************************   // N/2是對應coderate為1/2
	{
		//cout << L._matrix[log_Col][decoding_info.non_frozen[i]] << "," << R._matrix[log_Col][decoding_info.non_frozen[i]] << endl;
		//cout << L[logN][non_frozen[i]] + R[logN][non_frozen[i]] << " ";
		//cout << L._matrix[log_Col][decoding_info.non_frozen[i]] + R._matrix[log_Col][decoding_info.non_frozen[i]] << ",";

		if (L._matrix[log_Col][decoding_info.non_frozen[i]] + R._matrix[log_Col][decoding_info.non_frozen[i]] >= 0)
		{
			decoding_info.estimated_codeword[i] = 0; // 結果
		}
		else
		{
			decoding_info.estimated_codeword[i] = 1;
		}
		//cout << message_output[i]<<endl;
		/*Yin
		if (message[i] != message_output[i])
		{
			error_bits++;
		}
		*/
		//cout << error_bits;
		//cout << (int)decoding_info.estimated_codeword[i] << ",";

	}
	//cout << endl << endl;
}

void Polar_code_BP_Decoder_Astar_Decoding(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {

	int Row = G.Row_number, Col = G.Col_number;
	int log_Col = log2(Col);

	vector<int>u(Col, 0);
	vector<int>x(Col, 0);

	MATRIX <double>R;
	R.Building_Empty_Matrix(log_Col + 1, Col);
	MATRIX <double>L;
	L.Building_Empty_Matrix(log_Col + 1, Col);

	// Data 進入Polar Decoder
	for (int i = 0; i < Col; i++) //*******************************************************
	{
		L._matrix[0][i] = 2 * decoding_info.rx_signal_seq[i] / decoding_info.var;
	}

	//Decoding Start
	for (int i = 0; i < Col; i++)
	{
		if (decoding_info.frozen[i] == 0)
		{
			R._matrix[log_Col][i] = DBL_MAX;
		}//	cout << R[logN][i] << endl;
	}

	int upper, lower;
	short int reg = 0, triger = 1;

	for (int k = 0; k < BP_Decode_Iter; k++)
	{
		//cout << "Iteration:" << k << endl;
		for (int i = 0; i < log_Col; i++)
		{

			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					L._matrix[i + 1][upper] = min_sum(L._matrix[i][lower] + R._matrix[i + 1][lower], L._matrix[i][upper]);
					L._matrix[i + 1][lower] = L._matrix[i][lower] + min_sum(L._matrix[i][upper], R._matrix[i + 1][upper]);
				}
			}
		}

		for (int i = 0; i < log_Col; i++)
		{
			for (int i = log_Col; i >= 0; i--)
			{

				int GG = pow(2, i);
				for (int g = 0; g < GG; g++)
				{
					int p = pow(2, i + 1);

					for (int j = 0; j < Col / p; j++)
					{
						upper = g * (Col / GG) + j;
						lower = g * (Col / GG) + (Col / p) + j;

						R._matrix[i][upper] = min_sum(R._matrix[i + 1][lower] + L._matrix[i][lower], R._matrix[i + 1][upper]);
						R._matrix[i][lower] = min_sum(R._matrix[i + 1][upper], L._matrix[i][upper]) + R._matrix[i + 1][lower];
					}
				}
			}
		}
	}

	vector<double> Temp_R(Col, 0);
	for (int i = 0; i < Row; i++) //*******************************************************   // N/2是對應coderate為1/2
	{
		
		if (L._matrix[log_Col][decoding_info.non_frozen[i]] + R._matrix[log_Col][decoding_info.non_frozen[i]] >= 0)
		{
			decoding_info.estimated_codeword[i] = 0; // 結果
		}
		else
		{
			decoding_info.estimated_codeword[i] = 1;
		}
	}

	
}

void Inner_Polar_BP_Decoder_Outer_Astar_Decoder(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	int Row = G.Row_number, Col = G.Col_number;
	int log_Col = log2(Col);
	vector<int>u(Col, 0);
	vector<int>x(Col, 0);
	vector<double> R_temp = (decoding_info.rx_signal_seq);
	//cout << decoding_info.rx_signal_seq.at(0) << decoding_info.rx_signal_seq.at(1) << decoding_info.rx_signal_seq.at(2) << endl;
	//cout << R_temp.at(0) << R_temp.at(1) << R_temp.at(2) << endl;

	MATRIX <double>R;
	R.Building_Empty_Matrix(log_Col + 1, Col);
	MATRIX <double>L;
	L.Building_Empty_Matrix(log_Col + 1, Col);

	// Data 進入Polar Decoder
	for (int i = 0; i < Col; i++) //*******************************************************
	{
		L._matrix[0][i] = 2 * decoding_info.rx_signal_seq[i] / decoding_info.var;
	}
	//Decoding Start
	for (int i = 0; i < Col; i++)
	{
		if (decoding_info.frozen[i] == 0)
		{
			R._matrix[log_Col][i] = DBL_MAX;
		}
	}

	int upper, lower;

	for (int k = 0; k < BP_Decode_Iter; k++)
	{
		for (int i = 0; i < log_Col; i++)
		{
			int GG = pow(2, i);
			for (int g = 0; g < GG; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < Col / p; j++)
				{
					upper = g * (Col / GG) + j;
					lower = g * (Col / GG) + (Col / p) + j;

					L._matrix[i + 1][upper] = min_sum(L._matrix[i][lower] + R._matrix[i + 1][lower], L._matrix[i][upper]);
					L._matrix[i + 1][lower] = L._matrix[i][lower] + min_sum(L._matrix[i][upper], R._matrix[i + 1][upper]);
				}
			}
		}

		for (int i = 0; i < log_Col; i++)
		{
			for (int i = log_Col; i >= 0; i--)
			{

				int GG = pow(2, i);
				for (int g = 0; g < GG; g++)
				{
					int p = pow(2, i + 1);

					for (int j = 0; j < Col / p; j++)
					{
						upper = g * (Col / GG) + j;
						lower = g * (Col / GG) + (Col / p) + j;

						R._matrix[i][upper] = min_sum(R._matrix[i + 1][lower] + L._matrix[i][lower], R._matrix[i + 1][upper]);
						R._matrix[i][lower] = min_sum(R._matrix[i + 1][upper], L._matrix[i][upper]) + R._matrix[i + 1][lower];
					}
				}
			}
		}
	}
	double AverageLLR = 0;

	for (int i = 0; i < Row; i++) //*******************************************************   // N/2是對應coderate為1/2
	{
		decoding_info.rx_signal_seq.at(i) = (L._matrix[log_Col][decoding_info.non_frozen[i]] + R._matrix[log_Col][decoding_info.non_frozen[i]])*decoding_info.var / 2;
		//decoding_info.rx_signal_seq.at(i) *= decoding_info.Channel_Parameter.at(i);
		//cout << decoding_info.Channel_Parameter.at(i) << endl;
		AverageLLR += abs(decoding_info.rx_signal_seq.at(i)) * 2 / decoding_info.var;
	}

	// Test Stopping Criterion
	/*
	short int reg = 0, triger = 1;
	int breaker = 1;
	for (int i = 0; i < Col; i++)
	{
		if (decoding_info.rx_signal_seq.at(i) >= 0)
		{
			u[i] = 0;
		}
		else
		{
			u[i] = 1;
		}

		if (R._matrix[0][i] >= 0)
		{
			x[i] = 0;
		}
		else
		{
			x[i] = 1;
		}
	}

	for (int i = 0; i < Col; i++)
	{
		for (int j = 0; j < Row; j++)
		{
			reg ^= (G._matrix[j][i] & u[j]);
			//	cout << RM[j][i];
		}
		//cout << reg%2;
		if (reg % 2 != x[i])
		{
			breaker = 0;
			break;
		}
		reg = 0;
	}
	if (breaker == 1) {
		cout << "Successed!" << endl;
	}
	else {
		//R_temp
		decoding_info.rx_signal_seq.clear();
		decoding_info.rx_signal_seq = R_temp;
		cout << "Failed!" << endl;
	}

	
	// Test Stopping Criterion End
	*/

	//system("pause");
	/*
	for (int i = 0; i < Row; i++)
	{
		if (decoding_info.Channel_Parameter.at(i) == 2) cout << "*";
		cout << setprecision(3) << decoding_info.rx_signal_seq.at(i)<<" ";
	}
	system("pause");*/
	
	decoding_info.rx_signal_seq.resize(Row);
	//AverageLLR /= Row;
	//decoding_info.Ave_LLR += AverageLLR;
	//for (int i = 0; i < 20; i++) cout << setprecision(2) << decoding_info.rx_signal_seq.at(i) << " ";
	//cout << endl;
}


double min_sum(double a, double b)
{
	if (a*b >= 0)
		return min(abs(a), abs(b));
	else
		return (-1)*min(abs(a), abs(b));
}

void Majority_Soft_decoding(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	//cout << decoding_info.var << endl;
	vector<double> yx(decoding_info.rx_signal_seq.size(), 0);
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		yx.at(i) = tanh(2 * decoding_info.rx_signal_seq.at(i) / decoding_info.var);
		//cout << yx.at(i) << " ";
		//yx.at(i) = decoding_info.rx_signal_seq.at(i);
	}
	//cout << endl;

	if (decoding_info.Code_Number == RM_32_6) {
		
		int r = 1, m = 5;
		// r = 1
		vector<__int8> S(16, 0);
		vector<double> Final_Rx(6, 0);
		double W;
		// r = 1, a1 
		__int8 row = 5, S_2 = pow(2, 0), i_2;
		S = { 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W+= (yx.at(S.at(i))*yx.at(i_2));
		}
		Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		// r = 1, a2 
		--row;
		S_2 = pow(2, 1);
		S = { 0,1,4,5,8,9,12,13, 16,17,20,21,24,25,28,29 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += (yx.at(S.at(i))*yx.at(i_2));
		}
		Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		// r = 1, a3 
		--row;
		S_2 = pow(2, 2);
		S = { 0,1,2,3,8,9,10,11, 16,17,18,19,24,25,26,27 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += (yx.at(S.at(i))*yx.at(i_2));
		}
		Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		// r = 1, a4 
		--row;
		S_2 = pow(2, 3);
		S = { 0,1,2,3,4,5,6,7, 16,17,18,19,20,21,22,23 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += (yx.at(S.at(i))*yx.at(i_2));
		}
		Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		// r = 1, a5 
		--row;
		S_2 = pow(2, 4);
		S = { 0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += (yx.at(S.at(i))*yx.at(i_2));
		}
		Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		// r = 0
		for (int i = 0; i < 32; ++i) {
			for (int j = 1; j < 6; ++j) {
				if (G._matrix[j][i] == 1 && decoding_info.estimated_codeword.at(j) == 1) yx.at(i) *= (-1);
				//if (G._matrix[j][i] == 1 && decoding_info.estimated_codeword.at(j) == 1)  yx.at(i) *= Final_Rx.at(j);
				//cout << Final_Rx.at(j) << endl;
				//system("pause");
			}
		}
		--row; 
		W = 0;
		for (int i = 0; i < 32; ++i) {
			W += yx.at(i);
		}
		//cout << W << " " << endl;
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
	//	cout << "END of Majority" << endl;
	}

	else if (decoding_info.Code_Number == RM_16_11) {  // RM ( 16,11 )
		int r = 2, m = 4;
		__int8 row = 10;
		//vector<double> Final_Rx(11, 0);
		// r = 2
		vector<__int8> S(4, 0);
		double W,W1,W2,W3,W4;
		// r = 2, a12
		W1 = 2 * atan(yx.at(0) * yx.at(1) * yx.at(2) * yx.at(3));
		W2 = 2 * atan(yx.at(4) * yx.at(5) * yx.at(6) * yx.at(7));
		W3 = 2 * atan(yx.at(8) * yx.at(9) * yx.at(10) * yx.at(11));
		W4 = 2 * atan(yx.at(12) * yx.at(13) * yx.at(14) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		//Final_Rx.at(row) = W;
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 2, a13
		--row;
		W1 = 2 * atan(yx.at(0) * yx.at(1) * yx.at(4) * yx.at(5));
		W2 = 2 * atan(yx.at(2) * yx.at(3) * yx.at(6) * yx.at(7));
		W3 = 2 * atan(yx.at(8) * yx.at(9) * yx.at(12) * yx.at(13));
		W4 = 2 * atan(yx.at(10) * yx.at(11) * yx.at(14) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 2, a23
		--row;
		W1 = 2 * atan(yx.at(0) * yx.at(2) * yx.at(4) * yx.at(6));
		W2 = 2 * atan(yx.at(1) * yx.at(3) * yx.at(5) * yx.at(7));
		W3 = 2 * atan(yx.at(8) * yx.at(10) * yx.at(12) * yx.at(14));
		W4 = 2 * atan(yx.at(9) * yx.at(11) * yx.at(13) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 2, a14
		--row;
		W1 = 2 * atan(yx.at(0) * yx.at(1) * yx.at(8) * yx.at(9));
		W2 = 2 * atan(yx.at(2) * yx.at(3) * yx.at(10) * yx.at(11));
		W3 = 2 * atan(yx.at(4) * yx.at(5) * yx.at(12) * yx.at(13));
		W4 = 2 * atan(yx.at(6) * yx.at(7) * yx.at(14) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 2, a24
		--row;
		W1 = 2 * atan(yx.at(0) * yx.at(2) * yx.at(8) * yx.at(10));
		W2 = 2 * atan(yx.at(1) * yx.at(3) * yx.at(9) * yx.at(11));
		W3 = 2 * atan(yx.at(4) * yx.at(6) * yx.at(12) * yx.at(14));
		W4 = 2 * atan(yx.at(5) * yx.at(7) * yx.at(13) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 2, a34
		--row;
		W1 = 2 * atan(yx.at(0) * yx.at(4) * yx.at(8) * yx.at(12));
		W2 = 2 * atan(yx.at(1) * yx.at(5) * yx.at(9) * yx.at(13));
		W3 = 2 * atan(yx.at(2) * yx.at(6) * yx.at(10) * yx.at(14));
		W4 = 2 * atan(yx.at(3) * yx.at(7) * yx.at(11) * yx.at(15));
		W = (W1 + W2 + W3 + W4) / 4;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";

		// r = 2 -> r = 1
		for (int i = 0; i < 16; ++i) {
			for (int j = 5; j < 11; ++j) {
				if (G._matrix[j][i] == 1 && decoding_info.estimated_codeword.at(j) == 1) yx.at(i) *= (-1);
			}
		}

		// r = 1
		//vector<double> Final_Rx(6, 0);
		// r = 1, a1 
		S.resize(8);
		--row;
		int S_2 = pow(2, 0), i_2;
		S = { 0,2,4,6,8,10,12,14 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += 2 * atan((yx.at(S.at(i))*yx.at(i_2)));
		}
		//Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		W /= S.size();
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 1, a2 
		--row;
		S_2 = pow(2, 1);
		S = { 0,1,4,5,8,9,12,13 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += 2 * atan((yx.at(S.at(i))*yx.at(i_2)));
		}
		//Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		W /= S.size();
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 1, a3 
		--row;
		S_2 = pow(2, 2);
		S = { 0,1,2,3,8,9,10,11 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += 2 * atan((yx.at(S.at(i))*yx.at(i_2)));
		}
		//Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		W /= S.size();
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		// r = 1, a4 
		--row;
		S_2 = pow(2, 3);
		S = { 0,1,2,3,4,5,6,7 };
		W = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			W += 2 * atan((yx.at(S.at(i))*yx.at(i_2)));
		}
		//Final_Rx.at(row) = W / S.size();
		//cout << W << " ";
		W /= S.size();
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << " ";
		
		// r = 1 -> r = 0
		for (int i = 0; i < 16; ++i) {
			for (int j = 1; j < 5; ++j) {
				if (G._matrix[j][i] == 1 && decoding_info.estimated_codeword.at(j) == 1) yx.at(i) *= (-1);
				//if (G._matrix[j][i] == 1 && decoding_info.estimated_codeword.at(j) == 1)  yx.at(i) *= Final_Rx.at(j);
				//cout << Final_Rx.at(j) << endl;
				//system("pause");
			}
		}
		// r = 0
		--row;
		W = 0;
		for (int i = 0; i < 16; ++i) {
			W += 2* atan(yx.at(i));
		}
		//cout << W << " " << endl;
		W /= 16;
		decoding_info.rx_signal_seq.at(row) = (W);
		if (W > 0) decoding_info.estimated_codeword.at(row) = 0;
		else decoding_info.estimated_codeword.at(row) = 1;
		//cout << setprecision(2) << W << endl;
	}

}


void Majority_Soft_decoding(MATRIX<__int8> &G, vector<double> &Rx, double var) {
	//cout << decoding_info.var << endl;
	vector<long double> yx(Rx.size(), 0);
	for (int i = 0; i < Rx.size(); ++i) {
		yx.at(i) = tanh(2 * Rx.at(i) / var);
	}

	// RM ( 16,11 )
	int r = 2, m = 4;
	__int8 row = 10;
	//vector<double> Final_Rx(11, 0);
	// r = 2
	vector<__int8> S(4, 0);
	long double W, W1, W2, W3, W4;
	// r = 2, a12
	W1 = RM_LLR_order2(yx.at(0), yx.at(1), yx.at(2), yx.at(3));
	W2 = RM_LLR_order2(yx.at(4), yx.at(5), yx.at(6), yx.at(7));
	W3 = RM_LLR_order2(yx.at(8), yx.at(9), yx.at(10), yx.at(11));
	W4 = RM_LLR_order2(yx.at(12), yx.at(13), yx.at(14), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));
	// r = 2, a13
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(1), yx.at(4), yx.at(5));
	W2 = RM_LLR_order2(yx.at(2), yx.at(3), yx.at(6), yx.at(7));
	W3 = RM_LLR_order2(yx.at(8), yx.at(9), yx.at(12), yx.at(13));
	W4 = RM_LLR_order2(yx.at(10), yx.at(11), yx.at(14), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));
	// r = 2, a23
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(2), yx.at(4), yx.at(6));
	W2 = RM_LLR_order2(yx.at(1), yx.at(3), yx.at(5), yx.at(7));
	W3 = RM_LLR_order2(yx.at(8), yx.at(10), yx.at(12), yx.at(14));
	W4 = RM_LLR_order2(yx.at(9), yx.at(11), yx.at(13), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));
	// r = 2, a14
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(1), yx.at(8), yx.at(9));
	W2 = RM_LLR_order2(yx.at(2), yx.at(3), yx.at(10), yx.at(11));
	W3 = RM_LLR_order2(yx.at(4), yx.at(5), yx.at(12), yx.at(13));
	W4 = RM_LLR_order2(yx.at(6), yx.at(7), yx.at(14), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));
	// r = 2, a24
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(2), yx.at(8), yx.at(10));
	W2 = RM_LLR_order2(yx.at(1), yx.at(3), yx.at(9), yx.at(11));
	W3 = RM_LLR_order2(yx.at(4), yx.at(6), yx.at(12), yx.at(14));
	W4 = RM_LLR_order2(yx.at(5), yx.at(7), yx.at(13), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));
	// r = 2, a34
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(4), yx.at(8), yx.at(12));
	W2 = RM_LLR_order2(yx.at(1), yx.at(5), yx.at(9), yx.at(13));
	W3 = RM_LLR_order2(yx.at(2), yx.at(6), yx.at(10), yx.at(14));
	W4 = RM_LLR_order2(yx.at(3), yx.at(7), yx.at(11), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));

	// r = 2 -> r = 1
	for (int i = 0; i < 16; ++i) {
		for (int j = 5; j < 11; ++j) {
			if (G._matrix[j][i] == 1 && Rx.at(j) < 0) yx.at(i) *= (-1);
		}
	}

	// r = 1
	//vector<double> Final_Rx(6, 0);
	// r = 1, a1 
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(4), yx.at(8), yx.at(13));
	W2 = RM_LLR_order2(yx.at(1), yx.at(5), yx.at(9), yx.at(12));
	W3 = RM_LLR_order2(yx.at(2), yx.at(6), yx.at(10), yx.at(15));
	W4 = RM_LLR_order2(yx.at(3), yx.at(7), yx.at(11), yx.at(14));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));

	// r = 1, a2 
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(4), yx.at(8), yx.at(14));
	W2 = RM_LLR_order2(yx.at(1), yx.at(5), yx.at(9), yx.at(15));
	W3 = RM_LLR_order2(yx.at(2), yx.at(6), yx.at(10), yx.at(12));
	W4 = RM_LLR_order2(yx.at(3), yx.at(7), yx.at(11), yx.at(13));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));

	// r = 1, a3 
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(5), yx.at(10), yx.at(11));
	W2 = RM_LLR_order2(yx.at(1), yx.at(4), yx.at(8), yx.at(9));
	W3 = RM_LLR_order2(yx.at(2), yx.at(7), yx.at(14), yx.at(15));
	W4 = RM_LLR_order2(yx.at(3), yx.at(6), yx.at(12), yx.at(13));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));

	// r = 1, a4 
	--row;
	W1 = RM_LLR_order2(yx.at(0), yx.at(4), yx.at(5), yx.at(9));
	W2 = RM_LLR_order2(yx.at(1), yx.at(6), yx.at(7), yx.at(8));
	W3 = RM_LLR_order2(yx.at(2), yx.at(11), yx.at(12), yx.at(13));
	W4 = RM_LLR_order2(yx.at(3), yx.at(10), yx.at(14), yx.at(15));
	W = (W1 + W2 + W3 + W4) / 4;
	Rx.at(row) = ((W));

	// r = 1 -> r = 0
	for (int i = 0; i < 16; ++i) {
		for (int j = 1; j < 5; ++j) {
			if (G._matrix[j][i] == 1 && Rx.at(j) < 0) yx.at(i) *= (-1);
		}
	}
	// r = 0
	--row;
	W = 0;
	for (int i = 0; i < 16; ++i) {
		W += atanh(yx.at(i));
	}
	/*
	if (W > 1000 || W < -1000) {
		cout << W << endl;
		for (int i = 0; i < yx.size(); ++i) {
			cout << yx.at(i) << ", ";
		}
		cout << endl;
		for (int i = 0; i < yx.size(); ++i) {
			cout << atanh(yx.at(i)) << ", ";
		}
		cout << endl;

		system("pause");
	}*/
	//cout << W << " " << endl;
	W /= 16;
	Rx.at(row) = ((W));
	//cout << setprecision(2) << W << endl;

}

void Majority_Soft_decoding_ver2(MATRIX<__int8> &G, MATRIX<__int8> &G_, vector<double> &Rx, double var) {
	int method = 2;   //     1: Semi-soft method
	                  //     2: All soft method
	// For (64,44) RM code

	vector<long double> yx(Rx.size(), 0);
	//double W;
	for (int i = 0; i < Rx.size(); ++i) {
		yx.at(i) = tanh(Rx.at(i) / var);   // Result : LLR/2
	}
	// RM ( 64,42 )
	int r = 3, m = 6, n = pow(2, m);        // initialize r, m, n (length of codeword)
	int k = 0;                              // initialize length of message 
	vector<__int8> order_length_record;
	for (int i = r; i > -1; --i) {
		int temp = combination(m, r - i);
		k += temp;
		order_length_record.push_back(temp);
	}
	vector<long double>message_LLR(k, 0);

	// r = 3
	int k_temp = k - order_length_record.at(3);
	for (k_temp; k_temp < k; k_temp++) {
		for (int Set = 0; Set < 8; Set++) {
			double temp = 1;
			for (int index = 0; index < 8; index++) {
				temp *= yx.at(G_._matrix[k_temp][index + Set * 8]);
			}
			//cout <<  setw(4) << temp << " ";
			message_LLR.at(k_temp) += 2 * atanh(temp);
		}
		//cout << endl;
		//message_LLR.at(k_temp) /= 8;
	}
	//cout << endl;
	// r = 2
	k_temp = k - order_length_record.at(3) - order_length_record.at(2);
	int k_upper_B = k - order_length_record.at(3);
	vector<__int8> Col_Position(4, 0);

	if (method == 1) {
		for (int row = k - order_length_record.at(3); row < k; row++) {
			for (int col = 0; col < Rx.size(); col++) {
				if (G._matrix[row][col] == 1 && message_LLR.at(row) < 0) yx.at(col) *= (-1);
			}
		}
	}
	for (k_temp; k_temp < k_upper_B; k_temp++) {
		for (int Set = 0; Set < 16; Set++) {
			double temp = 1;
			for (int index = 0; index < 4; index++) {
				Col_Position.at(index) = G_._matrix[k_temp][index + Set * 4];
				temp *= yx.at(Col_Position.at(index));
			}
			if (method == 2) {
				// Test r = 3 interference
				vector<__int8> Candidates;
				for (int row = k - order_length_record.at(3); row < k; ++row) {
					if (G._matrix[row][Col_Position.at(0)] ^ G._matrix[row][Col_Position.at(1)] ^ G._matrix[row][Col_Position.at(2)] ^ G._matrix[row][Col_Position.at(3)] == 1)
						Candidates.push_back(row);
				}
				for (int i = 0; i < Candidates.size(); ++i) {
					temp *= tanh(message_LLR.at(Candidates.at(i)) / 2);
				}
			}
			message_LLR.at(k_temp) += 2 * atanh(temp);
		}
		//message_LLR.at(k_temp) /= 16;
	}

	// r = 1
	k_temp = 1;
    k_upper_B = k - order_length_record.at(3) - order_length_record.at(2);
	Col_Position.resize(2);

	if (method == 1) {
		for (int row = k - order_length_record.at(3) - order_length_record.at(2); row < k - order_length_record.at(3); row++) {
			for (int col = 0; col < Rx.size(); col++) {
				if (G._matrix[row][col] == 1 && message_LLR.at(row) < 0) yx.at(col) *= (-1);
			}
		}
	}
	for (k_temp; k_temp < k_upper_B; k_temp++) {
		for (int Set = 0; Set < 32; Set++) {
			double temp = 1;
			for (int index = 0; index < 2; index++) {
				Col_Position.at(index) = G_._matrix[k_temp][index + Set * 2];
				temp *= yx.at(Col_Position.at(index));
			}
			if (method == 2) {
				// Test r = 3 , r = 2 interference
				vector<__int8> Candidates;
				for (int row = k - order_length_record.at(3) - order_length_record.at(2); row < k; ++row) {
					if (G._matrix[row][Col_Position.at(0)] ^ G._matrix[row][Col_Position.at(1)] == 1)
						Candidates.push_back(row);
				}
				for (int i = 0; i < Candidates.size(); ++i) {
					temp *= tanh(message_LLR.at(Candidates.at(i)) / 2);
				}
			}
			message_LLR.at(k_temp) += 2 * atanh(temp);
		}
		//message_LLR.at(k_temp) /= 32;
	}


	// r = 0

	if (method == 1) {
		for (int row = 1; row < k - order_length_record.at(3) - order_length_record.at(2); row++) {
			for (int col = 0; col < Rx.size(); col++) {
				if (G._matrix[row][col] == 1 && message_LLR.at(row) < 0) yx.at(col) *= (-1);
			}
		}
	}
	for (int Set = 0; Set < 64; Set++) {
		double temp = yx.at(Set);
		if (method == 2) {
			// Test r = 3 , r = 2, r = 1 interference
			vector<__int8> Candidates;
			for (int row = 1; row < k; ++row) {
				if (G._matrix[row][Set] == 1)
					Candidates.push_back(row);
			}
			for (int i = 0; i < Candidates.size(); ++i) {
				temp *= tanh(message_LLR.at(Candidates.at(i)) / 2);
			}
		}
		message_LLR.at(0) += 2 * atanh(temp);
	}

	
	message_LLR.at(0) /= (64);
	
	for (int i = 1; i < k - order_length_record.at(3) - order_length_record.at(2); ++i) {
		message_LLR.at(i) /= (32);
	}
	for (int i = k - order_length_record.at(3) - order_length_record.at(2); i < k - order_length_record.at(3); ++i) {
		message_LLR.at(i) /= (16);
	}
	for (int i = k - order_length_record.at(3); i < k ; ++i) {
		message_LLR.at(i) /= 8;
	}
	
	// End Decoding
	for (int i = 0; i < message_LLR.size(); ++i) {
		Rx.at(i) = message_LLR.at(i);
	}
}

void Majority_Soft_decoding_ver3(MATRIX<__int8> &G, MATRIX<__int8> &G_, vector<double> &Rx, double var) {
	// For (64,44) RM code

	vector<long double> yx(Rx.size(), 0);

	// Pre-processing
	for (int i = 0; i < Rx.size(); ++i) {
		yx.at(i) = tanh(Rx.at(i) / var);   // Result : LLR/2
	}
	// RM ( 64,42 )
	int r = 3, m = 6, n = pow(2, m);        // initialize r, m, n (length of codeword)
	int k = 0;                              // initialize length of message 
	vector<__int8> order_length_record;
	for (int i = r; i > -1; --i) {
		int temp = combination(m, r - i);
		k += temp;
		order_length_record.push_back(temp);
	}
	vector<long double>message_LLR(k, 0);

	// Decoding Start
	long double temp;
	for (int row = k - 1; row > -1; --row) {
		//cout << row << ": " << endl;
		temp = 1;
		for (int set_traversal = 0; set_traversal < G_._matrix.at(row).size(); ++set_traversal) {
			//cout << "(" << temp << "," << message_LLR.at(row) << ")";
			if (G_._matrix[row][set_traversal] == CHAR_MAX) {
				message_LLR.at(row) += (2 * atanh(temp));
				temp = 1;
			}
			else if (G_._matrix[row][set_traversal] < 0) {
				temp *= tanh(message_LLR.at(G_._matrix[row][set_traversal] * (-1)) / 2);
				//cout << (message_LLR.at(G_._matrix[row][set_traversal] * (-1)) / 2) << "* ";
			}
			else temp *= yx.at(G_._matrix[row][set_traversal]);
		}
		//cout << endl;
		//cout << "R: " << message_LLR.at(row) << endl;
	}
//	cout << endl << endl;

	/*
	message_LLR.at(0) /= (4);

	for (int i = 1; i < k - order_length_record.at(3) - order_length_record.at(2); ++i) {
		message_LLR.at(i) /= (2);
	}
	for (int i = k - order_length_record.at(3) - order_length_record.at(2); i < k - order_length_record.at(3); ++i) {
		//message_LLR.at(i) /= (16);
	}
	*/
	// End Decoding
	for (int i = 0; i < message_LLR.size(); ++i) {
		Rx.at(i) = message_LLR.at(i);
		//cout << message_LLR.at(i) << " ";
	}
	//system("pause");
}


void Majority_Hard_decoding(MATRIX<__int8> &G, DECODING_INFO &decoding_info) {
	//cout << decoding_info.var << endl;
	vector<__int8> yx(decoding_info.rx_signal_seq.size(), 0);
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		if (decoding_info.rx_signal_seq.at(i) > 0) yx.at(i) = 0;
		else yx.at(i) = 1;
	}

	if (decoding_info.Code_Number == RM_32_6) {

		int r = 1, m = 5;
		// r = 1
		vector<__int8> S(16, 0);
		//vector<__int8> Final_Rx(6, 0);
		__int8 Num_0, Num_1;
		// r = 1, a1 
		__int8 row = 5, S_2 = pow(2, 0), i_2;
		Num_0 = 0;
		Num_1 = 0;
		S = { 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 };
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			if(yx.at(S.at(i))^yx.at(i_2)==1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " ";
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		// r = 1, a2 
		--row;
		S_2 = pow(2, 1);
		S = { 0,1,4,5,8,9,12,13, 16,17,20,21,24,25,28,29 };
		Num_0 = 0;
		Num_1 = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			if (yx.at(S.at(i)) ^ yx.at(i_2) == 1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " ";
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		// r = 1, a3 
		--row;
		S_2 = pow(2, 2);
		S = { 0,1,2,3,8,9,10,11, 16,17,18,19,24,25,26,27 };
		Num_0 = 0;
		Num_1 = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			if (yx.at(S.at(i)) ^ yx.at(i_2) == 1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " ";
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		// r = 1, a4 
		--row;
		S_2 = pow(2, 3);
		S = { 0,1,2,3,4,5,6,7, 16,17,18,19,20,21,22,23 };
		Num_0 = 0;
		Num_1 = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			if (yx.at(S.at(i)) ^ yx.at(i_2) == 1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " ";
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		// r = 1, a5 
		--row;
		S_2 = pow(2, 4);
		S = { 0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15 };
		Num_0 = 0;
		Num_1 = 0;
		for (int i = 0; i < S.size(); ++i) {
			i_2 = S.at(i) + S_2;
			if (yx.at(S.at(i)) ^ yx.at(i_2) == 1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " ";
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		// r = 0
		for (int i = 0; i < 32; ++i) {
			for (int j = 1; j < 6; ++j) {
				if (decoding_info.estimated_codeword.at(j) == 1 && G._matrix[j][i] == 1) yx.at(i) ^= (1);
			}
		}
		--row;
		Num_0 = 0;
		Num_1 = 0;
		for (int i = 0; i < 32; ++i) {
			if (yx.at(i) == 1) ++Num_1;
			else ++Num_0;
		}
		//cout << W << " " << endl;
		if (Num_1 > Num_0) decoding_info.estimated_codeword.at(row) = 1;
		else decoding_info.estimated_codeword.at(row) = 0;
		//	cout << "END of Majority" << endl;
	}
}

void Majority_Hard_decoding(MATRIX<__int8> &G, MATRIX<__int8> &G_,DECODING_INFO &decoding_info) {
	//cout << decoding_info.var << endl;
	vector<__int8> yx(decoding_info.rx_signal_seq.size(), 0);
	for (int i = 0; i < decoding_info.rx_signal_seq.size(); ++i) {
		if (decoding_info.rx_signal_seq.at(i) > 0) yx.at(i) = 0;
		else yx.at(i) = 1;
	}
	if (decoding_info.Code_Number == RM_64_42) {
		int r = 3, m = 6, n = pow(2, m);        // initialize r, m, n (length of codeword)
		int k = 0;                              // initialize length of message 
		vector<__int8> order_length_record;
		for (int i = r; i > -1; --i) {
			int temp = combination(m, r - i);
			k += temp;
			order_length_record.push_back(temp);
		}

		// r = 3
		int k_temp = k - order_length_record.at(3);
		for (k_temp; k_temp < k; k_temp++) {
			for (int Set = 0; Set < 8; Set++) {
				__int8 temp = 0;
				for (int index = 0; index < 8; index++) {
					temp ^= yx.at(G_._matrix[k_temp][index + Set * 8]);
				}
				decoding_info.estimated_codeword.at(k_temp) += temp;
			}
			if (decoding_info.estimated_codeword.at(k_temp) > 4)decoding_info.estimated_codeword.at(k_temp) = 1;
			else decoding_info.estimated_codeword.at(k_temp) = 0;
		}
		// renew the above order
		for (int i = 0; i < n; ++i) {
			for (int j = k - order_length_record.at(3); j < k; ++j) {
				if (decoding_info.estimated_codeword.at(j) == 1 && G._matrix[j][i] == 1) yx.at(i) ^= (1);
			}
		}

		// r = 2
		k_temp = k - order_length_record.at(3) - order_length_record.at(2);
		int k_upper_B = k - order_length_record.at(3);
		for (k_temp; k_temp < k_upper_B; k_temp++) {
			for (int Set = 0; Set < 16; Set++) {
				__int8 temp = 0;
				for (int index = 0; index < 4; index++) {
					temp ^= yx.at(G_._matrix[k_temp][index + Set * 4]);
				}
				decoding_info.estimated_codeword.at(k_temp) += temp;
			}
			if (decoding_info.estimated_codeword.at(k_temp) > 8)decoding_info.estimated_codeword.at(k_temp) = 1;
			else decoding_info.estimated_codeword.at(k_temp) = 0;
		}
		// renew the above order
		for (int i = 0; i < n; ++i) {
			for (int j = k - order_length_record.at(3) - order_length_record.at(2); j < k_upper_B; ++j) {
				if (decoding_info.estimated_codeword.at(j) == 1 && G._matrix[j][i] == 1) yx.at(i) ^= (1);
			}
		}

		// r = 1
		k_temp = 1;
		k_upper_B = k - order_length_record.at(3) - order_length_record.at(2);
		for (k_temp; k_temp < k_upper_B; k_temp++) {
			for (int Set = 0; Set < 32; Set++) {
				__int8 temp = 0;
				for (int index = 0; index < 2; index++) {
					temp ^= yx.at(G_._matrix[k_temp][index + Set * 2]);
				}
				decoding_info.estimated_codeword.at(k_temp) += temp;
			}
			if (decoding_info.estimated_codeword.at(k_temp) > 16)decoding_info.estimated_codeword.at(k_temp) = 1;
			else decoding_info.estimated_codeword.at(k_temp) = 0;
		}
		// renew the above order
		for (int i = 0; i < n; ++i) {
			for (int j = 1; j < k_upper_B; ++j) {
				if (decoding_info.estimated_codeword.at(j) == 1 && G._matrix[j][i] == 1) yx.at(i) ^= (1);
			}
		}

		// r = 0
		__int8 temp = 0;
		for (k_temp = 0; k_temp < 64; k_temp++) {
			temp += yx.at(k_temp);
		}
		if (temp > 32) decoding_info.estimated_codeword.at(0) = 1;
		else decoding_info.estimated_codeword.at(0) = 0;
	}

}
double RM_LLR_order2(double L1, double L2, double L3, double L4) {
	//double numerator = 1 + exp(L1 + L2) + exp(L1 + L3) + exp(L1 + L4) + exp(L2 + L3) + exp(L2 + L4) + exp(L3 + L4) + exp(L1 + L2 + L3 + L4);
	//double denominator = exp(L1) + exp(L2) + exp(L3) + exp(L4) + exp(L1 + L2 + L3) + exp(L1 + L2 + L4) + exp(L2 + L3 + L4) + exp(L1 + L3 + L4);
	//return log(numerator / denominator);
	return 2*(atanh(L1*L2*L3*L4));
}

double RM_LLR_order1(double L1, double L2) {
	//return log((1 + exp(L1 + L2)) / (exp(L1) + exp(L2)));
	return (2 * atanh(L1*L2));
}

