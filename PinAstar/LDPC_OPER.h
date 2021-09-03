#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "DefineClass.h"
#include "DefineParameter.h"
#include "AstarDecode.h"
//#include "LinearBlockCodes.h"
using namespace std;

void ReadFile_H_Matrix_ver2(string name, MATRIX<__int8> &H);
void SPA(MATRIX<__int8> &H, DECODING_INFO &decoding_info);
void SPA_Partial(MATRIX<__int8> &H, DECODING_INFO &decoding_info);
void SPA(MATRIX<__int8> &H, vector<double> &Rx, double var);
void LDPC_VC_RBP(MATRIX<__int8> &H, DECODING_INFO &decoding_info);
void LDPC_Modified_RBP(MATRIX<__int8> &H, DECODING_INFO &decoding_info);

class LDPC_FUNC
{
private:
	//void ReadFile_H_Matrix(string, MATRIX<__int8> &);

public:
	//MATRIX<__int8> &H;
	int fullrank = 0;
	LDPC_FUNC();

	void Col_Exchange(MATRIX<__int8> &, int, int);
	void Row_Exchange(MATRIX<__int8> &, int, int);
	void Col_Addition(MATRIX<__int8> &, int col1, int col2); // value of col1 add on col2 
	void Row_Addition(MATRIX<__int8> &, int row1, int row2); // value of row1 add on row2  
	void Sort_Matrix_Col_Forward(MATRIX<__int8> &, vector<int> &permutation_seq); // Sort Matrix's column according to the permutation seqeunce forwardly
	void Sort_Sequence_Backward(vector<double> &Rx, vector<int> &permutation_seq); //Sort Rx's column according to the permutation seqeunce reversely
	void Sort_Sequence_Backward(vector<__int8> &Rx, vector<int> &permutation_seq);
	void Sort_Sequence_Forward(vector<double> &Rx, vector<int> &permutation_seq); //Sort Rx's column according to the permutation seqeunce reversely
	void Sort_Sequence_Forward(vector<__int8> &Rx, vector<int> &permutation_seq);

	vector<double> Seq_Add(vector<double> cc, vector<double> nn);
	void Show_Matrix(MATRIX<__int8> &Matrix, int row_min, int row_max, int col_min, int col_max);
	vector<double> PAM(vector<size_t> codew);
	double Variance(double SNR_value);
	vector<double> Noise_Generator(int num, double mean, double var);

	vector<size_t> Message_Generator(size_t);
	vector<size_t> Codeword_Generator(MATRIX<__int8> &Matrix, vector<size_t> &);
	vector<__int8> Codeword_Generator(MATRIX<__int8> &, vector<__int8> &);
	vector<size_t> Estimation(vector <double> LLR);
	int BER_Counting(vector<size_t> &, vector<size_t> &, size_t);

	void Print_Matrix_Col(MATRIX<__int8>, int);
	void Print_Matrix_Row(MATRIX<__int8>, int);
	void H_Test(MATRIX<__int8>);
	void G_H_Multiple_Test(MATRIX<__int8>, MATRIX<__int8>);

	void ReadFile_H_Matrix(string, MATRIX<__int8> &);
	void H_G_convertor_G_I_P(MATRIX<__int8>, vector<int> &, MATRIX<__int8> &);
	void H_G_convertor_G_I_P_ver2(MATRIX<__int8> H, vector<int> &permutation_seq, MATRIX<__int8> & G);
	void H_G_convertor_G_P_I(MATRIX<__int8>, vector<int> &, MATRIX<__int8> &);
	//friend Astar_OuterCode;
};


class Astar_OuterCode {

public:
	double var;
	int AstarCounter = 0;

	//vector < vector<__int8> > Message_Seq;
	vector < vector <__int8> > Message_Seq, CodeWord_Seq, Error_Seq;
	vector < vector <double> > Received_Seq;
	vector<__int8> LDPC_Message_Seq, Inner_CodeWord_Seq; //LDPC_CodeWord_Seq           //  上次進度:要把LDPC全部換成Inner !
	vector<double> Tx_Signal_Seq_series, Rx_Signal_Seq_series;          // series of A* codes
	int MessageLength,CodeLength,CodeAmount,LDPC_RowNumber,LDPC_ColNumber;

	Astar_OuterCode(int messagelength,int codelength,int AmountOfCodes,int LDPCRow,int LDPCCol);

	// General part
	void All_Generate_Message(void(*Func)(size_t , vector<__int8> &));
	void All_Systematic_Linear_Block_Code_Encoder(void(*Func)(MATRIX<__int8> &, vector<__int8> &, vector<__int8> &),MATRIX<__int8> &);
	// 1. LDPC part
	void Seq_Combine();
	void Seq_Seperate();
	void Astar_Initialization();
	// 2. Trellis part
	void ParityCheck_Encoder();
	void ProductCode_Seq_Combine();

	void ProductCode_Seq_Seperate();
	void Trellis_Decoding();
	void Trellis_Decoding(vector<int>& count);

	// 3.Reverse Trellis
	void ParityCheck_Encoder_Rev();
	void ProductCode_Seq_Combine_Rev();
	void ProductCode_Seq_Seperate_Rev();
};

// Binary-Tree for codeword candidates

class BinaryTree;

class TreeNode {

public:
	TreeNode *upchild;
	TreeNode *downchild;
	TreeNode *parent;
	//double Metric;
	double MetricTotal;
	int BinaryBit;
	int BinaryBitTotal;
	int depth;

	//double Total_Metric;

	TreeNode() : parent(NULL), upchild(NULL), downchild(NULL), MetricTotal(0), BinaryBit(-1), depth(0), BinaryBitTotal(0) { BinaryBit = DBL_MAX; };
	TreeNode(int node, int level) : parent(NULL), upchild(NULL), downchild(NULL), BinaryBit(-1), MetricTotal(0), BinaryBitTotal(0) {
		//Metric = number;
		depth = level;
		BinaryBit = node;
	}

	friend class BinaryTree;
};

class BinaryTree {

public:
	//int NodeCount = 0;
	//double Var;
	TreeNode *root;
	double BestMetric = DBL_MAX;
	TreeNode *BestNodePtr = NULL;
	vector <int> HardSequence;
	vector <double> RecSequence;
	vector<double> MetricTable_0;
	vector<double> MetricTable_1;
	vector<int> ResultHardSeq;
	vector<double> ResultRecSeq;

	BinaryTree(vector<int> HardSeq, vector<double> RecSeq);
	void BinaryTree_Construction(double var);
	void ConstructTree(TreeNode* ptr, int depth);
	void DeleteTree(TreeNode* ptr);
	void ConstructAstar();

};

int Cycle4_Check(MATRIX<__int8> &H);
int Cycle8_Check(MATRIX<__int8> &H);