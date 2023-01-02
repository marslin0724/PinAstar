#pragma once
#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <functional>   // std::bit_xor
#include <math.h>
#include "FileOperation.h"
#include "DefineClass.h"
#include "AStarDecode.h"
#include "POLAR_OPER.h"
#include "LDPC_OPER.h"
#include "TURBO_OPER.h"

#define Golay_24_12		1
#define QR_48_24		2
#define QR_80_40		3
#define QR_90_45		4
#define QR_98_49		5
#define	QR_104_52		6

#define BCH_128_64		7
#define BCH_128_99		8
#define BCH_256_45		9
#define BCH_256_99		10
#define BCH_256_131		11	

#define LDPC_96_48		12
#define Polar_120_40	13
#define ReadFromTxt		14
#define RandomCode		15
#define BPSK			16

#define BCH_128_78               17
#define QR_192_96                18
#define RM_16_11                 19
#define RM_32_6                  20
#define RM_64_42                 21
#define Hamming_Code_15_11       22
#define LDPC_16_12               23

//PoHan
#define double_RandomCode		24
#define Parity_Check			25
#define Parity_CheckII			26
#define Parity_CheckIII			27
#define Random_ParityCheck		28
#define ReedSolomonConcatenated	29

//Pin
#define BCH_128_36				30
#define RSGF2_64_36				31
#define RS_16_9_concat_8_4		32
#define HM_8_4					33
#define BCH_64_36_Rep			34
#define HM2_1					35
#define LDPC_16_8				36
#define RS_16_9_concat_8_4_interleaver 37
#define RSGF2_64_36_Rep			38
#define RS16_15					39
#define RS16_15_concate_HM8_4	40
#define	RS16_14					41
#define RS16_14_concate_HM8_4	42
#define RS16_13					43
#define RS16_13_concate_HM8_4	44
#define RS_16_9_concat_8_4_new	45
#define C6_5					46
#define RS_32_19_concat_6_5		47
#define RSconcate_192_96		48
#define RS16_9_concat_16_8		49
#define RS16_9_concat_24_12		50
#define RS16_8_concat_16_8		51
#define RS16_10_concat_16_8		52
#define RS16_11_concat_16_8		53
#define BCH_64_36				54
#define RS16_10_concat_8_4		55

// 這邊的LDPC的數字是codeword (n,k), 直接用SPA解(或後面接Astar)
#define LDPC_256_192                90
#define LDPC_128_64                 91
#define LDPC_192_96                 92
#define LDPC_256_128                93
#define LDPC_512_256                94
#define LDPC_224_128                95
#define LDPC_160_128                96
#define LDPC_192_128                97

#define Polar_Code                  100

// Concatenation Codes
#define Polar_256_192_AStar_192_128     101
#define LDPC_256_192_AStar_192_128      102
#define Hamming_255_187_Astar_187_128   103  // Inner: 17 (15,11) Hamming Code, Outer: RC (187,128) ; Total: (255,128)
#define RM_256_192_Astar_192_128        104  // Inner: 16 (16,11) Hamming Code, Outer: RC (176,128) ; Total: (256,128)
#define RM_256_168_Astar_168_128        105  // Inner:  4 (64,42) Hamming Code, Outer: RC (168,128) ; Total: (256,128)
#define Turbo_256_128_with_2LDPC        106  // p1: LDPC (192,128), p2: LDPC(192,128) 
#define Turbo_256_128_with_2LDPC_Punc   107  // p1: LDPC (192-x,128), p2: LDPC(192+x,128)  
#define Turbo_256_128_with_2LDPC_ver2   108  // p1: LDPC (160,128), p2: LDPC(224,128) 
#define Turbo_320_128_with_2LDPC        109  // p1: LDPC (224,128), p2: LDPC(224,128) 
#define RandomCode_with_LDPC            110  
#define Turbo_256_128_RC_LDPC_192_128   111 
#define Polar_256_128_Astar             112
#define Polar_with_Astar_middle_index   113

using std::vector;

void BCH_Generator_Matrix(size_t G1[], MATRIX<__int8> &G);

void Extended_BCH_Generator_Matrix(size_t G1[], MATRIX<__int8> &G);

void MatrixForm_to_Generator(MATRIX<__int8> &G);

//void Show_All_Codeword(MATRIX<__int8> G);

void ReadFile_GeneratorMatrix(string name, MATRIX<__int8> &G);

void CreateRandomCode(size_t message_length, size_t codeword_length);

void CreateRandomCode(size_t row_start_position, size_t row_end_position, double probability, MATRIX<__int8> &H);

void FillRandomBits(size_t row_start_position, MATRIX<__int8> &H);

void GenerateParityMatrix(MATRIX<__int8> &G);

void Systematic_Linear_Block_Code_Encoder(
	MATRIX<__int8> &G,
	vector<__int8> &message_seq,
	vector<__int8> &codeword_seq
);

void HammingCode_Generator_Matrix(MATRIX<__int8> &G, MATRIX<__int8> &H, vector<__int8> &generator_polynomial, __int8 codelength);

void RM_16_11_Generator_Matrix(MATRIX<__int8> &G);

void RM_32_6_Generator_Matrix(MATRIX<__int8> &G);

void RM_order3_Generator_Matrix(int m, MATRIX<__int8> &G, MATRIX<__int8>& G_);

// Some tools
double combination(double n, double k);
vector<__int8> Vector_Mulitiplication(vector<__int8> vec1, vector<__int8> vec2);
vector<__int8> Vector_Mulitiplication(vector<__int8> vec1, vector<__int8> vec2, vector<__int8> vec3);

vector<__int8> Return_S_order1(int i1);
vector<__int8> Return_E_order1(int i1, int m);
vector<__int8> Return_Sc_order1(vector<__int8> E);
vector<__int8> DefineSets_order1(int i1, int m);

vector<__int8> Return_S_order2(int i1, int i2);
vector<__int8> Return_E_order2(int i1, int i2, int m);
vector<__int8> Return_Sc_order2(vector<__int8> E);
vector<__int8> DefineSets_order2(int i1, int i2, int m);

vector<__int8> Return_S_order3(int i1, int i2, int i3);
vector<__int8> Return_E_order3(int i1, int i2, int i3, int m);
vector<__int8> Return_Sc_order3(vector<__int8> E);
vector<__int8> DefineSets_order3(int i1, int i2, int i3, int m);

// Short LDPC
void LDPC_Code_16_12(MATRIX<__int8> &H);
void Convert_SystematicG_to_H(MATRIX<__int8> &H, MATRIX<__int8> &G);
inline void Matrix_Mul(MATRIX<__int8> & M1, MATRIX<__int8> & M2, MATRIX<__int8>& res);
inline void GH_test(MATRIX<__int8> & M1, MATRIX<__int8> & M2);
inline void repeatition(MATRIX<__int8> & G);
inline void convert_H_G(MATRIX<__int8> & H, MATRIX<__int8> & G);
inline void exchang_column(MATRIX<__int8> & G, int level);