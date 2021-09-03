#pragma once

#include <iostream>
#include "DefineClass.h"

void GJ_Elimination(MATRIX<__int8> &input);

void Sort_Matrix_Col(
	MATRIX<__int8>		&input_matrix,
	vector<size_t>		&permutation_seq, 
	MATRIX<__int8>		&output_matrix
);

void Sort_Matrix_Row(
	MATRIX<__int16>		&input_matrix,
	vector<size_t>		&permutation_seq,
	MATRIX<__int16>		&output_matrix
);
void H_G_convertor_G_P_I(MATRIX<__int8> & H, vector<int> &permutation_seq, MATRIX<__int8> & G);
void Col_Exchange(MATRIX<__int8> &Matrix, int col1, int col2);
void Row_Exchange(MATRIX<__int8> &Matrix, int row1, int row2);
void Row_Addition(MATRIX<__int8> &Matrix, int row1, int row2);
void Col_Addition(MATRIX<__int8> &Matrix, int col1, int col2);
void Sort_Matrix_Col_Forward(MATRIX<__int8> &Matrix, vector<int> &permutation_seq);
void G_H_Multiple_Test(MATRIX<__int8> G, MATRIX<__int8> H);
void Transpose_Matrix(MATRIX<__int8> &H, MATRIX<__int8> &Transpose_H);
