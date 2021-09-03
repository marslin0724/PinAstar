#pragma once

#include <cstdlib>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <cmath>
#include "AWGN.h"
#include "DefineClass.h"
#include "AStarDecode.h"
#include "DefineParameter.h"
#include "MatrixOperation.h"
#include "LinearBlockCodes.h"
#include "FileOperation.h"
#include "LDPC_OPER.h"
#include "TURBO_OPER.h"

bool ML_LowerBound(
	vector<double> &received_signal,
	vector<double> &transmitted_signal,
	vector<__int8> &estimated_codeword);

double vectorDistance(vector<double> Input1, vector <double> Input2) {
	double ret = 0.0;
	vector<double>::iterator
		first1 = Input1.begin(),
		last1 = Input1.end(),
		first2 = Input2.begin();
	while (first1 != last1) {
		double dist = (*first1++) - (*first2++);
		ret += dist * dist;
	}
	return ret > 0.0 ? sqrt(ret) : 0.0;
}