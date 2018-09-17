#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <string>
#include "HmmComputation.h"

using namespace std;

double** Create2DimensionalArray(const char *line, int r, int c) {
	double** vec = new double*[r];
	for (int i = 0; i < r; i++) {
		vec[i] = new double[c];
	}

	char *lineValues;
	lineValues = strtok(strdup(line), " ");

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			if (lineValues != NULL) {
				vec[i][j] = atof(lineValues);
				lineValues = strtok(NULL, " ");
			}
		}
	}

	return vec;
}

int main(int argv, char *argc[]) {

	Matrix *A;
	Matrix *B;
	Matrix *pi;
	Matrix *seq;
	HmmComputation *hmmComputation;

	string lineInput;
	for (int i = 0; i < 4; i++)
	{
		getline(cin, lineInput);

		int rows = 0; 
		int cols = 0;
		
		if (i < 3) {
			rows = (int)lineInput[0] - 48;
			cols = (int)lineInput[2] - 48;
		}
		else if (i == 3)
		{
			rows = 1;
			cols = (int)lineInput[0] - 48;
		}

		double** values;
		switch (i)
		{
		case 0:
			A = new Matrix(rows, cols);
			values = Create2DimensionalArray(lineInput.substr(4).c_str(), rows, cols);
			A->SetValues(values);
			free(values);
			break;
		case 1:
			B = new Matrix(rows, cols);
			values = Create2DimensionalArray(lineInput.substr(4).c_str(), rows, cols);
			B->SetValues(values);
			free(values);
			break;
		case 2:
			pi = new Matrix(rows, cols);
			values = Create2DimensionalArray(lineInput.substr(4).c_str(), rows, cols);
			pi->SetValues(values);
			free(values);
			break;
		case 3:
			seq = new Matrix(1, cols);
			hmmComputation = new HmmComputation(cols, seq);
			values = Create2DimensionalArray(lineInput.substr(2).c_str(), 1, cols);
			seq->SetValues(values);
			free(values);
			break;
		default:
			throw invalid_argument("Failed to determine which HMM matrix to initialize.");
			break;
		}
	}

	double sequenceProbability = hmmComputation->ComputeSequenceProb(B, pi, A);
	cout << sequenceProbability << endl;

	return 0;
}
