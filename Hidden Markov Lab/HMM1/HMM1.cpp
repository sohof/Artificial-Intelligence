#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <string>
#include "Matrix.h"
using namespace std;

double** Create2DimensionalArray(const char *line, int r, int c){
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

	string lineInput;
	for(int i = 0; i < 3; i++)
	{
		getline(cin, lineInput);

		int rows = (int)lineInput[0] - 48;
		int cols = (int)lineInput[2] - 48;
		double** values;
		switch(i)
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
			default:
			throw invalid_argument("Failed to determine which HMM matrix to initialize.");
			break;
		}
	}

	Matrix* stateDist = pi->MultiplyMat(*A);
	Matrix* observationDist = stateDist->MultiplyMat(*B);
	observationDist->PrintMatrix();
return 0;
}
