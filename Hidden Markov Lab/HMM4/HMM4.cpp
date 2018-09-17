#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <string>
#include "HmmComputation.h"

using namespace std;

double** CreateSequenceMatrix(char *str, int * rowArg, int * colArg)
{
	const char delim[2]= " ";
	char * dimensionCol;
	dimensionCol = strtok(str,delim);

	int cols = atof(dimensionCol);
	int rows = 1; // sequence matrix has only one row
	*rowArg = rows; 
	*colArg = cols;
	double** vec = new double*[rows];
	for (int i = 0; i < rows; i++) {
		vec[i] = new double[cols];
	}
	char *lineValues;
	lineValues = strtok(NULL, delim);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (lineValues != NULL) {
				vec[i][j] = atof(lineValues);
				lineValues = strtok(NULL, delim);
			}
		}
	}
	return vec;
}

double** Create2DimensionalArray(char *str, int * rowArg, int * colArg) {

	const char delim[2]= " ";
	char * dimensionRow;
	char * dimensionCol;
	dimensionRow = strtok(str,delim);
	dimensionCol = strtok(NULL,delim);

	int rows = atof(dimensionRow);
	int cols = atof(dimensionCol);

	*rowArg = rows;
	*colArg = cols;
	double** vec = new double*[rows];
	for (int i = 0; i < rows; i++) {
		vec[i] = new double[cols];
	}

	char *lineValues;
	lineValues = strtok(NULL, delim);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (lineValues != NULL) {
				vec[i][j] = atof(lineValues);
				lineValues = strtok(NULL, delim);
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

		double** values;
		switch (i)
		{
		case 0:
			values = Create2DimensionalArray(&lineInput[0], &rows, &cols);
			A = new Matrix(rows, cols);
			A->SetValues(values);
			free(values);
			break;
		case 1:
			values = Create2DimensionalArray(&lineInput[0], &rows, &cols);
			B = new Matrix(rows, cols);
			B->SetValues(values);
			free(values);
			break;
		case 2:
			values = Create2DimensionalArray(&lineInput[0], &rows, &cols);
			pi = new Matrix(rows, cols);
			pi->SetValues(values);
			free(values);
			break;
		case 3:
			values = CreateSequenceMatrix(&lineInput[0], &rows, &cols);
			seq = new Matrix(rows, cols);
			seq->SetValues(values);
			hmmComputation = new HmmComputation(cols, A->rows, seq);//cols here is the nr of Observation
			free(values);
			break;
		default:
			throw invalid_argument("Failed to determine which HMM matrix to initialize.");
			break;
		}
	}

	/*	Kod för HMM1    MÅSTE ÄNDRA FOR LOOPEN TILL ATT BARA GÅ TILL 3 !!
	Matrix * stateDist = pi->MultiplyMat(*A);
	Matrix* observationDist = stateDist->MultiplyMat(*B);
	observationDist->PrintMatrix(false);*/

	//hmmComputation->ComputeSequenceProb(A,B,pi);
	//cout << "----------------------------------------" << endl;
	hmmComputation->BaumWelch(A,B,pi);

	//Clean up
	delete(A);
	delete(B);
	delete(pi);
	delete(seq);
	delete(hmmComputation);
	return 0;
}
