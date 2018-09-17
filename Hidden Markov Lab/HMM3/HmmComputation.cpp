#include "HmmComputation.h"
#include <stdexcept>
#include <algorithm>
#include <stdio.h>
#include <iostream>
using namespace std;

HmmComputation::HmmComputation(int nrOfEmissions, Matrix * seq)
{
	NrOfEmissions = nrOfEmissions;
	Seq = seq;
}

HmmComputation::HmmComputation()
{
}

void HmmComputation::ComputeMostLikelyHiddenStateSequence(Matrix *b, Matrix *pi, Matrix*a)
{
	int nrOfStates = b->rows;
	//cout << "NrOfEmissions = " << NrOfEmissions << " brows = " << b->rows; 
	int *result = new int[NrOfEmissions];
	Delta = new Matrix(nrOfStates, NrOfEmissions); //b.rows = nr of hidden states
	DeltaIdx = new Matrix(nrOfStates, NrOfEmissions);
	//The first column in deltaidx should be zero since it is time zero and there are 
	//no prev. states.

	CreateInitialDeltaVector(pi, b, Seq->matrix[0][0]);
	// after init.delta vector we start computing second deltavector for second observ.value
	// second column in Delta matrix
	for (int i = 1; i < NrOfEmissions ; i++)
	{ 
		CreateDeltaVector(a, b, Seq->matrix[0][i], i);
	}
	/*
	cout << " Printing delta matrix " << endl;
	Delta->PrintMatrix(true);
	cout << " Finished Printing delta matrix " << endl;
	cout << " Printing DeltaIdx matrix " << endl;
	DeltaIdx->PrintMatrix(true); */

	double delta_T =0; 	int argMax =-1; 
	for (int i = 0; i<nrOfStates; i++){

		double delta_tmp_T = Delta->matrix[i][NrOfEmissions -1];
		if (delta_tmp_T >= delta_T) {
			delta_T = delta_tmp_T;
			argMax = i;
		}
	}
	result[0]=argMax;
	int col_idx = NrOfEmissions -1; 
	for (int i =1; i<NrOfEmissions; ++i){
		int row = result[i-1];
		result[i] = DeltaIdx->matrix[row][col_idx];
		col_idx--;
	}
	for (int i=NrOfEmissions-1; i>=0; --i){
		cout << result[i] << " ";
	}
	cout << endl;
}

void HmmComputation::CreateInitialDeltaVector(Matrix *pi, Matrix *b, int seqValue)
{
	//cout << " Printing seqValue, i.e first observation value = " << seqValue <<endl;
	int nrOfStates = b->rows;
	for (int i = 0; i < nrOfStates; i++)
	{
		double probability = b->matrix[i][seqValue] * pi->matrix[0][i];
		//cout << " Probabilty = " << probability <<endl;
		Delta->matrix[i][0] = probability;
	}
}
// Create the deltavector and the delta idx vector for the given observation 
void HmmComputation::CreateDeltaVector(Matrix *a, Matrix *b, int seqValue, int emissionNr)
{

	int nrOfStates = b->rows;

	for (int i = 0; i < nrOfStates; i++) // creating delta_t(i) for all states
	{
		double deltaProb =0; double tmpDelta =0;
		int stateThatMaximized = -1;
		//cout << "For delta_1:State " << i << endl;   

		for (int j = 0; j < nrOfStates; j++) //figuring out which state j maximizes delta_t(i)
		{
			tmpDelta = a->matrix[j][i] * Delta->matrix[j][emissionNr - 1] * b->matrix[i][seqValue];
			if(tmpDelta > deltaProb) 
			{
			 	deltaProb = tmpDelta; 
			 	stateThatMaximized = j; 
			}
		//	cout << "aji value = " <<a->matrix[j][i]<<". d_t-1("<<i<< ") = "<<Delta->matrix[j][emissionNr - 1]
		//	 <<". b_" <<i<<"("<<seqValue<<") = "<<b->matrix[i][seqValue] <<". Probability = "<<deltaProb<<endl;
		}
		//cout << "------------------------------"<<endl;
		Delta->matrix[i][emissionNr] = deltaProb;
		DeltaIdx->matrix[i][emissionNr] = stateThatMaximized;

	}
}

double HmmComputation::ComputeSequenceProb(Matrix *b, Matrix *pi, Matrix *a)
{
	Alpha = new Matrix(NrOfEmissions, pi->columns);
	int seqValue = Seq->matrix[0][0];
	CreateInitialAlphaPassVector(pi, b, seqValue);

	for (int i = 1; i < NrOfEmissions; i++)
	{
		seqValue = Seq->matrix[0][i];
		CreateAlphaPassVector(b, a, i, seqValue);
	}
	return Alpha->SumAtRow(NrOfEmissions - 1);
}

void HmmComputation::CreateInitialAlphaPassVector(Matrix *pi, Matrix *b, int seqCol)
{
	for (int j = 0; j < pi->columns; j++)
	{
		Alpha->matrix[0][j] = pi->matrix[0][j] * b->matrix[j][seqCol];
	}
}

void HmmComputation::CreateAlphaPassVector(Matrix *b, Matrix *a, int alphaPassRow, int seqCol)
{
	for (int j = 0; j < a->columns; j++)
	{
		double sum = 0;
		for (int k = 0; k < a->rows; k++)
		{
			sum += Alpha->matrix[alphaPassRow - 1][k] * a->matrix[k][j];
		}
		Alpha->matrix[alphaPassRow][j] = sum * b->matrix[j][seqCol];
	}
}