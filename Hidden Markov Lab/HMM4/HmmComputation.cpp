#include "HmmComputation.h"
//#include <stdexcept>
//#include <algorithm>
#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;

HmmComputation::HmmComputation(int nrEmissions, int nrStates, Matrix * seq)
{
	nrOfEmissions = nrEmissions;
	nrOfStates = nrStates;
	Seq = seq;
}

HmmComputation::HmmComputation()
{
	nrOfEmissions = 0;
	nrOfStates = 0;
	Seq = NULL;

}

void HmmComputation::ComputeMostLikelyHiddenStateSequence(Matrix *b, Matrix *pi, Matrix*a)
{
	int nrOfStates = b->rows;
	int *result = new int[nrOfEmissions];
	Delta = new Matrix(nrOfStates, nrOfEmissions); //b.rows = nr of hidden states
	DeltaIdx = new Matrix(nrOfStates, nrOfEmissions);
	//The first column in deltaidx should be zero since it is time zero and there are 
	//no prev. states.

	CreateInitialDeltaVector(pi, b, Seq->matrix[0][0]);
	// after init.delta vector we start computing second deltavector for second observ.value
	// second column in Delta matrix
	for (int i = 1; i < nrOfEmissions ; i++)
	{ 
		CreateDeltaVector(a, b, Seq->matrix[0][i], i);
	}

	double delta_T =0; 	int argMax =-1; 
	for (int i = 0; i<nrOfStates; i++){

		double delta_tmp_T = Delta->matrix[i][nrOfEmissions -1];
		if (delta_tmp_T >= delta_T) {
			delta_T = delta_tmp_T;
			argMax = i;
		}
	}
	result[0]=argMax;
	int col_idx = nrOfEmissions -1; 
	for (int i =1; i<nrOfEmissions; ++i){
		int row = result[i-1];
		result[i] = DeltaIdx->matrix[row][col_idx];
		col_idx--;
	}
	for (int i=nrOfEmissions-1; i>=0; --i){
		cout << result[i] << " ";
	}
	cout << endl;
}

void HmmComputation::CreateInitialDeltaVector(Matrix *pi, Matrix *b, int seqValue)
{
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

		for (int j = 0; j < nrOfStates; j++) //figuring out which state j maximizes delta_t(i)
		{
			tmpDelta = a->matrix[j][i] * Delta->matrix[j][emissionNr - 1] * b->matrix[i][seqValue];
			if(tmpDelta > deltaProb) 
			{
			 	deltaProb = tmpDelta; 
			 	stateThatMaximized = j; 
			}
		}
		Delta->matrix[i][emissionNr] = deltaProb;
		DeltaIdx->matrix[i][emissionNr] = stateThatMaximized;

	}
}

double HmmComputation::ComputeSequenceProb(Matrix *a, Matrix *b, Matrix *pi)
{
	Alpha = new Matrix(nrOfStates, nrOfEmissions);
	Scales = new Matrix(1,nrOfEmissions);

	AlphaPassScaled(a,b,pi);

	double res =1;	
	for (int t=0; t<nrOfEmissions; ++t){

		res *= Scales->matrix[0][t];
	}

	return (1/res); //Alpha->SumAtColumn(nrOfEmissions - 1);
}

void HmmComputation::BaumWelch(Matrix * a, Matrix * b, Matrix * pi) 
{
	int iters = 0;
	int maxIters  = 10000;
	double oldLogProb = -10000000;
	double logProb = -1000000;

	// Set up all the matrices needed
	Alpha = new Matrix(nrOfStates,nrOfEmissions);
	Scales = new Matrix(1,nrOfEmissions);
	Beta = new Matrix(nrOfStates, nrOfEmissions);

	DiGamma = new Matrix*[nrOfEmissions-1]; //Digamma only goes to T-2. I.e 
	for (int t=0; t<(nrOfEmissions-1); ++t) {
			DiGamma[t] = new Matrix(nrOfStates,nrOfStates);
		} 
	Gamma = new Matrix(nrOfStates,nrOfEmissions);
	
	do 
	{
		if (logProb > oldLogProb) {
			oldLogProb =logProb;
		}
		else 
			break;
		// Alfa pass must be performed before beta pass, the Scaled Matrix must to be 
		AlphaPassScaled(a, b, pi);
		BetaPassScaled(a, b, pi);
		GammaPass(a, b);
		ReEstimateModel(a,b,pi);
		logProb = ComputeLogSequenceProb();

		iters++;

	} while ((iters < maxIters) && (logProb > oldLogProb));

	a->PrintMatrix(true);
	b->PrintMatrix(true);
	cout << "Nr of iterations " << iters << endl;
	cout << "logprob value =  " << logProb << endl;
	cout << "Oldlogprob value =  " << oldLogProb << endl;

	delete(Alpha);
	delete(Beta);
	delete(Scales);
}

void HmmComputation::AlphaPassScaled(Matrix * a, Matrix * b, Matrix * pi)
{
	// compute the initial vector alpha_0(i) for time=0 and the scale to use
	int o_time_0 = Seq->matrix[0][0]; // observation at time t=0.
	double c_0=0; //scale factor for time t=0
	for (int i = 0; i < nrOfStates; i++)
	{
		Alpha->matrix[i][0] = pi->matrix[0][i] * b->matrix[i][o_time_0];
		c_0 += Alpha->matrix[i][0];
	}

	Scales->matrix[0][0]= 1.0 / c_0;  // Save the c_0 scale calculated for time t=0.

	// scale the initial vector alpha_0(i)
	for (int i = 0; i < nrOfStates; i++)
	{
		Alpha->matrix[i][0] = Alpha->matrix[i][0] * Scales->matrix[0][0];
	}

	// Compute all the other vectors alpha_t(i) for time t=1 to T-1. nrOfEmmission = time T
	for (int t = 1; t < nrOfEmissions; t++) {
		int o_time_t = Seq->matrix[0][t]; // observation at time t
		double c_t= 0; // scale factor time t 

		for (int i = 0; i < nrOfStates; i++)
		{
			Alpha->matrix[i][t]=0; 
			for (int j = 0; j < nrOfStates; j++)
			{
				Alpha->matrix[i][t] += Alpha->matrix[j][t-1] * a->matrix[j][i];
			}

			Alpha->matrix[i][t] = Alpha->matrix[i][t] * b->matrix[i][o_time_t];

			c_t+= Alpha->matrix[i][t];
		}

		Scales->matrix[0][t]= 1.0 / c_t; // Save the c_t scale calculated for time t.
		// scale the current alpha_t(i) vector 
		for (int i = 0; i < nrOfStates; i++)
		{
			Alpha->matrix[i][t] = Alpha->matrix[i][t] * Scales->matrix[0][t];
		}
	}
}
void HmmComputation::BetaPassScaled(Matrix * a, Matrix * b, Matrix * pi)
{
	// Create BetaPass Initial vector. Let B_T-1(i) scaled by C_T-1

	for(int i =0; i<nrOfStates; ++i) {
		Beta->matrix[i][nrOfEmissions - 1] = Scales->matrix[0][nrOfEmissions - 1];
	}

	//Beta pass works backward from 2nd last column. At each beta vector column/time i  we look
	// at the observation for colum i+1, i.e time t+1
	//	for(int t = nrOfEmissions - 2; t >= 0;  --t) {  

	for(int t = nrOfEmissions - 2; t >= 0;  t--) {  

    	int o_tplus1 = Seq->matrix[0][t + 1];
    	for(int i=0; i < nrOfStates; ++i) {
			
			Beta->matrix[i][t] = 0;
			for (int j=0; j < nrOfStates; ++j)	{

				Beta->matrix[i][t]+= a->matrix[i][j]* b->matrix[j][o_tplus1]* Beta->matrix[j][t+1];
			}
			//scale Beta_t(i) with same scale factor as Alfa_t(i)
			Beta->matrix[i][t] = Scales->matrix[0][t] * Beta->matrix[i][t];
		}
	}
}

void HmmComputation::CreateInitialAlphaPassVector(Matrix *b, Matrix *pi, int alphaCol, int seqCol)
{
	for (int j = 0; j < nrOfStates; j++)
	{
		Alpha->matrix[j][alphaCol] = pi->matrix[0][j] * b->matrix[j][seqCol];
	}
}

void HmmComputation::CreateAlphaPassVector(Matrix *a, Matrix *b, int alphaCol, int seqCol)
{
	for (int i = 0; i < nrOfStates; i++)
	{
		double sum = 0;
		for (int j = 0; j < nrOfStates; j++)
		{
			sum += a->matrix[j][i] * Alpha->matrix[j][alphaCol - 1];
		}
		Alpha->matrix[i][alphaCol] = sum * b->matrix[i][seqCol];
	}
}

void HmmComputation::CreateInitialBetaVector()
{
	int lastColum = nrOfEmissions - 1;
	for(int i =0; i<nrOfStates; ++i){

		Beta->matrix[i][lastColum] = 1;
	}
}

void HmmComputation::CreateBetaPassVector(Matrix * a, Matrix * b, int betaCol, int seqCol)
{

	for(int i=0; i < nrOfStates; ++i) {
		double sum =0;

		for (int j=0; j < nrOfStates; ++j){

			sum+= Beta->matrix[j][betaCol+1] * b->matrix[j][seqCol] * a->matrix[i][j];
		}
	Beta->matrix[i][betaCol] = sum;
	}
}

void HmmComputation::GammaPass(Matrix * a, Matrix * b)
{
	int t_minus1 = nrOfEmissions-1;
	for (int t = 0; t < t_minus1; t++) { // the for loop goes from 0 to T-2

		double denom =0; 

		//These two for loops, seems to be just a way to calc. P(O|lambda). See Rabiner p.264
		for (int i=0; i<nrOfStates; i++) {

			for (int j =0; j<nrOfStates; j++) {

			int o_tplus1 = Seq->matrix[0][t + 1];
			denom += Alpha->matrix[i][t] * a->matrix[i][j] * b->matrix[j][o_tplus1] * Beta->matrix[j][t+1];
			}
		}

		for (int i=0; i < nrOfStates; i++){

			Gamma->matrix[i][t] = 0; 

			for (int j = 0; j < nrOfStates; j++) {

				int o_tplus1 = Seq->matrix[0][t + 1];
			    DiGamma[t]->matrix[i][j]= (Alpha->matrix[i][t]*a->matrix[i][j]*b->matrix[j][o_tplus1]*Beta->matrix[j][t+1]) /denom;
				Gamma->matrix[i][t] += DiGamma[t]->matrix[i][j]; 
			}
		}
	}

	//Special Case for Gamma_T-1(i)
	double denom=0;
	for (int i =0; i < nrOfStates; i++){
		denom += Alpha->matrix[i][nrOfEmissions-1];
	}
	for (int i =0; i < nrOfStates; i++){
		Gamma->matrix[i][nrOfEmissions-1] = Alpha->matrix[i][nrOfEmissions-1] / denom;
	}
}

void HmmComputation::ReEstimateModel(Matrix * a, Matrix * b, Matrix * pi){

	//Re-estimate pi
	//cout << " nrOfStates before entering Re-est Model = "<< nrOfStates << endl;

	for(int i =0; i<nrOfStates; i++) {
		pi->matrix[0][i]= Gamma->matrix[i][0];
	}
	// cout << " Finished Re-estimating pi  " << endl;
	// pi->PrintMatrix(true);

	//Re-estimate A

	int t_minus1 = nrOfEmissions-1;
	for(int i =0; i < nrOfStates; i++) {

		for (int j=0; j < nrOfStates; j++) {

			double numer =0; 
			double denom =0;

			for (int t=0; t < nrOfEmissions-1; t++) {

				numer += DiGamma[t]->matrix[i][j]; 
				denom += Gamma->matrix[i][t];
			}
			a->matrix[i][j]= numer / denom;
		}
	}

	//cout << " Finished Re-estimating A " << endl;
	//a->PrintMatrix(true);



	//Re-estimate B.
	//Important 2nd for loop has to run to the NR of obser.symbols which is not always 
	// the same as the nr of states!! 
	int nrOfObservationsSymbols = b->columns;
	for(int i =0; i<nrOfStates; i++) {

		for (int j=0; j<nrOfObservationsSymbols; j++) {

			double numer =0; 
			double denom =0;
			for (int t=0; t< nrOfEmissions-1 ; t++) {  // looping form t=0 to T-2. Last one T-1 is not included. 

				if ((Seq->matrix[0][t]) == j){ // Our indicator function! 
					numer += Gamma->matrix[i][t]; 
				}
				denom += Gamma->matrix[i][t];
			}
			b->matrix[i][j]= numer / denom;
		}
	}
	//cout << " Finished Re-estimating B " << endl;
	//b->PrintMatrix(true);
}


double HmmComputation::ComputeLogSequenceProb()
{
	double logProb=0;
	for (int t=0; t<nrOfEmissions; t++) {

		logProb += log(Scales->matrix[0][t]);
	}
	return -1*logProb;
}



