#include "HmmComputation.h"
#include <stdexcept>
using namespace std;

HmmComputation::HmmComputation(int nrOfEmissions, Matrix * seq)
{
	NrOfEmissions = nrOfEmissions;
	Seq = seq;
}

HmmComputation::HmmComputation()
{
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