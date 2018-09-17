#include "Matrix.h"

class HmmComputation
{
public:
	HmmComputation(int nrOfEmissions, Matrix *seq);
	HmmComputation();
	void ComputeMostLikelyHiddenStateSequence(Matrix * b, Matrix * pi, Matrix * a);
	void CreateInitialDeltaVector(Matrix * pi, Matrix * b, int seqValue);
	void CreateDeltaVector(Matrix * a, Matrix * b, int seqValue, int row);
	double ComputeSequenceProb(Matrix *b, Matrix *pi, Matrix *a);
	void CreateInitialAlphaPassVector(Matrix * pi, Matrix * b, int seqCol);
	void CreateAlphaPassVector(Matrix * b, Matrix * a, int row, int seqCol);
	int NrOfEmissions;
	Matrix *Seq;
	Matrix *Alpha;
	Matrix *Delta;
	Matrix *DeltaIdx;
};