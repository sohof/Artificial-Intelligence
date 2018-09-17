#include "Matrix.h"

class HmmComputation
{
public:
	HmmComputation(int nrOfEmissions, int nrOfStates, Matrix * seq);
	//seq is the obs.sequence, a "1 x nrOfEmissions" matrix
	HmmComputation();
	void ComputeMostLikelyHiddenStateSequence(Matrix * b, Matrix * pi, Matrix * a);
	double ComputeSequenceProb(Matrix *a, Matrix *b, Matrix *pi);

	double ComputeLogSequenceProb();
	void CreateInitialDeltaVector(Matrix * pi, Matrix * b, int seqValue);
	void CreateDeltaVector(Matrix * a, Matrix * b, int seqValue, int row);

	void CreateInitialAlphaPassVector(Matrix * pi, Matrix * b, int alphaColumn ,int seqCol);
	void CreateAlphaPassVector(Matrix * b, Matrix * a, int alphaColum, int seqCol);

	void CreateInitialBetaVector();
	void CreateBetaPassVector(Matrix * b, Matrix * a, int betaColumn, int seqColum);

	void AlphaPassScaled(Matrix * a, Matrix * b, Matrix * pi);
	void BetaPassScaled(Matrix * a, Matrix * b, Matrix * pi);

	void BaumWelch(Matrix * a, Matrix * b, Matrix * pi);
	void GammaPass(Matrix * a, Matrix * b);

	void ReEstimateModel(Matrix * a, Matrix * b, Matrix * pi);

	int nrOfEmissions;
	int nrOfStates;
	Matrix *Seq;
	Matrix *Alpha;
	Matrix *Beta;
	Matrix *Scales; // vector of scaling factors
	Matrix *Gamma;
	Matrix ** DiGamma;
	Matrix *Delta;
	Matrix *DeltaIdx;
};