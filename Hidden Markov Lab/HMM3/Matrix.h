class Matrix {

public:
	Matrix(int rows, int columns);
	Matrix * MultiplyMat(Matrix b);
	void SetValues(double **values);
	double SumAtRow(int row);
	void PrintMatrix(bool skipDimensions);
	double **matrix;
	int rows;
	int columns;

};