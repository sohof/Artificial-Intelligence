class Matrix {

public:
	Matrix(int rows, int columns);
	Matrix * MultiplyMat(Matrix b);
	void SetValues(double **values);
	double SumAtRow(int row);
	void PrintMatrix(bool skip);
	double **matrix;
	int rows;
	int columns;

};