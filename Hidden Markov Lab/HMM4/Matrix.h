class Matrix {

public:
	Matrix(int rows, int columns);
	Matrix * MultiplyMat(Matrix b);
	void SetValues(double **values);
	double SumAtRow(int row);
	double SumAtColumn(int column);
	void PrintMatrix(bool skipDimensions);
	void PrintOutput();
	double **matrix;
	int rows;
	int columns;

};