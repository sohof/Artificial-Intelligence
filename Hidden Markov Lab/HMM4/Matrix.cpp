#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <stdio.h>
#include "Matrix.h"
using namespace std;

Matrix::Matrix(int r, int c)
{
	rows =r;
	columns =c;

	matrix = new double*[rows];

	for (int i=0; i<rows; i++)
	{
		matrix[i]= new double[columns];
	}

	for (int i=0; i<rows; i++ ){
		for (int j=0; j<columns; j++)
			matrix[i][j] = 0.0;
	}
}
void Matrix::SetValues(double** values){

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < columns; j++){
			matrix[i][j] = values[i][j];
		}
	}
}

double Matrix::SumAtRow(int row)
{
	double sum = 0;
	for (int j = 0; j < columns; j++)
	{
		sum += matrix[row][j];
	}
	return sum;
}
double Matrix::SumAtColumn(int column)
{
	double sum = 0;
	for (int j = 0; j < rows; j++)
	{
		sum += matrix[j][column];
	}
	return sum;
}

Matrix * Matrix::MultiplyMat(Matrix b) {
	if (columns!= b.rows) 
		throw std::invalid_argument("Matrix dimensions incompatible"); 

	Matrix * result = new Matrix(rows,b.columns);

	for (int i=0; i<rows; i++){
		for (int j=0; j<b.columns; j++) {

			double sum = 0;
			for (int k =0; k< b.rows; k++) {
				sum += matrix[i][k]* b.matrix[k][j];
			}
			result->matrix[i][j] += sum;
		}
	}

	return result;
}
void Matrix::PrintMatrix(bool skipDimensions = false){
	if(!skipDimensions)
		std::cout << rows << " " << columns << " ";

	for (int i=0; i<rows; i++ )
	{
		for (int j=0; j<columns; j++) 
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
    cout << endl;
}
void Matrix::PrintOutput(){

	cout << rows << " " << columns << " ";

	for (int i=0; i<rows; i++ )
	{
		for (int j=0; j<columns; j++) 
		{
			cout << matrix[i][j] << " ";
		}
	}
    cout << endl;
}





