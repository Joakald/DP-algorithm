#pragma once

#include <vector>
#include <assert.h>

template<class T>
class Matrix
{
public:
	Matrix(int in_rows, int in_columns)
		:
		rows(in_rows),
		columns(in_columns)
	{
		elements = std::vector<T>(rows * columns, T);
	}
	Matrix(int in_rows, int in_columns, T in_t)
		:
		rows(in_rows),
		columns(in_columns)
	{
		elements = std::vector<T>(rows * columns, in_t);
	}
public:
	T Get(int row, int column)
	{
		assert(row >= 0);
		assert(row < rows);
		assert(column >= 0);
		assert(column < columns);

		return elements[row * columns + column];
	}
	void Set(int row, int column, T value)
	{
		assert(row >= 0);
		assert(row < rows);
		assert(column >= 0);
		assert(column < columns);

		elements[row * columns + column] = value;
	}
public:
	int rows;
	int columns;
	std::vector<T> elements;
};
typedef Matrix<int> MatrixInt;
typedef Matrix<double> MatrixDouble;