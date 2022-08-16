#include "Complex_2D_Array.h"

Complex_2D_Array::Complex_2D_Array(int m, int n)
{
	row = m;
	column = n;

	array_2d.setlength(row, column);
}

void Complex_2D_Array::Reset() {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < column; j++) {
			array_2d(i, j).x = 0;
			array_2d(i, j).y = 0;
		}
	}
}

Complex_2D_Array::~Complex_2D_Array()
{
}
