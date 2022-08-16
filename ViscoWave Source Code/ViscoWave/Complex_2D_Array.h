#pragma once
#include "ap.h"

#ifndef Complex_2D_Array_H
#define Complex_2D_Array_H

class Complex_2D_Array {
private:
	int row, column;

public:
	alglib::complex_2d_array array_2d;
	Complex_2D_Array(int m, int n);
	void Reset();
	~Complex_2D_Array();
};

#endif
